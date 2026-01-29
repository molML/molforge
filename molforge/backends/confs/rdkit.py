"""
RDKit ETKDG conformer generation backend.

Uses RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry) method
for conformer generation. Stores molecules in memory for efficient processing.
"""

from typing import Dict, Any, Iterator, Tuple
import time
import pickle
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem

from .base import ConformerBackend
from ...utils.actortools.multiprocess import calculate_chunk_params


def _generate_conformer_chunk(
    chunk: list[Tuple[str, str]],
    max_confs: int,
    random_seed: int,
    rms_threshold: float,
    use_random_coords: bool,
    use_uff: bool,
    max_iterations: int
) -> list[Tuple[str, Dict[str, Any]]]:
    """
    Generate conformers for a chunk of SMILES.

    Module-level function for efficient parallel processing with chunking.
    Processes multiple molecules in a single worker to reduce overhead.

    Args:
        chunk: List of (smiles, name) tuples
        max_confs: Maximum number of conformers to generate
        random_seed: Random seed for reproducibility
        rms_threshold: RMS threshold for conformer pruning
        use_random_coords: Whether to use random coordinates
        use_uff: Whether to use UFF optimization
        max_iterations: Maximum iterations for UFF optimization

    Returns:
        List of (name, result_dict) tuples for each molecule in chunk
    """
    results = []
    for smiles, name in chunk:
        result = _generate_single_conformer(
            smiles, name,
            max_confs, random_seed, rms_threshold,
            use_random_coords, use_uff, max_iterations
        )
        results.append(result)
    return results


def _generate_single_conformer(
    smiles: str,
    name: str,
    max_confs: int,
    random_seed: int,
    rms_threshold: float,
    use_random_coords: bool,
    use_uff: bool,
    max_iterations: int
) -> Tuple[str, Dict[str, Any]]:
    """
    Generate conformers for a single SMILES using RDKit ETKDG.

    Module-level function for efficient parallel processing.
    Does not access instance attributes, making it safe for multiprocessing.

    Args:
        smiles: SMILES string
        name: Molecule identifier
        max_confs: Maximum number of conformers to generate
        random_seed: Random seed for reproducibility
        rms_threshold: RMS threshold for conformer pruning
        use_random_coords: Whether to use random coordinates
        use_uff: Whether to use UFF optimization
        max_iterations: Maximum iterations for UFF optimization

    Returns:
        Tuple of (name, result_dict) where result_dict contains:
            - success: bool
            - n_conformers: int
            - status: str (empty for success, error description for failure)
            - molecule: Chem.Mol (only if successful)
            - rotors: int
            - elapsed_time: float
    """
    start_time = time.time()

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (name, {
                'success': False,
                'n_conformers': 0,
                'status': 'Invalid SMILES',
                'rotors': 0,
                'elapsed_time': time.time() - start_time
            })

        # Add hydrogens
        mol = Chem.AddHs(mol)
        mol.SetProp("_Name", name)

        # Count rotatable bonds
        rotors = Chem.Lipinski.NumRotatableBonds(mol)

        # Configure ETKDG parameters
        etkdg_params = AllChem.ETKDGv3()
        etkdg_params.randomSeed = random_seed
        etkdg_params.numThreads = 1  # Each worker processes one molecule
        etkdg_params.pruneRmsThresh = rms_threshold
        etkdg_params.useRandomCoords = use_random_coords

        # Generate conformers
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=max_confs,
            params=etkdg_params
        )

        if len(conf_ids) == 0:
            return (name, {
                'success': False,
                'n_conformers': 0,
                'status': 'Conformer embedding failed',
                'rotors': rotors,
                'elapsed_time': time.time() - start_time
            })

        # Optional UFF optimization
        if use_uff:
            for conf_id in conf_ids:
                try:
                    AllChem.UFFOptimizeMolecule(
                        mol,
                        confId=conf_id,
                        maxIters=max_iterations
                    )
                except Exception:
                    pass  # Continue if optimization fails

        elapsed = time.time() - start_time

        return (name, {
            'success': True,
            'n_conformers': len(conf_ids),
            'status': '',
            'molecule': mol,
            'rotors': rotors,
            'elapsed_time': round(elapsed, 2)
        })

    except Exception as e:
        return (name, {
            'success': False,
            'n_conformers': 0,
            'status': str(e),
            'rotors': 0,
            'elapsed_time': time.time() - start_time
        })


class RDKitBackend(ConformerBackend):
    """
    RDKit ETKDG conformer generation backend.

    Features:
    - In-memory storage (no file I/O)
    - ETKDG conformer generation
    - Optional UFF optimization
    - Configurable pruning threshold
    """
    __backend_name__ = 'RDK'

    # Backend metadata
    description = "Generate conformers using RDKit ETKDGv3 algorithm"
    required_params = ['SMILES_column']
    optional_params = [
        'names_column', 'max_confs', 'n_jobs', 'rms_threshold',
        'random_seed', 'verbose', 'dropna'
    ]

    def __init__(self, params, logger=None, context=None):
        """Initialize RDKit backend."""
        super().__init__(params, logger, context)

        # Storage for generated conformers (in input order)
        self._molecules: list[Chem.Mol] = []
        self._output_file: str = None  # Path to pickle file

        self.log(
            f"RDKit ETKDG initialized: max_confs={params.max_confs}, "
            f"rms={params.rms_threshold}, optimize={params.use_uff}"
        )

    @property
    def _report_file(self) -> str:
        """Path to report CSV file."""
        if self.context and hasattr(self.context, 'output_dir'):
            return str(Path(self.context.output_dir) / "RDKit_rpt.csv")
        return "RDKit_rpt.csv"

    def generate_conformers(self, smiles_list: list[str], names_list: list[str]) -> Dict[str, Dict[str, Any]]:
        """
        Generate conformers using RDKit ETKDG with parallel processing and progress reporting.

        Args:
            smiles_list: List of SMILES strings
            names_list: List of molecule identifiers

        Returns:
            Dict mapping names to generation results
        """
        # Clear previous results
        self._molecules.clear()

        # Prepare chunks (use smaller chunk size for conformer generation)
        max_chunk_size = 50  # Smaller chunks enable frequent progress updates
        n_processes, chunk_size, _ = calculate_chunk_params(len(smiles_list), max_chunk_size)

        # Create chunks of (smiles, name) pairs
        pairs = list(zip(smiles_list, names_list))
        chunks = [pairs[i:i + chunk_size] for i in range(0, len(pairs), chunk_size)]
        n_chunks = len(chunks)

        self.log(
            f"Generating conformers for {len(smiles_list):,} molecules with {n_processes} processes "
            f"({n_chunks} chunks of {chunk_size})."
        )

        # Process chunks in parallel
        results_dict = {}
        chunk_results = [None] * n_chunks
        start_time = time.time()
        completed_chunks = 0

        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            # Submit all chunk tasks
            future_to_chunk = {
                executor.submit(
                    _generate_conformer_chunk,
                    chunk,
                    self.params.max_confs,
                    self.params.random_seed,
                    self.params.rms_threshold,
                    self.params.use_random_coords,
                    self.params.use_uff,
                    self.params.max_iterations
                ): i
                for i, chunk in enumerate(chunks)
            }

            # Process as completed with progress reporting
            for future in as_completed(future_to_chunk, timeout=3600):
                chunk_idx = future_to_chunk[future]
                try:
                    chunk_results[chunk_idx] = future.result(timeout=300)
                    completed_chunks += 1

                    # Report progress after each chunk
                    self._report_progress(completed_chunks, n_chunks, start_time, chunk_size)

                except Exception as e:
                    self.log(f"Chunk {chunk_idx} failed: {e}", level='ERROR')
                    # Fill with failures for this chunk
                    chunk_results[chunk_idx] = [
                        (name, {
                            'success': False,
                            'n_conformers': 0,
                            'status': f'Chunk processing error: {str(e)}',
                            'rotors': 0,
                            'elapsed_time': 0.0
                        })
                        for _, name in chunks[chunk_idx]
                    ]
                    completed_chunks += 1

        # Flatten results into dict
        for chunk_result in chunk_results:
            if chunk_result:
                for name, result in chunk_result:
                    results_dict[name] = result

        # Store molecules and build report in input order
        report_data = []
        for smiles, name in zip(smiles_list, names_list):
            result = results_dict[name]

            # Build report entry (matches OMEGA format)
            report_data.append({
                'Molecule': smiles,
                'Title': name,
                'Rotors': result.get('rotors', 0),
                'Conformers': result['n_conformers'],
                'ElapsedTime(s)': result.get('elapsed_time', 0.0),
                'Status': result['status']
            })

            # Store successful molecules
            if result['success'] and 'molecule' in result:
                mol = result['molecule']
                # Ensure name property is set (may be lost during pickling/unpickling)
                if not mol.HasProp('_Name'):
                    mol.SetProp('_Name', name)
                self._molecules.append(mol)

        success_count = sum(1 for r in results_dict.values() if r['success'])
        self.log(f"RDKit generation complete: {success_count}/{len(smiles_list)} succeeded")

        # Write molecules to pickle file
        self._output_file = self._write_pickle_file()

        # Write report to CSV file
        self._write_report_file(report_data)

        # Clear molecules from memory to free RAM
        self._molecules.clear()

        return results_dict

    def _report_progress(self, completed_chunks: int, total_chunks: int,
                         start_time: float, chunk_size: int):
        """
        Report conformer generation progress.

        Args:
            completed_chunks: Number of chunks completed
            total_chunks: Total number of chunks
            start_time: Start time of generation
            chunk_size: Approximate molecules per chunk
        """
        elapsed = time.time() - start_time
        progress_pct = 100 * completed_chunks / total_chunks
        molecules_processed = completed_chunks * chunk_size

        if completed_chunks > 0:
            estimated_total = elapsed * total_chunks / completed_chunks
            eta = estimated_total - elapsed
            rate = molecules_processed / elapsed

            self.log(
                f"Progress: {completed_chunks}/{total_chunks} chunks "
                f"({progress_pct:.1f}%) | "
                f"Elapsed: {elapsed:.1f}s | "
                f"ETA: {eta:.1f}s | "
                f"Rate: {rate:.0f} mol/s"
            )

    def _write_pickle_file(self) -> str:
        """
        Write molecules to pickle file using RDKit binary format.

        Returns:
            Path to pickle file
        """
        # Determine output path (follow OpenEye pattern)
        if self.context and hasattr(self.context, 'output_dir'):
            output_path = Path(self.context.output_dir) / "RDKit_conformers.pkl"
        else:
            output_path = Path("RDKit_conformers.pkl")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Store molecules as (binary, name) tuples to preserve names reliably
        mol_data = []
        for mol in self._molecules:
            mol_binary = mol.ToBinary()
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
            mol_data.append((mol_binary, name))

        # Write molecule data to pickle
        with open(output_path, 'wb') as f:
            pickle.dump(mol_data, f)

        self.log(f"Wrote {len(self._molecules)} molecules to {output_path}")
        return str(output_path)

    def _write_report_file(self, report_data: list[Dict[str, Any]]) -> None:
        """
        Write conformer generation report to CSV file.

        Args:
            report_data: List of report entries (dicts)
        """
        # Create DataFrame from report data
        report_df = pd.DataFrame(report_data)

        # Write to CSV file
        Path(self._report_file).parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(self._report_file, index=False)

        self.log(f"Wrote report to {self._report_file}")

    def get_endpoint(self) -> str:
        """
        Return path to pickle file (homogeneous with OpenEye backend).

        Returns:
            Path to pickle file containing RDKit molecules
        """
        if self._output_file is None:
            raise RuntimeError("No conformers generated yet. Call generate_conformers() first.")
        return self._output_file

    def extract_molecules(self) -> Iterator[Chem.Mol]:
        """
        Extract molecules from pickle file.

        Yields:
            RDKit molecules with conformers
        """
        if self._output_file is None:
            self.log("No conformers generated yet. Call generate_conformers() first.", level='ERROR')
            return

        # Check file exists
        if not Path(self._output_file).exists():
            self.log("Output file not found", level='ERROR')
            return

        # Read molecule data (binary, name) tuples from pickle file
        try:
            with open(self._output_file, 'rb') as f:
                mol_data = pickle.load(f)
        except Exception as e:
            self.log(f"Failed to read output file: {e}", level='ERROR')
            return

        # Convert from binary format and restore names
        for mol_binary, name in mol_data:
            try:
                mol = Chem.Mol(mol_binary)
                if name and not mol.HasProp('_Name'):
                    mol.SetProp('_Name', name)
                yield mol
            except Exception as e:
                self.log(f"Failed to load molecule from binary ({name}): {e}", level='WARNING')
                continue

    def get_report_dataframe(self) -> pd.DataFrame:
        """Get conformer generation report as DataFrame (matches OMEGA format)."""
        if not Path(self._report_file).exists():
            self.log("Report file not found", level='WARNING')
            return pd.DataFrame(columns=['Molecule', 'Title', 'Rotors', 'Conformers', 'ElapsedTime(s)', 'Status'])

        try:
            return pd.read_csv(self._report_file)
        except Exception as e:
            self.log(f"Failed to read report file: {e}", level='ERROR')
            return pd.DataFrame(columns=['Molecule', 'Title', 'Rotors', 'Conformers', 'ElapsedTime(s)', 'Status'])

    def get_successful_names(self) -> list[str]:
        """Get list of successfully generated molecule names from pickle file."""
        if self._output_file is None:
            self.log("Cannot extract names from file: no output file generated yet", level='WARNING')
            return []

        if not Path(self._output_file).exists():
            self.log("Cannot extract names from file: output file does not exist", level='WARNING')
            return []

        if Path(self._output_file).stat().st_size == 0:
            self.log("Cannot extract names from file: output file is empty", level='WARNING')
            return []

        # Read molecule names from pickle file
        try:
            with open(self._output_file, 'rb') as f:
                mol_data = pickle.load(f)

            # Extract names from (binary, name) tuples
            names = [name for _, name in mol_data]
            return names

        except Exception as e:
            self.log(f"Cannot extract names from file: failed to read output file: {e}", level='ERROR')
            return []

    def clear(self):
        """Clear stored molecules to free memory."""
        self._molecules.clear()
