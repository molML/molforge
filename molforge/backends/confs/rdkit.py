"""
RDKit ETKDG conformer generation backend.

Uses RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry) method
for conformer generation. Writes output in chunks for memory-efficient streaming.
"""

from typing import Dict, Any, Iterator, List, Tuple
import time
import pickle
import json
import gc
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem

from .base import ConformerBackend
from ...utils.constants import OUTPUT_CHUNK_SIZE, DEFAULT_N_JOBS


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
    - Chunked output for memory-efficient streaming
    - ETKDG conformer generation with parallel processing
    - Optional UFF optimization
    - Configurable pruning threshold
    """
    __backend_name__ = 'RDK'

    # Backend metadata
    description = "Generate conformers using RDKit ETKDGv3 algorithm"
    required_params = ['SMILES_column']
    optional_params = [
        'names_column', 'max_confs', 'rms_threshold', 'dropna', 'timeout',
        'use_random_coords', 'random_seed', 'num_threads',
        'use_uff', 'max_iterations',
    ]

    def __init__(self, params, logger=None, context=None):
        """Initialize RDKit backend."""
        super().__init__(params, logger, context)

        self._confs_dir: str = None

        self.log(
            f"RDKit ETKDG initialized: max_confs={params.max_confs}, "
            f"rms={params.rms_threshold}, optimize={params.use_uff}"
        )

    @property
    def _report_file(self) -> str:
        """Path to report CSV file."""
        if self._confs_dir:
            return str(Path(self._confs_dir) / "RDKit_rpt.csv")
        if self.context and hasattr(self.context, 'output_dir'):
            return str(Path(self.context.output_dir) / "conformers" / "RDKit_rpt.csv")
        return "RDKit_rpt.csv"

    def _setup_confs_dir(self) -> str:
        """Create and return the conformers output directory."""
        if self.context and hasattr(self.context, 'output_dir'):
            confs_dir = Path(self.context.output_dir) / "conformers"
        else:
            confs_dir = Path("conformers")

        confs_dir.mkdir(parents=True, exist_ok=True)
        return str(confs_dir)

    def generate_conformers(self, smiles_list: list[str], names_list: list[str]) -> Dict[str, Dict[str, Any]]:
        """
        Generate conformers using RDKit ETKDG with parallel processing and progress reporting.

        Writes successful molecules to chunked pickle files in input order.

        Args:
            smiles_list: List of SMILES strings
            names_list: List of molecule identifiers

        Returns:
            Dict mapping names to generation results
        """
        self._confs_dir = self._setup_confs_dir()

        n_mols = len(smiles_list)

        self.log(
            f"Generating conformers for {n_mols:,} molecules with {DEFAULT_N_JOBS} workers."
        )

        # Submit one molecule per future for dynamic load balancing
        results_dict = {}
        start_time = time.time()
        completed = 0

        with ProcessPoolExecutor(max_workers=DEFAULT_N_JOBS) as executor:
            future_to_name = {
                executor.submit(
                    _generate_single_conformer,
                    smiles, name,
                    self.params.max_confs,
                    self.params.random_seed,
                    self.params.rms_threshold,
                    self.params.use_random_coords,
                    self.params.use_uff,
                    self.params.max_iterations
                ): name
                for smiles, name in zip(smiles_list, names_list)
            }

            for future in as_completed(future_to_name, timeout=self.params.timeout):
                name = future_to_name[future]
                try:
                    result_name, result = future.result(timeout=self.params.timeout)
                    results_dict[result_name] = result
                except Exception as e:
                    self.log(f"Molecule '{name}' failed: {e}", level='ERROR')
                    results_dict[name] = {
                        'success': False,
                        'n_conformers': 0,
                        'status': f'Worker error: {str(e)}',
                        'rotors': 0,
                        'elapsed_time': 0.0
                    }

                completed += 1
                if completed % max(1, n_mols // 10) == 0 or completed == n_mols:
                    elapsed = time.time() - start_time
                    rate = completed / elapsed if elapsed > 0 else 0
                    eta = (n_mols - completed) / rate if rate > 0 else 0
                    self.log(
                        f"Progress: {completed:,}/{n_mols:,} ({100 * completed / n_mols:.0f}%) | "
                        f"{elapsed:.1f}s elapsed | ETA: {eta:.0f}s | {rate:.1f} mol/s"
                    )

        # Write molecules and report in input order
        report_data = []
        mol_buffer = []
        manifest = []
        output_chunk_idx = 0

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

            # Buffer successful molecules for chunked output
            if result['success'] and 'molecule' in result:
                mol = result['molecule']
                if not mol.HasProp('_Name'):
                    mol.SetProp('_Name', name)
                mol_buffer.append(mol)

                # Flush when buffer reaches output chunk size
                if len(mol_buffer) >= OUTPUT_CHUNK_SIZE:
                    chunk_entry = self._write_chunk(mol_buffer, output_chunk_idx)
                    manifest.append(chunk_entry)
                    mol_buffer = []
                    output_chunk_idx += 1
                    gc.collect()

        # Flush remaining molecules
        if mol_buffer:
            chunk_entry = self._write_chunk(mol_buffer, output_chunk_idx)
            manifest.append(chunk_entry)

        # Write manifest and report
        self._write_manifest(manifest)
        self._write_report_file(report_data)

        success_count = sum(1 for r in results_dict.values() if r['success'])
        total_mols = sum(entry['count'] for entry in manifest)
        total_time = time.time() - start_time
        rate = len(smiles_list) / total_time if total_time > 0 else 0
        self.log(
            f"RDKit generation complete: {success_count}/{len(smiles_list)} succeeded, "
            f"wrote {total_mols:,} molecules in {len(manifest)} chunk(s) "
            f"({total_time:.1f}s, {rate:.1f} mol/s)"
        )

        return results_dict

    def _write_chunk(self, molecules: List[Chem.Mol], chunk_idx: int) -> Dict[str, Any]:
        """
        Write a chunk of molecules to pickle file.

        Args:
            molecules: List of RDKit molecules with conformers
            chunk_idx: Chunk index for filename

        Returns:
            Manifest entry dict with file, names, and count
        """
        filename = f"chunk_{chunk_idx:04d}.pkl"
        chunk_path = Path(self._confs_dir) / filename

        mol_data = []
        names = []
        for mol in molecules:
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else ''
            mol_data.append((mol.ToBinary(), name))
            names.append(name)

        with open(chunk_path, 'wb') as f:
            pickle.dump(mol_data, f)

        return {'file': filename, 'names': names, 'count': len(mol_data)}

    def _write_manifest(self, manifest: List[Dict[str, Any]]) -> None:
        """
        Write manifest file for fast name lookup and integrity checking.

        Args:
            manifest: List of chunk entries with file, names, count
        """
        manifest_path = Path(self._confs_dir) / "manifest.json"
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)

    def _write_report_file(self, report_data: list[Dict[str, Any]]) -> None:
        """
        Write conformer generation report to CSV file.

        Args:
            report_data: List of report entries (dicts)
        """
        report_df = pd.DataFrame(report_data)
        Path(self._report_file).parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(self._report_file, index=False)
        self.log(f"Wrote report to {self._report_file}")

    def get_endpoint(self) -> str:
        """
        Return path to conformers directory.

        Returns:
            Path to directory containing chunked conformer output
        """
        if self._confs_dir is None:
            raise RuntimeError("No conformers generated yet. Call generate_conformers() first.")
        return self._confs_dir

    def extract_molecules(self) -> Iterator[Chem.Mol]:
        """
        Stream molecules from chunked pickle files.

        Loads one chunk at a time to avoid loading all molecules into memory.

        Yields:
            RDKit molecules with conformers, in input order
        """
        if self._confs_dir is None:
            self.log("No conformers generated yet. Call generate_conformers() first.", level='ERROR')
            return

        confs_dir = Path(self._confs_dir)
        chunk_files = sorted(confs_dir.glob("chunk_*.pkl"))

        if not chunk_files:
            self.log("No chunk files found", level='ERROR')
            return

        for chunk_path in chunk_files:
            try:
                with open(chunk_path, 'rb') as f:
                    mol_data = pickle.load(f)
            except Exception as e:
                self.log(f"Failed to read {chunk_path.name}: {e}", level='ERROR')
                continue

            for mol_binary, name in mol_data:
                try:
                    mol = Chem.Mol(mol_binary)
                    if name and not mol.HasProp('_Name'):
                        mol.SetProp('_Name', name)
                    yield mol
                except Exception as e:
                    self.log(f"Failed to load molecule from binary ({name}): {e}", level='WARNING')
                    continue

            del mol_data
            gc.collect()

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
        """
        Get list of successfully generated molecule names from manifest.

        Reads the lightweight manifest file and verifies each chunk file
        exists and is non-empty, without deserializing molecule data.

        Returns:
            List of molecule names that have conformers on disk
        """
        if self._confs_dir is None:
            self.log("Cannot extract names: no output generated yet", level='WARNING')
            return []

        manifest_path = Path(self._confs_dir) / "manifest.json"
        if not manifest_path.exists():
            self.log("Cannot extract names: manifest not found", level='WARNING')
            return []

        try:
            with open(manifest_path, 'r') as f:
                manifest = json.load(f)
        except Exception as e:
            self.log(f"Cannot extract names: failed to read manifest: {e}", level='ERROR')
            return []

        names = []
        for entry in manifest:
            chunk_path = Path(self._confs_dir) / entry['file']

            if not chunk_path.exists() or chunk_path.stat().st_size == 0:
                self.log(f"Chunk file missing or empty: {entry['file']}", level='WARNING')
                continue

            names.extend(entry['names'])

        return names
