"""
Torsion analysis actor plugin.

Analyzes torsional properties for molecules with conformers.
Requires phd-tools package for calculations.
"""

from typing import List, Tuple, Dict, Optional, Literal, Iterator, Any
import pickle
import multiprocessing as mp
from pathlib import Path
import gc
import os

import pandas as pd
from rdkit import Chem
from joblib import Parallel, delayed

# Attempt to import private toolkit
try:
    from phd_tools.chemistry.torsion import TorsionCalculator
    from phd_tools.chemistry.savebonds import CanonicalBondProperties
    from phd_tools.chemistry.convert import MoleculeConverter, HAS_OPENEYE
    TORSION_AVAILABLE = True
except ImportError:
    TORSION_AVAILABLE = False
    TorsionCalculator = None
    CanonicalBondProperties = None
    MoleculeConverter = None
    HAS_OPENEYE = False

if HAS_OPENEYE:
    from openeye import oechem


# Plugin imports
from molforge.actors.base import BaseActor
from molforge.actors.protocol import ActorOutput
from molforge.actors.params.base import BaseParams
from molforge.utils.constants import MAX_CHUNK_SIZE

from dataclasses import dataclass


# ============================================================================
# CONSTANTS
# ============================================================================

# Number of worker processes (leave one core for system)
N_WORKERS = max(1, mp.cpu_count() - 1)

# Chunk size for parallel processing (from molforge constants)
CHUNK_SIZE = MAX_CHUNK_SIZE


# ============================================================================
# PARAMETERS
# ============================================================================

@dataclass
class CalculateTorsionsParams(BaseParams):
    """
    Configuration for torsion analysis.

    Parameters control computational settings and analysis thresholds.
    See phd-tools documentation for detailed parameter descriptions.
    """

    aggregation_method: Literal['min', 'max', 'mean'] = 'mean'
    account_for_symmetry: bool = True
    symmetry_radius: int = 3
    ignore_colinear_bonds: bool = True
    n_confs_threshold: int = 50

    def _validate_params(self) -> None:
        """Validate parameter values."""
        self._validate_policy(
            'aggregation_method',
            self.aggregation_method,
            ['min', 'max', 'mean']
        )

        if self.n_confs_threshold < 1:
            raise ValueError("Conformer threshold must be at least 1")


# ============================================================================
# WORKER FUNCTION (module-level for pickling)
# ============================================================================

def _process_batch_worker(
    batch: List[Tuple[int, Any]],
    calculator_config: Dict,
) -> List[Tuple[int, Dict, Dict, bool, str]]:
    """
    Worker function for parallel batch processing.

    Processes a batch of molecules in a single worker, amortizing the
    overhead of process spawning and serialization across many molecules.

    Args:
        batch: List of (row_idx, molecule) tuples
        calculator_config: Configuration dict for TorsionCalculator

    Returns:
        List of (row_idx, canonical_mapping, stats, success, torsion_smiles) tuples
    """
    # Instantiate calculator once per worker (reused for all molecules in batch)
    calculator = TorsionCalculator(**calculator_config)

    results = []
    for row_idx, mol in batch:
        if mol is None:
            results.append((row_idx, {}, _empty_stats(), False, ""))
            continue

        try:
            # Compute torsions
            bond_variance, stats = calculator.compute(mol)

            # Get RDKit topology for canonicalization
            if HAS_OPENEYE and isinstance(mol, oechem.OEMol):
                rd_topology = MoleculeConverter.openeye_to_rdkit_topology(mol)
            else:
                rd_topology = mol

            # Canonicalize bond indices
            canonical_mapping = CanonicalBondProperties.bond_idx_to_canonical(
                rd_topology, bond_variance
            )

            # Generate SMILES with explicit hydrogens
            torsion_smiles = Chem.MolToSmiles(rd_topology)#, allHsExplicit=True)

            results.append((row_idx, canonical_mapping, stats, True, torsion_smiles))

        except Exception as e:
            error_stats = _empty_stats()
            error_stats['warnings'] = [f"{type(e).__name__}: {str(e)}"]
            results.append((row_idx, {}, error_stats, False, ""))

    return results


def _empty_stats() -> Dict:
    """Return empty statistics dict."""
    return {
        'n_confs': 0,
        'low_confs': False,
        'n_ring_torsions': 0,
        'n_rotor_torsions': 0,
        'n_torsions': 0,
        'mean_variance': 0.0,
        'warnings': []
    }


# ============================================================================
# ACTOR
# ============================================================================

class CalculateTorsions(BaseActor):
    """
    Torsion analysis actor.

    Analyzes torsional properties for molecules with multiple conformers.
    Outputs bond-level metrics and conformer statistics.

    Uses chunked parallel processing for large datasets. Results are saved
    incrementally to avoid memory issues with millions of molecules.

    Requires:
        - GenerateConfs actor output (conformers)
        - phd-tools package (private)

    Output columns:
        - torsion_smiles: SMILES with explicit hydrogens for molecule reconstruction
        - torsion_mapping: Canonical bond flexibility mapping {(rank1, rank2): variance} (serialized)
        - torsion_success: Boolean flag for successful torsion calculation
        - n_confs: Number of conformers analyzed
        - n_torsions: Total torsions detected
        - n_rotor_torsions: Rotatable bond torsions
        - n_ring_torsions: Ring torsions
        - mean_variance: Average metric across bonds
        - low_confs: Flag for low conformer count
        - warnings: Analysis warnings
    """

    __step_name__ = 'torsion'
    __param_class__ = CalculateTorsionsParams

    OUTPUT_COLUMNS = [
        'torsion_smiles', 'torsion_mapping', 'torsion_success', 'n_confs', 'low_confs',
        'n_ring_torsions', 'n_rotor_torsions', 'n_torsions',
        'mean_variance', 'warnings'
    ]

    @property
    def required_columns(self) -> List[str]:
        """Required columns from GenerateConfs."""
        return ['conformer_success']

    @property
    def output_columns(self) -> List[str]:
        """Columns added by torsion analysis."""
        return self.OUTPUT_COLUMNS

    @property
    def forge_endpoint(self) -> str:
        """Endpoint for MolForge integration. Points to torsion mapping column."""
        return 'torsion_mapping'

    @property
    def torsions_dir(self) -> str:
        """Path to torsions directory in current run directory."""
        dir_path = self._get_run_path("torsions")
        os.makedirs(dir_path, exist_ok=True)
        return dir_path

    def __post_init__(self):
        """Initialize torsion analysis."""
        if not TORSION_AVAILABLE:
            raise ImportError(
                "This actor requires the 'phd-tools' package.\n"
                "\n"
                "Installation:\n"
                "  pip install git+ssh://git@github.com/LukeRossen/phd-tools.git\n"
                "\n"
                "Note: Repository access required. Contact maintainer."
            )

        # Store calculator config for worker processes
        self._calculator_config = {
            'aggregation_method': self.aggregation_method,
            'account_for_symmetry': self.account_for_symmetry,
            'symmetry_radius': self.symmetry_radius,
            'ignore_colinear_bonds': self.ignore_colinear_bonds,
            'n_confs_threshold': self.n_confs_threshold,
        }

        self.log(
            f"Torsion analysis initialized "
            f"(workers: {N_WORKERS}, chunk_size: {CHUNK_SIZE:,})"
        )

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Analyze torsional properties for molecules with conformers.

        Processes molecules in parallel chunks to handle large datasets
        efficiently while managing memory.

        Args:
            data: DataFrame with conformer_success column

        Returns:
            DataFrame with torsion analysis columns added
        """
        df = data.copy()

        gc_actor = self.get_actor('confs')
        if gc_actor is None:
            raise ValueError(
                "GenerateConfs actor not found. "
                "Ensure GenerateConfs runs before CalculateTorsions."
            )

        # Get molecule generator (streams molecules, no memory load)
        mol_generator = gc_actor.extract_molecules()

        # Process molecules in chunks and save incrementally
        torsion_stats = self._process_chunked(
            mol_generator,
            df['conformer_success'],
            len(df),
        )

        # Assign results to DataFrame
        for col, values in torsion_stats.items():
            df[col] = values

        return df

    def _process_chunked(
        self,
        molecules: Iterator,
        success_mask: pd.Series,
        total_rows: int,
    ) -> Dict[str, List]:
        """
        Process molecules in parallel chunks with incremental output.

        Args:
            molecules: Generator yielding successful molecules
            success_mask: Boolean series indicating successful conformers
            total_rows: Total number of rows to process

        Returns:
            Dictionary mapping column names to lists of values
        """
        # Initialize result collectors
        results = {col: [] for col in self.OUTPUT_COLUMNS}

        # Track chunk statistics
        chunk_idx = 0
        processed_count = 0

        # Collect molecules with their row indices
        chunk_data = []  # List of (row_idx, mol or None)

        for row_idx, has_conformers in enumerate(success_mask):
            if has_conformers:
                try:
                    mol = next(molecules)
                    chunk_data.append((row_idx, mol))
                except StopIteration:
                    self.log(f"Generator exhausted at row {row_idx}", level='ERROR')
                    chunk_data.append((row_idx, None))
            else:
                # No conformer - will use empty result
                chunk_data.append((row_idx, None))

            # Process chunk when full or at end
            if len(chunk_data) >= CHUNK_SIZE or row_idx == total_rows - 1:
                if chunk_data:
                    chunk_results = self._process_chunk(chunk_data, chunk_idx)
                    self._append_results(results, chunk_results)

                    processed_count += len(chunk_data)
                    self.log(
                        f"Processed chunk {chunk_idx}: "
                        f"{processed_count:,}/{total_rows:,} molecules"
                    )

                    chunk_data = []
                    chunk_idx += 1
                    gc.collect()

        self.log(f"Completed: {processed_count:,} molecules in {chunk_idx} chunks")
        return results

    def _process_chunk(
        self,
        chunk_data: List[Tuple[int, Any]],
        chunk_idx: int,
    ) -> List[Tuple]:
        """
        Process a chunk of molecules in parallel.

        Splits molecules into N_WORKERS batches for efficient parallelization,
        minimizing serialization overhead by processing many molecules per worker.

        Args:
            chunk_data: List of (row_idx, molecule) tuples
            chunk_idx: Chunk index for logging/saving

        Returns:
            List of result tuples, ordered by row_idx
        """
        # Separate molecules that need processing from empty ones
        to_process = [(idx, mol) for idx, mol in chunk_data if mol is not None]
        empty_indices = {idx for idx, mol in chunk_data if mol is None}

        # Process non-empty molecules in parallel using batch-per-worker
        processed_results = {}

        if to_process:
            self.log(f"Splitting {len(to_process)} items into {N_WORKERS} batches of size ~{len(to_process) // N_WORKERS}")
            # Split into N_WORKERS batches for efficient parallelization
            batches = self._split_into_batches(to_process, N_WORKERS)

            # Each worker processes a full batch (thousands of molecules)
            batch_results = Parallel(n_jobs=N_WORKERS)(
                delayed(_process_batch_worker)(batch, self._calculator_config)
                for batch in batches
            )

            # Flatten and map results by row_idx
            for batch_result in batch_results:
                for result in batch_result:
                    row_idx = result[0]
                    processed_results[row_idx] = result

        # Combine results in original order
        ordered_results = []
        chunk_molecules = []
        chunk_mappings = []

        for row_idx, _ in chunk_data:
            if row_idx in empty_indices:
                # Empty result for failed/missing conformers
                result = (row_idx, {}, _empty_stats(), False, "")
            else:
                result = processed_results.get(
                    row_idx,
                    (row_idx, {}, _empty_stats(), False, "")
                )

            ordered_results.append(result)

            # Collect successful results for pickle
            if result[3]:  # success flag
                _, mapping, _, _, smiles = result
                # Parse SMILES to mol for pickle
                params = Chem.SmilesParserParams()
                params.removeHs = False
                mol_obj = Chem.MolFromSmiles(smiles, params)
                if mol_obj is not None:

                    # Convert canonical mapping back to bond indices for storage
                    bond_map = CanonicalBondProperties.canonical_to_bond_idx(
                        mol_obj, mapping
                    )
                    chunk_molecules.append(mol_obj)
                    chunk_mappings.append(bond_map)

        # Save chunk to pickle
        if chunk_molecules:
            self._save_chunk(chunk_molecules, chunk_mappings, chunk_idx)

        return ordered_results

    def _append_results(
        self,
        results: Dict[str, List],
        chunk_results: List[Tuple],
    ) -> None:
        """
        Append chunk results to main results dictionary.

        Args:
            results: Main results dictionary to append to
            chunk_results: List of result tuples from chunk processing
        """
        for _, mapping, stats, success, smiles in chunk_results:
            results['torsion_smiles'].append(smiles)
            results['torsion_mapping'].append(mapping.__repr__())
            results['torsion_success'].append(success)
            results['n_confs'].append(stats['n_confs'])
            results['low_confs'].append(stats['low_confs'])
            results['n_ring_torsions'].append(stats['n_ring_torsions'])
            results['n_rotor_torsions'].append(stats['n_rotor_torsions'])
            results['n_torsions'].append(stats['n_torsions'])
            results['mean_variance'].append(stats['mean_variance'])
            results['warnings'].append(stats['warnings'])

    @staticmethod
    def _split_into_batches(
        items: List[Any],
        n_batches: int,
    ) -> List[List[Any]]:
        """
        Split items into n approximately equal batches.

        Args:
            items: List of items to split
            n_batches: Number of batches to create

        Returns:
            List of batches (lists)
        """
        if not items:
            return []

        # Ensure at least 1 batch
        n_batches = max(1, min(n_batches, len(items)))
        batch_size = (len(items) + n_batches - 1) // n_batches
    
        return [
            items[i:i + batch_size]
            for i in range(0, len(items), batch_size)
        ]

    def _save_chunk(
        self,
        molecules: List[Chem.Mol],
        mappings: List[Dict],
        chunk_idx: int,
    ) -> None:
        """
        Save chunk results to pickle file.

        Args:
            molecules: List of RDKit molecules
            mappings: List of bond variance mappings
            chunk_idx: Chunk index for filename
        """
        chunk_path = os.path.join(self.torsions_dir, f"chunk_{chunk_idx:04d}.pkl")

        with open(chunk_path, 'wb') as f:
            pickle.dump((molecules, mappings), f)

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with torsion mapping endpoint."""
        if 'n_torsions' in data.columns and 'torsion_success' in data.columns:
            success_df = data[data['torsion_success']]
            total_torsions = success_df['n_torsions'].sum() if len(success_df) > 0 else 0
            mean_metric = success_df['mean_variance'].mean() if len(success_df) > 0 else 0.0
            n_success = len(success_df)
        else:
            total_torsions = 0
            mean_metric = 0.0
            n_success = 0

        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'n_molecules': len(data),
                'n_success': n_success,
                'total_torsions': int(total_torsions),
                'mean_metric': float(mean_metric),
                'n_chunks': len(list(Path(self.torsions_dir).glob("chunk_*.pkl"))),
            },
            endpoint=self.forge_endpoint
        )

    # ========================================================================
    # RESULT ACCESS METHODS
    # ========================================================================

    def extract_results(self) -> Iterator[Tuple[Chem.Mol, Dict[int, float]]]:
        """
        Stream torsion results from cached pickle chunks.

        Yields results one at a time to handle large datasets without
        loading all into memory.

        Yields:
            (molecule, bond_mapping) tuples where molecules are RDKit Mol
            objects and mappings are {bond_idx: variance} dicts

        Raises:
            FileNotFoundError: If no chunk files exist.
        """
        chunk_files = sorted(Path(self.torsions_dir).glob("chunk_*.pkl"))

        if not chunk_files:
            raise FileNotFoundError(
                f"No torsion results found in {self.torsions_dir}.\n"
                "Ensure the torsion actor has been run in this pipeline."
            )

        self.log(f"Streaming results from {len(chunk_files)} chunks")

        total_yielded = 0
        for chunk_path in chunk_files:
            with open(chunk_path, 'rb') as f:
                molecules, mappings = pickle.load(f)

            for mol, mapping in zip(molecules, mappings):
                yield mol, mapping
                total_yielded += 1

            # Free memory after each chunk
            del molecules, mappings
            gc.collect()

        self.log(f"Streamed {total_yielded:,} molecules")

    def extract_results_batch(
        self,
        batch_size: Optional[int] = None,
    ) -> Iterator[Tuple[List[Chem.Mol], List[Dict[int, float]]]]:
        """
        Stream torsion results in batches.

        Args:
            batch_size: Number of molecules per batch. If None, yields
                one chunk at a time (efficient for chunk-aligned processing).

        Yields:
            (molecules_batch, mappings_batch) tuples
        """
        chunk_files = sorted(Path(self.torsions_dir).glob("chunk_*.pkl"))

        if not chunk_files:
            raise FileNotFoundError(f"No torsion results found in {self.torsions_dir}")

        if batch_size is None:
            # Yield whole chunks
            for chunk_path in chunk_files:
                with open(chunk_path, 'rb') as f:
                    molecules, mappings = pickle.load(f)
                yield molecules, mappings
        else:
            # Yield custom batch sizes
            mol_buffer = []
            map_buffer = []

            for chunk_path in chunk_files:
                with open(chunk_path, 'rb') as f:
                    molecules, mappings = pickle.load(f)

                for mol, mapping in zip(molecules, mappings):
                    mol_buffer.append(mol)
                    map_buffer.append(mapping)

                    if len(mol_buffer) >= batch_size:
                        yield mol_buffer, map_buffer
                        mol_buffer = []
                        map_buffer = []

                del molecules, mappings
                gc.collect()

            # Yield remaining
            if mol_buffer:
                yield mol_buffer, map_buffer

    def get_results_count(self) -> int:
        """
        Get total count of successful results without loading data.

        Returns:
            Number of molecules with torsion results
        """
        count = 0
        for chunk_path in Path(self.torsions_dir).glob("chunk_*.pkl"):
            with open(chunk_path, 'rb') as f:
                molecules, _ = pickle.load(f)
                count += len(molecules)

        return count
