"""
Torsion analysis actor plugin.

Analyzes torsional properties for molecules with conformers.
Requires phd-tools package for calculations.
"""

from typing import List, Tuple, Dict, Optional, Literal, Iterator, Any
import pickle
from pathlib import Path
import json
import time
import gc
import os

import pandas as pd
from rdkit import Chem
from joblib import Parallel, delayed

# Attempt to import private toolkit
try:
    from phd_tools.chemistry.torsion import TorsionCalculator
    from phd_tools.chemistry.savebonds import CanonicalBondProperties
    TORSION_AVAILABLE = True
except ImportError:
    TORSION_AVAILABLE = False
    TorsionCalculator = None
    CanonicalBondProperties = None


# Plugin imports
from molforge.actors.base import BaseActor
from molforge.actors.protocol import ActorOutput
from molforge.actors.params.base import BaseParams
from molforge.utils.constants import MAX_CHUNK_SIZE, DEFAULT_MP_THRESHOLD, DEFAULT_N_JOBS

from dataclasses import dataclass


# ============================================================================
# CONSTANTS
# ============================================================================

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
    min_confs: Optional[int] = None
    """Hard minimum conformer count. Molecules below this are treated as
    failures and removed from results. None disables this filter.
    Distinct from n_confs_threshold, which only flags an informational warning."""
    dropna: bool = True
    """Drop rows where torsion analysis failed. Failures include: torsion
    computation errors (malformed molecule, no torsions), SMILES reconstruction
    failures, bond canonicalization errors, and min_confs violations.
    Rows without upstream conformers are also dropped."""

    def _validate_params(self) -> None:
        """Validate parameter values."""
        self._validate_policy(
            'aggregation_method',
            self.aggregation_method,
            ['min', 'max', 'mean']
        )

        if self.n_confs_threshold < 1:
            raise ValueError("Conformer threshold must be at least 1")

        if self.min_confs is not None and self.min_confs < 1:
            raise ValueError("min_confs must be at least 1")


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
    # Prevent PyTorch OpenMP thread oversubscription in forked workers
    import torch
    torch.set_num_threads(1)

    # Instantiate calculator once per worker (reused for all molecules in batch)
    calculator = TorsionCalculator(**calculator_config)

    results = []
    for row_idx, mol in batch:
        if mol is None:
            results.append((row_idx, {}, _empty_stats(), False, ""))
            continue

        try:
            # Compute torsions (rd_topology returned to avoid duplicate conversion)
            bond_variance, stats, rd_topology = calculator.compute(mol)

            # Canonicalize bond indices
            canonical_mapping = CanonicalBondProperties.bond_idx_to_canonical(
                rd_topology, bond_variance
            )

            # Generate SMILES for molecule reconstruction
            torsion_smiles = Chem.MolToSmiles(rd_topology)

            results.append((row_idx, canonical_mapping, stats, True, torsion_smiles))

        except Exception as e:
            error_stats = _empty_stats()
            error_stats['warnings'] = [f"Torsion computation: {type(e).__name__}: {str(e)}"]
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
        - torsion_mapping: JSON with 'canonical_atom_ranks' and 'circular_variance'
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
            f"(workers: {DEFAULT_N_JOBS}, chunk_size: {CHUNK_SIZE:,})"
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

        # Intake summary
        n_with_confs = int(df['conformer_success'].sum())
        n_without_confs = len(df) - n_with_confs
        self.log(
            f"Received {len(df):,} molecules "
            f"({n_with_confs:,} with conformers, {n_without_confs:,} without)."
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

        # Report: torsion-specific failures (not upstream)
        failed = df[~df['torsion_success'] & df['conformer_success']]
        if len(failed) > 0:
            failure_types = failed['warnings'].apply(
                lambda w: w[0].split(':')[0] if w else 'Unknown'
            )
            summary = failure_types.value_counts()
            self.log(
                f"Failure summary ({len(failed):,} torsion failures):\n"
                f"{summary.to_frame('count').to_markdown()}",
                level='WARNING'
            )

        # Report: low conformer warnings (informational, not failures)
        successful = df[df['torsion_success']]
        n_low_confs = int(successful['low_confs'].sum()) if len(successful) > 0 else 0
        if n_low_confs > 0:
            self.log(
                f"{n_low_confs:,} successful molecules flagged with low conformer "
                f"count (< {self.n_confs_threshold}). Results may be less "
                f"statistically reliable for these molecules."
            )

        # Drop failed rows if requested
        if self.dropna:
            initial_count = len(df)
            df = df[df['torsion_success']]
            dropped = initial_count - len(df)
            if dropped > 0:
                self.log(f"Dropped {dropped:,} rows with failed torsion analysis")

        return df

    def _process_chunked(
        self,
        molecules: Iterator,
        success_mask: pd.Series,
        total_rows: int,
    ) -> Dict[str, List]:
        """
        Process molecules in chunks with optional parallelism.

        Uses single-process mode for small datasets (< DEFAULT_MP_THRESHOLD)
        and parallel processing for larger ones.

        Args:
            molecules: Generator yielding successful molecules
            success_mask: Boolean series indicating successful conformers
            total_rows: Total number of rows to process

        Returns:
            Dictionary mapping column names to lists of values
        """
        results = {col: [] for col in self.OUTPUT_COLUMNS}
        n_with_conformers = int(success_mask.sum())
        use_mp = n_with_conformers >= DEFAULT_MP_THRESHOLD
        total_chunks = (total_rows + CHUNK_SIZE - 1) // CHUNK_SIZE

        # Dispatch log (matches curate actor format)
        if use_mp:
            self.log(
                f"Processing {total_rows:,} molecules ({n_with_conformers:,} with conformers) "
                f"with {DEFAULT_N_JOBS} workers ({total_chunks} chunks of {CHUNK_SIZE:,})."
            )
        else:
            self.log(
                f"Processing {total_rows:,} molecules ({n_with_conformers:,} with conformers) "
                f"in single process."
            )

        chunk_idx = 0
        processed_count = 0
        chunk_data = []
        pipeline_start = time.time()
        chunk_times = []

        for row_idx, has_conformers in enumerate(success_mask):
            if has_conformers:
                try:
                    mol = next(molecules)
                    chunk_data.append((row_idx, mol))
                except StopIteration:
                    self.log(f"Generator exhausted at row {row_idx}", level='ERROR')
                    chunk_data.append((row_idx, None))
            else:
                chunk_data.append((row_idx, None))

            # Process chunk when full or at end
            if len(chunk_data) >= CHUNK_SIZE or row_idx == total_rows - 1:
                if chunk_data:
                    chunk_start = time.time()
                    chunk_results = self._process_chunk(
                        chunk_data, chunk_idx, use_mp
                    )
                    chunk_time = time.time() - chunk_start
                    chunk_times.append(chunk_time)

                    self._append_results(results, chunk_results)
                    processed_count += len(chunk_data)

                    # Per-chunk progress
                    elapsed = time.time() - pipeline_start
                    rate = processed_count / elapsed if elapsed > 0 else 0
                    self.log(
                        f"Chunk {chunk_idx + 1}/{total_chunks} | "
                        f"{processed_count:,}/{total_rows:,} "
                        f"({100 * processed_count / total_rows:.1f}%) | "
                        f"Rate: {rate:.0f} mol/s | Chunk: {chunk_time:.1f}s"
                    )

                    chunk_data = []
                    chunk_idx += 1
                    gc.collect()

        # Completion summary
        total_time = time.time() - pipeline_start
        avg_chunk = sum(chunk_times) / len(chunk_times) if chunk_times else 0
        self.log(
            f"Completed: {processed_count:,} molecules in {chunk_idx} chunks "
            f"({total_time:.1f}s total, avg {avg_chunk:.1f}s/chunk)"
        )
        return results

    def _process_chunk(
        self,
        chunk_data: List[Tuple[int, Any]],
        chunk_idx: int,
        use_mp: bool = True,
    ) -> List[Tuple]:
        """
        Process a chunk of molecules, optionally in parallel.

        When use_mp is True, splits molecules into DEFAULT_N_JOBS batches
        for parallel processing. Otherwise processes sequentially in-process.

        Args:
            chunk_data: List of (row_idx, molecule) tuples
            chunk_idx: Chunk index for logging/saving
            use_mp: Whether to use multiprocessing

        Returns:
            List of result tuples, ordered by row_idx
        """
        to_process = [(idx, mol) for idx, mol in chunk_data if mol is not None]
        empty_indices = {idx for idx, mol in chunk_data if mol is None}

        processed_results = {}

        if to_process:
            if use_mp:
                batches = self._split_into_batches(to_process, DEFAULT_N_JOBS)
                batch_results = Parallel(n_jobs=DEFAULT_N_JOBS)(
                    delayed(_process_batch_worker)(batch, self._calculator_config)
                    for batch in batches
                )
            else:
                # Single-process: call worker directly (no IPC overhead)
                batch_results = [
                    _process_batch_worker(to_process, self._calculator_config)
                ]

            for batch_result in batch_results:
                for result in batch_result:
                    row_idx = result[0]
                    processed_results[row_idx] = result

        # Combine results in original order, validating successful ones
        ordered_results = []
        chunk_molecules = []
        chunk_mappings = []

        smiles_params = Chem.SmilesParserParams()
        smiles_params.removeHs = False

        for row_idx, _ in chunk_data:
            if row_idx in empty_indices:
                result = (row_idx, {}, _empty_stats(), False, "")
            else:
                result = processed_results.get(
                    row_idx,
                    (row_idx, {}, _empty_stats(), False, "")
                )

            # Validate and pickle successful results
            if result[3]:  # success flag from worker
                _, mapping, stats, _, smiles = result

                # Check minimum conformer requirement (hard cutoff)
                if self.min_confs is not None and stats['n_confs'] < self.min_confs:
                    stats = dict(stats)
                    stats['warnings'] = stats.get('warnings', []) + [
                        f"Minimum conformers: {stats['n_confs']} < {self.min_confs}"
                    ]
                    result = (row_idx, {}, stats, False, "")
                    ordered_results.append(result)
                    continue

                mol_obj = Chem.MolFromSmiles(smiles, smiles_params)

                if mol_obj is None:
                    stats = dict(stats)
                    stats['warnings'] = stats.get('warnings', []) + [
                        'SMILES reconstruction: MolFromSmiles returned None'
                    ]
                    result = (row_idx, {}, stats, False, "")
                else:
                    try:
                        bond_map = CanonicalBondProperties.canonical_to_bond_idx(
                            mol_obj, mapping
                        )
                        chunk_molecules.append(mol_obj)
                        chunk_mappings.append(bond_map)
                    except (ValueError, KeyError) as e:
                        stats = dict(stats)
                        stats['warnings'] = stats.get('warnings', []) + [
                            f"Bond canonicalization: {e}"
                        ]
                        result = (row_idx, {}, stats, False, "")

            ordered_results.append(result)

        # Save chunk to pickle
        if chunk_molecules:
            self._save_chunk(chunk_molecules, chunk_mappings, chunk_idx)

        return ordered_results

    @staticmethod
    def _serialize_mapping(mapping: Dict) -> str:
        """
        Serialize canonical mapping to JSON for CSV-safe storage.

        Args:
            mapping: {(canonical_rank1, canonical_rank2): circular_variance}

        Returns:
            JSON string with keys 'canonical_atom_ranks' and 'circular_variance'
        """
        if not mapping:
            return '{}'
        return json.dumps({
            'canonical_atom_ranks': [list(k) for k in mapping.keys()],
            'circular_variance': list(mapping.values()),
        })

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
            results['torsion_mapping'].append(self._serialize_mapping(mapping))
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
