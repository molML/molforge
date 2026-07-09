"""
Molecule curation and standardization module using RDKit.

This module provides comprehensive molecule curation functionality including
desalting, tautomer canonicalization, neutralization, and stereochemistry handling.
"""

from typing import Tuple, List, Optional, Any, Callable, Union
from contextlib import contextmanager
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed
import signal
import time
import gc
import multiprocessing as mp

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.rdmolops import RemoveHs, GetFormalCharge
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator, CleanupParameters

from .params.curate import CurateMolParams
from .base import BaseActor
from .protocol import ActorOutput
from ..configuration.steps import Steps

from ..utils.mols.stereo import StereoHandler
from ..utils.actortools.multiprocess import multiprocess_worker, calculate_chunk_params
from ..utils.constants import DEFAULT_MP_THRESHOLD, MAX_CHUNK_SIZE


class CurateMol(BaseActor):
    """Molecule curation and standardization."""
    __step_name__ = Steps.CURATE
    """
    RDKit-based SMILES/molecule standardization and curation.
    
    This class provides a comprehensive pipeline for molecular curation including:
    - Salt removal
    - Tautomer canonicalization
    - Charge neutralization
    - Stereochemistry handling
    - Parallel processing for large datasets
    """
    
    # Class constants (imported from utils.constants)
    DEFAULT_MP_THRESHOLD = DEFAULT_MP_THRESHOLD
    MAX_CHUNK_SIZE = MAX_CHUNK_SIZE
    DEFAULT_PROGRESS_INTERVAL = 10
    DEFAULT_CHUNK_PROGRESS_INTERVAL = 50_000

    OUTPUT_COLUMNS = ['curated_smiles', 'curation_success', 'curation_step']
    
    # Neutralization patterns (from Hans de Winter's RDKit contribution)
    NEUTRALIZATION_PATTERNS = (
        ('[n+;H]', 'n'),                       # Imidazoles
        ('[N+;!H0]', 'N'),                     # Amines
        ('[$([O-]);!$([O-][#7])]', 'O'),       # Carboxylic acids and alcohols
        ('[S-;X1]', 'S'),                      # Thiols
        ('[$([N-;X2]S(=O)=O)]', 'N'),          # Sulfonamides
        ('[$([N-;X2][C,N]=C)]', 'N'),          # Enamines
        ('[n-]', '[nH]'),                      # Tetrazoles
        ('[$([S-]=O)]', 'S'),                  # Sulfoxides
        ('[$([N-]C=O)]', 'N'),                 # Amides
    )

    @property
    def required_columns(self) -> List[str]:
        """Required input columns."""
        return [self.SMILES_column]

    @property
    def output_columns(self) -> List[str]:
        """Output columns from curation."""
        return self.OUTPUT_COLUMNS

    @property
    def forge_endpoint(self) -> str:
        """Endpoint for MolForge integration. Points to curated SMILES column."""
        return self.OUTPUT_COLUMNS[0]

    def __post_init__(self):
        """Post-initialization setup after BaseActor.__init__."""
        self._log_configuration()
        self._initialize_components()
        self.output_cols = self.OUTPUT_COLUMNS

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Perform standardized molecule-level curation on a pandas DataFrame.

        Args:
            data: Input DataFrame with SMILES column

        Returns:
            Curated DataFrame with standardized molecules

        Raises:
            ValueError: If SMILES column not found in DataFrame
        """
        if self.SMILES_column not in data.columns:
            raise ValueError(
                f"'{self.SMILES_column}' not in dataframe. "
                f"Available columns: {list(data.columns)}"
            )

        curated_dataset = self.curate_dataframe(data)

        if self.check_duplicates:
            curated_dataset = self._post_process_duplicates(curated_dataset)

        return curated_dataset

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with endpoint."""
        return ActorOutput(
            data=data,
            success=True,
            endpoint=self.forge_endpoint,
            metadata={
                'curated_column': self.OUTPUT_COLUMNS[0],
                'n_curated': len(data),
                'n_successful': len(data[data['curation_success']]) if 'curation_success' in data.columns else None
            }
        )
    
    # ==================== Main Curation Methods ====================
    
    def curate_dataframe(self, df: pd.DataFrame, mp_threshold: Optional[int] = None,
        n_processes: Optional[int] = None, chunk_size: Optional[int] = None,
        progress_interval: Optional[int] = None) -> pd.DataFrame:
        """
        Curate DataFrame with automatic multiprocessing for large datasets.
        
        Args:
            df: Input DataFrame
            mp_threshold: Row count threshold for multiprocessing
            n_processes: Number of processes (default: CPU count - 1)
            chunk_size: Molecules per chunk (default: len(df) / n_processes)
            progress_interval: Progress report interval in seconds
            
        Returns:
            Curated DataFrame
        """
        df = df.copy()  # Avoid modifying original
        
        mp_threshold = mp_threshold or self.DEFAULT_MP_THRESHOLD
        progress_interval = progress_interval or self.DEFAULT_PROGRESS_INTERVAL
        
        if len(df) < mp_threshold:
            self.log(f"Processing {len(df):,} molecules with single process.")
            return self._curate_dataframe_single_process(df)
        
        self.log(f"Processing {len(df):,} molecules with multiprocessing.")
        return self._curate_dataframe_multiprocess(
            df, n_processes, chunk_size, progress_interval
        )
    
    def curate_molecule(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Curate a single molecule through the pipeline."""
        mol, success, _ = self._try_curate_molecule(mol, return_info=True)
        return mol if success else None
    
    def curate_smiles(self, smiles: str) -> Optional[str]:
        """Curate a single SMILES string through the pipeline."""
        smiles, success, _ = self._try_curate_smiles(smiles, return_info=True)
        return smiles if success else None
    
    def curate_molecule_to_smiles(self, mol: Chem.Mol) -> Optional[str]:
        """Curate a molecule and return as SMILES."""
        smiles, success, _ = self._try_curate_mol_to_smiles(mol, return_info=True)
        return smiles if success else None
    
    def curate_smiles_array(self, smiles_list: List[str]) -> List[Optional[str]]:
        """
        Curate a list of SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            List of curated SMILES (None for failures if dropna=False)
        """
        curated = [self.curate_smiles(smi) for smi in smiles_list]
        
        if self.dropna:
            curated = [smi for smi in curated if smi is not None]
        
        return curated
    
    # ==================== Curation Pipeline Methods ====================
    
    def _try_curate_molecule(self, mol: Chem.Mol, return_info: bool = True) -> Union[Chem.Mol, Tuple[Chem.Mol, bool, str]]:
        """
        Main curation pipeline for molecules.
        
        Args:
            mol: RDKit molecule object
            return_info: Whether to return metadata
            
        Returns:
            Molecule or tuple (molecule, success, last_step)
        """
        success = mol is not None
        step = 'start'
        original_smiles = self._get_smiles_safe(mol, 'Invalid Molecule')
        
        for step in self.mol_steps:
            if not success:
                break
            
            try:
                with self._timeout_context(self.step_timeout, step):
                    mol, success = self._execute_curation_step(mol, step)
                    
            except TimeoutError as e:
                self._handle_timeout(step, original_smiles, mol, e)
                mol, success = None, False
                break
                
            except Exception as e:
                self._handle_curation_error(step, original_smiles, mol, e)
                mol, success = None, False
                break
        
        if return_info:
            return mol, success, step if not success else 'finished'
        return mol
    
    def _try_curate_smiles(self, smiles: str, return_info: bool = False) -> Union[str, List[str], Tuple[Any, bool, str]]:
        """Curate SMILES string through the pipeline."""
        mol = Chem.MolFromSmiles(smiles)
        mol, success, step = self._try_curate_molecule(mol, return_info=True)
        
        if success:
            curated_smiles, success, step = self._mol_to_smiles_step(mol, return_info=True)
        else:
            curated_smiles = None
        
        if return_info:
            return curated_smiles, success, step
        return curated_smiles
    
    def _try_curate_mol_to_smiles(self, mol: Chem.Mol, return_info: bool = False) -> Union[str, Tuple[str, bool, str]]:
        """Curate molecule and convert to SMILES."""
        mol, success, step = self._try_curate_molecule(mol, return_info=True)
        
        if success:
            curated_smiles, success, step = self._mol_to_smiles_step(mol, return_info=True)
        else:
            curated_smiles = None
        
        if return_info:
            return curated_smiles, success, step
        return curated_smiles
    
    def _try_curate_smiles_array(self, smiles_list: List[str], **kwargs) -> List[Any]:
        """Curate array of SMILES with optional progress tracking."""
        if len(smiles_list) >= 10_000:
            self.log(f"Tracking SMILES progress (0/{len(smiles_list)}).", level="DEBUG")
            return self._with_progress_tracking(
                self._try_curate_smiles, 
                smiles_list, 
                **kwargs
            )
        
        return [self._try_curate_smiles(smi, **kwargs) for smi in smiles_list]
    
    # ==================== Individual Curation Steps ====================
    
    def _execute_curation_step(self, mol: Chem.Mol, step: str) -> Tuple[Chem.Mol, bool]:
        """Execute a specific curation step."""
        step_methods = {
            "desalt": self._desalt_molecule,
            "removeIsotope": self._remove_isotope,
            "removeHs": self._remove_hydrogens,
            "tautomers": self._canonicalize_tautomer,
            "neutralize": self._neutralize_molecule,
            "sanitize": self._sanitize_molecule,
            "computeProps": self._compute_properties,
            "handleStereo": self._handle_stereochemistry,
        }
        
        if step not in step_methods:
            raise ValueError(f"Unknown curation step: {step}")
        
        return step_methods[step](mol)
    
    def _desalt_molecule(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Remove salts from molecule."""
        if '.' not in Chem.MolToSmiles(mol):
            return mol, True
        
        if self.desalt_policy == 'remove':
            return None, False
        
        if self.desalt_policy == 'keep':
            mol = self.salt_remover.StripMol(mol, dontRemoveEverything=False)
            desalted_smiles = Chem.MolToSmiles(mol)
            
            if '.' in desalted_smiles:
                if self.brute_force_desalt:
                    # Keep largest fragment
                    mol = max(
                        Chem.GetMolFrags(mol, asMols=True), 
                        key=lambda x: x.GetNumAtoms()
                    )
                else:
                    return None, False
            
            if mol and mol.GetNumAtoms() > 0:
                return mol, True
        
        return None, False
    
    def _canonicalize_tautomer(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Canonicalize tautomers for consistent representation."""
        try:
            mol = self.tautomer_enumerator.Canonicalize(mol)
            return mol, True
        except Exception as e:
            if "Invariant Violation" in str(e):
                self.log(
                    f"Tautomer Invariant Violation for {Chem.MolToSmiles(mol)}", 
                    level='DEBUG'
                )
            else:
                self.log(f"Tautomer canonicalization error: {e}", level='DEBUG')
            return mol, False
    
    def _neutralize_molecule(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Neutralize charged groups in molecule."""
        if GetFormalCharge(mol) == 0 and self.neutralize_policy == 'keep':
            return mol, True
        
        # TODO: consider adding a policy for type of neutralization aproach (reaction-based/Cookbook-atomistic/rdMolStandardize.Uncharger)
        for reactant, product in self.neutralization_reactions:
            while mol.HasSubstructMatch(reactant):
                rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]
        
        return mol, GetFormalCharge(mol) == 0
    
    def _sanitize_molecule(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Sanitize molecule using RDKit."""
        sanitize_flags = Chem.SanitizeFlags.SANITIZE_ALL
        sanitize_fail = Chem.SanitizeMol(
            mol, 
            catchErrors=True, 
            sanitizeOps=sanitize_flags
        )
        
        if sanitize_fail:
            self.log(f"Sanitization failed: {sanitize_fail}", level='DEBUG')
            return mol, False
        
        return mol, True
    
    def _handle_stereochemistry(self, mol: Chem.Mol) -> Tuple[Union[Chem.Mol, List[Chem.Mol]], bool]:
        """Handle stereochemistry using configured StereoHandler."""
        mols, success = self.StereoHandler(mol)
        # Return single mol or list depending on enumerate policy
        if len(mols) == 1:
            return mols[0], success
        return mols, success
    
    def _remove_isotope(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Remove isotope information from atoms."""
        try:
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            return mol, True
        except Exception as e:
            self.log(f"Remove isotope error: {e}", level='DEBUG')
            return mol, False
    
    def _remove_hydrogens(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Remove explicit hydrogens from molecule."""
        # TODO: refactor to handle_hydrogens, allow for add or remove. 
        params = Chem.rdmolops.RemoveHsParameters()
        params.removeWithWedgedBond = True
        params.removeDefiningBondStereo = True
        params.updateExplicitCount = True
        
        mol = RemoveHs(mol, params)
        
        if mol:
            has_no_h = all(atom.GetSymbol() != 'H' for atom in mol.GetAtoms())
            return mol, has_no_h
        
        return mol, False
    
    def _compute_properties(self, mol: Chem.Mol) -> Tuple[Chem.Mol, bool]:
        """Compute molecular properties."""
        # TODO: either refactor to compute_charges, or make dynamic. currently acts mostly as a 'forcefield' curation step.
        try:
            Chem.rdPartialCharges.ComputeGasteigerCharges(
                mol, 
                nIter=20, 
                throwOnParamFailure=True
            )
            return mol, True
        except ValueError as e:
            self.log(f"Property computation error: {e}", level='DEBUG')
            return mol, False
    
    def _mol_to_smiles_step(self, mol: Union[Chem.Mol, List[Chem.Mol]], return_info: bool = False) -> Union[str, Tuple[str, bool, str]]:
        """Convert molecule(s) to SMILES string(s)."""
        canonical = 'canonical' in self.smiles_steps
        kekulize = 'kekulize' in self.smiles_steps
        if mol is None:
            curated_smiles = None
        elif isinstance(mol, list):
            # Handle stereoisomer enumeration
            curated_smiles = [
                Chem.MolToSmiles(
                    isomer,
                    canonical=canonical,
                    isomericSmiles=(self.stereo_policy != 'remove'),
                    kekuleSmiles=kekulize
                ) for isomer in mol
            ]
        else:
            curated_smiles = Chem.MolToSmiles(
                mol,
                canonical=canonical,
                isomericSmiles=(self.stereo_policy != 'remove'),
                kekuleSmiles=kekulize
            )
        
        if curated_smiles == '':
            curated_smiles = None
        
        success = curated_smiles is not None and len(curated_smiles) > 0
        
        if return_info:
            return curated_smiles, success, 'smiles'
        return curated_smiles
    
    # ==================== Multiprocessing Methods ====================
    
    def _curate_dataframe_multiprocess(self, df: pd.DataFrame, n_processes: Optional[int], 
                                       chunk_size: Optional[int], progress_interval: int) -> pd.DataFrame:
        """Process DataFrame using multiple processes."""
        
        if n_processes is None or chunk_size is None:
            calc_n_processes, calc_chunk_size, calc_n_chunks = calculate_chunk_params(len(df), self.MAX_CHUNK_SIZE)
            n_processes = n_processes or calc_n_processes
            chunk_size = chunk_size or calc_chunk_size
            n_chunks = (len(df) + chunk_size - 1) // chunk_size
        else:
            n_chunks = (len(df) + chunk_size - 1) // chunk_size
        
        self.log(
            f"Processing {len(df):,} molecules with {n_processes} processes "
            f"({n_chunks} chunks of {chunk_size:,})."
        )
        
        # Prepare chunks
        smiles_list = df[self.SMILES_column].tolist()
        smiles_chunks = [
            smiles_list[i:i + chunk_size] 
            for i in range(0, len(smiles_list), chunk_size)
        ]
        
        # Process chunks in parallel
        results = self._process_chunks_parallel(
            smiles_chunks, 
            n_processes, 
            progress_interval
        )
        
        # Combine results
        return self._combine_results(df, results)
    
    def _process_chunks_parallel(self, chunks: List[List[str]], n_processes: int, progress_interval: int) -> List[List[Tuple]]:
        """Process chunks in parallel with progress tracking."""
        total_chunks = len(chunks)
        chunk_results = [None] * total_chunks
        
        start_time = time.time()
        last_progress_time = start_time
        completed_chunks = 0
        
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            worker_func = partial(
                multiprocess_worker,
                actor_class=self.__class__,
                actor_method='_try_curate_smiles',
                actor_params=self._params,
                logger=self.logger,
                method_kwargs={'return_info': True},
            )
            
            future_to_chunk = {
                executor.submit(worker_func, chunk): i 
                for i, chunk in enumerate(chunks)
            }
            
            for future in as_completed(future_to_chunk, timeout=3600):
                try:
                    chunk_idx = future_to_chunk[future]
                    chunk_results[chunk_idx] = future.result(timeout=300)  # 5 min per chunk
                    completed_chunks += 1
                
                except Exception as e:
                    self.log(f"Failed to get result for chunk {chunk_idx}: {e}", level="ERROR")
                    chunk_results[chunk_idx] = [(None, False, 'retrieval_error')] * len(chunks[chunk_idx])
                    
                # Report progress
                if self._should_report_progress(
                    time.time(), last_progress_time, 
                    progress_interval, completed_chunks, total_chunks
                ):
                    self._report_progress(
                        completed_chunks, total_chunks,
                        start_time, len(chunks[0])
                    )
                    last_progress_time = time.time()
        
        # Flatten results
        all_results = []
        for chunk_result in chunk_results:
            all_results.extend(chunk_result)
        
        gc.collect()
        
        elapsed = time.time() - start_time
        rate = len(all_results) / elapsed if elapsed > 0 else 0
        self.log(
            f"Completed multiprocessing in {self._format_duration(elapsed)} | "
            f"Rate: {rate:.0f} mol/s"
        )
        
        return all_results
    
    def _curate_dataframe_single_process(self, df: pd.DataFrame) -> pd.DataFrame:
        """Process DataFrame using a single process."""
        curation_results = self._try_curate_smiles_array(
            df[self.SMILES_column], 
            return_info=True
        )
        
        return self._combine_results(df, curation_results)
    
    def _combine_results(self, df: pd.DataFrame, results: List[Tuple]) -> pd.DataFrame:
        """Combine curation results with original DataFrame."""
        results_df = pd.DataFrame(results, columns=self.output_cols)
        
        for col in self.output_cols:
            df[col] = results_df[col].values
        
        # Handle stereoisomer enumeration
        if (self.stereo_policy == 'enumerate' and 
            any(isinstance(smi, list) for smi, _, _ in results)):
            df = df.explode(self.output_cols[0])
            self.log(
                f"WARNING: Stereoisomers enumerated. "
                f"Dataset expanded: {len(results)} -> {len(df)} rows"
            )
        
        # Report success rate
        success_count = df['curation_success'].sum()
        self.log(
            f"Molecular curation: {success_count}/{len(df)} successful. "
            f"drop_na={self.dropna}"
        )
        
        if self.dropna:
            df = df[df.curation_success == True]
        
        return df
    
    # ==================== Post-processing Methods ====================
    
    def _post_process_duplicates(self, curated_dataset: pd.DataFrame) -> pd.DataFrame:
        """Handle duplicate SMILES after curation using ChEMBL actor method."""
        smiles_col = self.output_cols[0]

        # Check for duplicates
        unique_count = len(curated_dataset[smiles_col].unique())
        duplicate_count = curated_dataset.duplicated(smiles_col, keep=False).sum()

        if unique_count == len(curated_dataset) and duplicate_count == 0:
            self.log("No post-curation SMILES duplicates found.")
            return curated_dataset

        # Try to get ChEMBLCurator actor instance from context
        chembl_actor = self.get_actor(Steps.CHEMBL)

        if chembl_actor is not None:
            # Use ChEMBLCurator's handle_duplicate_smiles method directly
            self.log(f"Post-curation duplicates will be handled by ChEMBLCurator.handle_duplicate_smiles().")
            # TODO: keep most source-assigned stereo-isomer if duplication originates from stereo-assignment of racemic smiles.

            return chembl_actor.handle_duplicate_smiles(
                curated_dataset,
                SMILES_column=smiles_col
            )

        else:
            # Simple deduplication fallback
            duplicates_to_remove = curated_dataset.duplicated(
                smiles_col,
                keep=self.duplicates_policy
            ).sum()

            self.log(
                f"Found {duplicates_to_remove} duplicate SMILES. "
                f"Removing with policy: {self.duplicates_policy}",
                level='WARNING'
            )

            return curated_dataset.drop_duplicates(
                subset=[smiles_col],
                keep=self.duplicates_policy
            )
    
    # ==================== Helper Methods ====================
    
    def _initialize_components(self):
        """Initialize curation components based on configuration."""
        if 'handleStereo' in self.mol_steps:
            self._initialize_stereo_handler()

        if 'desalt' in self.mol_steps:
            self._initialize_salt_remover()

        if 'tautomers' in self.mol_steps:
            self._initialize_tautomer_enumerator()

        if 'neutralize' in self.mol_steps:
            self._initialize_neutralization_reactions()
    
    def _initialize_stereo_handler(self):
        """Initialize stereochemistry handler."""
        self.StereoHandler = StereoHandler(
            stereo_policy=self.stereo_policy,
            assign_policy=self.assign_policy,
            max_isomers=self.max_isomers,
            try_embedding=self.try_embedding,
            only_unassigned=self.only_unassigned,
            only_unique=self.only_unique,
            random_seed=self.random_seed,
            verbose=self.verbose,
            logger=self.logger,
            suppress_init_logs=self._suppress_init_logs,
        )

        if (self.stereo_policy == 'enumerate' and
            self.mol_steps[-1] != 'handleStereo'):
            self.log(
                "WARNING: Isomer enumeration enabled, but HandleStereo "
                "not at end of pipeline.",
                level='WARNING'
            )
    
    def _initialize_salt_remover(self):
        """Initialize salt remover."""
        self.salt_remover = SaltRemover()

        if (
            self.brute_force_desalt and 
            self.mol_steps[0] != 'desalt'):
            self.log(
                "brute_force_desalt=True, but desalt not at beginning of pipeline. "
                "May result in loss of metadata on MOL objects.",
                level='WARNING'
            )
    
    def _initialize_tautomer_enumerator(self):
        """Initialize tautomer enumerator."""
        cup = CleanupParameters()
        cup.maxTautomers = self.max_tautomers
        cup.maxTransforms = self.max_tautomers
        cup.tautomerReassignStereo = True
        cup.tautomerRemoveBondStereo = True
        cup.tautomerRemoveSp3Stereo = True
        
        self.tautomer_enumerator = TautomerEnumerator(cup)

        if self.max_tautomers > 256:
            self.log(
                f"Tautomer enumeration may slow down curation with "
                f"max_tautomers={self.max_tautomers}.",
                level='WARNING'
            )
    
    def _initialize_neutralization_reactions(self):
        """Initialize and compile neutralization reactions."""
        self.neutralization_reactions = []
        
        for smarts, smiles in self.NEUTRALIZATION_PATTERNS:
            try:
                reactant = Chem.MolFromSmarts(smarts)
                product = Chem.MolFromSmiles(smiles, False)
                
                if reactant and product:
                    self.neutralization_reactions.append((reactant, product))
            except Exception as e:
                self.log(
                    f"Failed to compile neutralization pattern {smarts}: {e}",
                    level='WARNING'
                )
    
    def _log_configuration(self):
        """Log configuration details."""
        self.log(f"CurateMol configured with MOL steps: {' > '.join(self.mol_steps)}")
        self.log(f"CurateMol configured with SMILES steps: {', '.join(self.smiles_steps)}")
    
    def _with_progress_tracking(self, func: Callable, items: List[Any], interval: Optional[int] = None, **kwargs) -> List[Any]:
        """Execute function on list with progress tracking."""
        interval = interval or self.DEFAULT_CHUNK_PROGRESS_INTERVAL
        total = len(items)
        start_time = time.time()
        results = []
        
        for i, item in enumerate(items):
            results.append(func(item, **kwargs))
            
            if (i + 1) % interval == 0 or (i + 1) == total:
                self._report_item_progress(i + 1, total, start_time)
        
        return results
    
    def _report_item_progress(self, completed: int, total: int, start_time: float):
        """Report progress for item processing."""
        elapsed = time.time() - start_time
        rate = completed / elapsed if elapsed > 0 else 0
        remaining_time = (total - completed) / rate if rate > 0 else 0
        
        self.log(
            f"Processed {completed:,}/{total:,} molecules "
            f"({100 * completed / total:.1f}%) | "
            f"Rate: {rate:.0f} mol/s | "
            f"ETA: {remaining_time / 60:.1f}m"
        )
    
    def _report_progress(self, completed_chunks: int, total_chunks: int, 
                         start_time: float, chunk_size: int):
        """Report multiprocessing progress."""
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
                f"Elapsed: {self._format_duration(elapsed)} | "
                f"ETA: {self._format_duration(eta) if eta > 0 else 'N/A'} | "
                f"Rate: {rate:.0f} mol/s"
            )
        else:
            self.log(
                f"Progress: {completed_chunks}/{total_chunks} chunks "
                f"({progress_pct:.1f}%) | "
                f"Elapsed: {self._format_duration(elapsed)}"
            )
    
    def _should_report_progress(self, current_time: float, last_report_time: float, 
                                interval: int, completed: int, total: int) -> bool:
        """Determine if progress should be reported."""
        time_exceeded = (current_time - last_report_time) >= interval
        is_complete = completed == total
        return time_exceeded or is_complete
    
    def _format_duration(self, seconds: float) -> str:
        """Format duration in human-readable format."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            return f"{seconds / 60:.1f}m"
        else:
            hours = seconds // 3600
            minutes = (seconds % 3600) / 60
            return f"{hours:.0f}h {minutes:.0f}m"
    
    def _get_smiles_safe(self, mol: Optional[Chem.Mol], default: str = 'Invalid Molecule') -> str:
        """Safely get SMILES string from molecule."""
        if mol is None:
            return default
        try:
            return Chem.MolToSmiles(mol)
        except Exception:
            return default
    
    def _handle_timeout(self, step: str, original_smiles: str,
        mol: Optional[Chem.Mol], error: TimeoutError):
        """Handle timeout errors during curation."""
        current_smiles = self._get_smiles_safe(mol)
        self.log(
            f"TIMEOUT: Step '{step}' timed out for molecule "
            f"(original: {original_smiles} | current: {current_smiles})",
            level="WARNING"
        )
        self.log(f"TIMEOUT ERROR: {error}", level="DEBUG")
    
    def _handle_curation_error(self, step: str, original_smiles: str,
        mol: Optional[Chem.Mol], error: Exception):
        """Handle general errors during curation."""
        current_smiles = self._get_smiles_safe(mol)
        self.log(
            f"Curation failed at step '{step}' for molecule "
            f"(original: {original_smiles} | current: {current_smiles})",
            level="WARNING"
        )
        self.log(f"ERROR: {error}", level="DEBUG")
    
    @contextmanager
    def _timeout_context(self, timeout_seconds: int, step_name: str):
        """
        Context manager for timing out curation steps.
        
        Useful for handling ultra-complex molecules that may cause
        RDKit operations to hang.
        """
        def timeout_handler(signum, frame):
            raise TimeoutError(
                f"Curation step '{step_name}' timed out after {timeout_seconds}s"
            )
        
        # Set up timeout signal
        old_handler = signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(timeout_seconds)
        
        try:
            yield
        finally:
            # Reset signal handling
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)