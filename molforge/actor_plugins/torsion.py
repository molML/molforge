"""
Torsion analysis actor plugin.

Analyzes torsional properties for molecules with conformers.
Requires phd-tools package for calculations.
"""

from typing import List, Tuple, Dict
import traceback

import pandas as pd

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

from dataclasses import dataclass
from typing import Optional, Literal


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

    device: str = 'cpu'
    account_for_symmetry: bool = True
    symmetry_radius: int = 3
    ignore_colinear_bonds: bool = True
    n_confs_threshold: int = 50
    numerical_epsilon: float = 1e-10
    aggregation_method: Optional[Literal['min', 'max', 'mean']] = None

    def _validate_params(self) -> None:
        """Validate parameter values."""
        if self.device != 'cpu':
            raise ValueError("CPU device required for torsion calculations")

        if self.aggregation_method is not None:
            self._validate_policy('aggregation_method', self.aggregation_method,
                                ['min', 'max', 'mean'])

        if self.n_confs_threshold < 1:
            raise ValueError("Conformer threshold must be at least 1")


# ============================================================================
# ACTOR
# ============================================================================

class CalculateTorsions(BaseActor):
    """
    Torsion analysis actor.

    Analyzes torsional properties for molecules with multiple conformers.
    Outputs bond-level metrics and conformer statistics.

    Requires:
        - GenerateConfs actor output (conformers)
        - phd-tools package (private)

    Output columns:
        - torsion_mapping: Canonical bond flexibility mapping {(rank1, rank2): variance}
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
        'torsion_mapping', 'torsion_success', 'n_confs', 'low_confs',
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

        self.log(f"Torsion analysis initialized (device: {self.device})")

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Analyze torsional properties for molecules with conformers.

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
        # Note: extract_molecules() only yields successful molecules
        mol_generator = gc_actor.extract_molecules()

        # Process molecules and collect results
        # Pass DataFrame success mask so we know when to pull from generator
        torsion_stats = self.process_molecules(mol_generator, df['conformer_success'])

        # Assign results to DataFrame
        for col, values in torsion_stats.items():
            df[col] = values

        return df

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with torsion mapping endpoint."""
        if 'n_torsions' in data.columns and 'torsion_success' in data.columns:
            success_df = data[data['torsion_success'] == True]
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
            },
            endpoint=self.forge_endpoint
        )

    def process_molecules(self, molecules, success_mask: pd.Series) -> Dict[str, List]:
        """
        Process molecules from generator and return column data.

        Args:
            molecules: Generator yielding successful molecules (RDKit or OpenEye)
            success_mask: Boolean series indicating which rows had successful conformers

        Returns:
            Dictionary mapping column names to lists of values
        """
        # Initialize result collectors
        results = {col: [] for col in self.OUTPUT_COLUMNS}

        # Iterate through success mask and pull from generator only for successful rows
        for row_idx, has_conformers in enumerate(success_mask):
            if has_conformers:
                # Pull next successful molecule from generator
                try:
                    mol = next(molecules)
                    torsion_mapping, stats, success = self._process_single_molecule(mol, row_idx)
                except StopIteration:
                    # Generator exhausted unexpectedly
                    self.log(f"Generator exhausted at row {row_idx}", level='ERROR')
                    torsion_mapping, stats, success = self._get_empty_result()
            else:
                # No conformer for this row - use empty results
                torsion_mapping, stats, success = self._get_empty_result()

            # Store results 
            results['torsion_mapping'].append(torsion_mapping.__repr__()) # Store as string
            results['torsion_success'].append(success)
            results['n_confs'].append(stats['n_confs'])
            results['low_confs'].append(stats['low_confs'])
            results['n_ring_torsions'].append(stats['n_ring_torsions'])
            results['n_rotor_torsions'].append(stats['n_rotor_torsions'])
            results['n_torsions'].append(stats['n_torsions'])
            results['mean_variance'].append(stats['mean_variance'])
            results['warnings'].append(stats['warnings'])

        return results

    def _process_single_molecule(self, mol, mol_idx: int = -1) -> Tuple[Dict, Dict, bool]:
        """
        Process single molecule and return canonical results.

        Args:
            mol: RDKit or OpenEye molecule with conformers
            mol_idx: Molecule index for logging

        Returns:
            Tuple of (canonical_torsion_mapping, statistics, success)
        """
        if mol is None:
            return self._get_empty_result()

        try:
            # Step 1: Calculate torsions (returns bond_idx -> variance)
            bond_variance, stats = TorsionCalculator.compute_flexibility(
                mol=mol,
                account_for_symmetry=self.account_for_symmetry,
                symmetry_radius=self.symmetry_radius,
                ignore_colinear_bonds=self.ignore_colinear_bonds,
                numerical_epsilon=self.numerical_epsilon,
                aggregation_method=self.aggregation_method,
                n_confs_threshold=self.n_confs_threshold,
                device=self.device
            )

            # Step 2: Get RDKit topology for canonicalization
            if HAS_OPENEYE and isinstance(mol, oechem.OEMol):
                rd_topology = MoleculeConverter.openeye_to_rdkit_topology(mol)
            else:
                rd_topology = mol

            # Step 3: Canonicalize bond indices to atom rank pairs
            canonical_mapping = CanonicalBondProperties.bond_idx_to_canonical(
                rd_topology, bond_variance
            )

            return canonical_mapping, stats, True

        except Exception as e:
            error_msg = f"{type(e).__name__}: {str(e)}" if str(e) else f"{type(e).__name__}"
            self.log(f"Error processing molecule {mol_idx}: {error_msg}", level='WARNING')
            self.log(f"Traceback: {traceback.format_exc()}", level='DEBUG')

            return self._get_empty_result()

    @staticmethod
    def _get_empty_result() -> Tuple[Dict, Dict, bool]:
        """
        Return empty torsion results for failed molecules.

        Returns:
            Tuple of (empty_mapping, empty_stats, success=False)
        """
        empty_mapping = {}
        empty_stats = {
            'n_confs': 0,
            'low_confs': False,
            'n_ring_torsions': 0,
            'n_rotor_torsions': 0,
            'n_torsions': 0,
            'mean_variance': 0.0,
            'warnings': []
        }
        return empty_mapping, empty_stats, False
