"""
Torsion analysis actor plugin.

Analyzes torsional properties for molecules with conformers.
Requires phd-tools package for calculations.
"""

from typing import List, Tuple, Dict, Iterator
import traceback

import pandas as pd
from rdkit import Chem

# Attempt to import private toolkit
try:
    from phd_tools.chemistry.torsion import TorsionCalculator
    from phd_tools.chemistry.savebonds import CanonicalBondPropertiesDB
    TORSION_AVAILABLE = True
except ImportError:
    TORSION_AVAILABLE = False
    TorsionCalculator = None
    CanonicalBondPropertiesDB = None


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
        'n_confs', 'low_confs', 'n_ring_torsions',
        'n_rotor_torsions', 'n_torsions', 'mean_variance', 'warnings'
    ]

    @property
    def required_columns(self) -> List[str]:
        """Required columns from GenerateConfs."""
        return ['conformer_success']

    @property
    def output_columns(self) -> List[str]:
        """Columns added by torsion analysis."""
        return self.OUTPUT_COLUMNS

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
        self.path_prefix = None
        self.db_path = None

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

        self.db_path = self._get_run_path("torsions.db")
        self.log(f'Results database: {self.db_path}')

        mol_generator = gc_actor.extract_molecules()

        _, torsion_stats = self.process_molecules(
            mol_generator,
            return_info=True,
            db_path=self.db_path
        )

        for col, values in torsion_stats.items():
            df[col] = values

        return df

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with database endpoint."""
        if 'n_torsions' in data.columns:
            total_torsions = data['n_torsions'].sum()
            mean_metric = data['mean_variance'].mean() if 'mean_variance' in data.columns else 0.0
        else:
            total_torsions = 0
            mean_metric = 0.0

        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'n_molecules': len(data),
                'total_torsions': int(total_torsions),
                'mean_metric': float(mean_metric),
                'database_path': self.db_path
            },
            endpoint=self.db_path
        )

    def read_results(self) -> Tuple[List[Chem.Mol], List[Dict[int, float]]]:
        """
        Read analysis results from database.

        Returns:
            Tuple of (molecules, bond_metrics)
        """
        if self.db_path is None:
            self.log("Database path not set. Run process() first.", level='ERROR')
            return [], []

        molecules, bond_metrics = CanonicalBondPropertiesDB.load_all_molecules(self.db_path)
        self.log(f"Loaded {len(molecules)} molecules from {self.db_path}")

        return molecules, bond_metrics

    def process_molecules(
        self,
        molecules: Iterator,
        return_info: bool = False,
        save_results: bool = True,
        db_path: str = 'default.db'
    ) -> Tuple[List[Dict], List[Dict]] | List[Dict]:
        """
        Process multiple molecules.

        Args:
            molecules: Iterator of RDKit molecules with conformers
            return_info: Return statistics if True
            save_results: Save to database if True
            db_path: Database file path

        Returns:
            Bond metrics for each molecule, optionally with statistics
        """
        bond_metrics = []

        if return_info:
            stats = {
                'n_confs': [], 'low_confs': [],
                'n_ring_torsions': [], 'n_rotor_torsions': [],
                'n_torsions': [], 'mean_variance': [],
                'warnings': [],
            }

        if save_results:
            self.log(f"Saving results to database: {db_path}")
            CanonicalBondPropertiesDB.initialize_database(db_path)
            saved = 0

        for mol_idx, mol in enumerate(molecules):
            if mol is None:
                if return_info:
                    for key in stats.keys():
                        stats[key].append(None)
                bond_metrics.append({})
                continue

            try:
                # Delegate to private calculator
                bond_data, mol_stats = TorsionCalculator.compute_flexibility(
                    mol=mol,
                    account_for_symmetry=self.account_for_symmetry,
                    symmetry_radius=self.symmetry_radius,
                    ignore_colinear_bonds=self.ignore_colinear_bonds,
                    numerical_epsilon=self.numerical_epsilon,
                    aggregation_method=self.aggregation_method,
                    n_confs_threshold=self.n_confs_threshold,
                    device=self.device
                )

                if return_info:
                    for key, val in mol_stats.items():
                        stats[key].append(val)

                bond_metrics.append(bond_data)

                if save_results:
                    mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{mol_idx}"
                    CanonicalBondPropertiesDB.add_molecule(
                        db_path, mol, bond_data, mol_name, overwrite=False
                    )
                    saved += 1

            except Exception as e:
                error_msg = f"{type(e).__name__}: {str(e)}" if str(e) else f"{type(e).__name__}"
                self.log(f"Error processing molecule {mol_idx}: {error_msg}", level='WARNING')
                self.log(f"Traceback: {traceback.format_exc()}", level='DEBUG')

                if return_info:
                    for key in stats.keys():
                        stats[key].append(None)
                bond_metrics.append({})

        if save_results:
            self.log(f"Saved {saved} entries to database")

        if return_info:
            return bond_metrics, stats
        else:
            return bond_metrics

    def process_molecule(self, mol: Chem.Mol, mol_idx: int = -1,
                        return_info: bool = False):
        """
        Process single molecule.

        Args:
            mol: RDKit molecule with conformers
            mol_idx: Molecule index (for compatibility)
            return_info: Return statistics if True

        Returns:
            Bond metrics dict, or (bond_metrics, stats, mol) if return_info=True
        """
        bond_data, stats = TorsionCalculator.compute_flexibility(
            mol=mol,
            account_for_symmetry=self.account_for_symmetry,
            symmetry_radius=self.symmetry_radius,
            ignore_colinear_bonds=self.ignore_colinear_bonds,
            numerical_epsilon=self.numerical_epsilon,
            aggregation_method=self.aggregation_method,
            n_confs_threshold=self.n_confs_threshold,
            device=self.device
        )

        if return_info:
            return bond_data, stats, mol
        else:
            return bond_data
