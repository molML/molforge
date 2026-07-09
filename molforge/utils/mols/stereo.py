"""
Stereochemistry handling utility for molecular curation.

Provides stereoisomer enumeration, assignment, and removal capabilities.
"""

from typing import List, Optional, Tuple, Any
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import RemoveStereochemistry
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import FindMolChiralCenters


class StereoHandler:
    """
    Helper class for stereochemistry curation.

    Supports multiple policies:
    - keep: Preserve existing stereochemistry
    - remove: Strip all stereochemistry
    - assign: Enumerate and select a single isomer
    - enumerate: Return all possible stereoisomers
    """

    def __init__(
        self,
        stereo_policy: str = "keep",
        assign_policy: str = "first",
        max_isomers: int = 32,
        try_embedding: bool = False,
        only_unassigned: bool = True,
        only_unique: bool = True,
        random_seed: int = 42,
        verbose: bool = False,
        logger: Optional[Any] = None,
        suppress_init_logs: bool = False
    ):
        """
        Initialize stereochemistry handler.

        Args:
            stereo_policy: How to handle stereocenters (keep/remove/assign/enumerate)
            assign_policy: Selection method when assigning (first/random/lowest)
            max_isomers: Maximum stereoisomers to enumerate
            try_embedding: Try 3D embedding for enumeration
            only_unassigned: Only enumerate unassigned centers
            only_unique: Filter duplicate isomers
            random_seed: Random seed for reproducibility
            verbose: Enable verbose logging
            logger: Optional logging function (message, level)
            suppress_init_logs: If True, suppress logging during __init__.
                               Useful for multiprocessing workers to reduce noise.
        """
        # Initialize tracking attributes FIRST (needed for log suppression)
        self._suppress_init_logs = suppress_init_logs
        self._in_init_phase = True

        self.policy = stereo_policy
        self.assign = assign_policy
        self.max_isomers = max_isomers
        self.try_embedding = try_embedding
        self.only_unassigned = only_unassigned
        self.only_unique = only_unique
        self.seed = random_seed
        self.verbose = verbose
        self.logger = logger

        np.random.seed(random_seed)

        self.opts = StereoEnumerationOptions(
            tryEmbedding=try_embedding,
            maxIsomers=max_isomers,
            onlyUnassigned=only_unassigned,
            unique=only_unique,
            rand=random_seed
        )

        # Log configuration only if verbose
        if verbose and logger:
            self.log(f"StereoHandler: policy='{stereo_policy}', assign='{assign_policy}'", "INFO")

        # Mark initialization complete
        self._in_init_phase = False

    @property
    def _nametag(self) -> str:
        from ...configuration.steps import Steps
        width = max(Steps.max_length(), 7)  # min 7 for backend tags
        step_name = "stereo"
        return f"[{step_name.upper():^{width}s}]"
    
    # ==================== Public API ====================
            
    def curate(self, mol: Chem.Mol) -> Tuple[List[Chem.Mol], bool]:
        """
        Main entry point for stereochemistry curation.

        Args:
            mol: Input RDKit molecule

        Returns:
            Tuple of (list of curated molecules, success flag)
        """
        chiral_centers = self._find_chiral_centers(mol)
        if chiral_centers is None:
            return [mol], False

        if self.policy == "keep" or len(chiral_centers) == 0:
            return [mol], True

        if self.policy == "remove":
            curated = self._remove_stereochemistry(mol)
            return [curated], curated is not None

        isomers = self._enumerate_isomers(mol)
        if not isomers:
            self.log(f"Enumeration failed for {Chem.MolToSmiles(mol)}", "WARNING")
            return [mol], False

        if self.policy == "enumerate":
            return isomers, True

        if self.policy == "assign":
            selected = self._select_isomer(isomers)
            return [selected], selected is not None

        return [mol], False

    def __call__(self, mol: Chem.Mol) -> Tuple[List[Chem.Mol], bool]:
        """Allow instance to be called directly."""
        return self.curate(mol)

    # ==================== Core Methods ====================

    def _find_chiral_centers(self, mol: Chem.Mol) -> Optional[List[Tuple]]:
        """Detect chiral centers. Returns None on failure."""
        try:
            return FindMolChiralCenters(mol, includeUnassigned=True)
        except Exception as e:
            self.log(f"Chiral center detection failed: {e}", "ERROR")
            return None

    def _enumerate_isomers(self, mol: Chem.Mol) -> List[Chem.Mol]:
        """Enumerate stereoisomers using RDKit."""
        try:
            return list(EnumerateStereoisomers(mol, self.opts))
        except Exception as e:
            self.log(f"Enumeration error: {e}", "ERROR")
            return []

    def _select_isomer(self, isomers: List[Chem.Mol]) -> Optional[Chem.Mol]:
        """Select single isomer based on assignment policy."""
        if not isomers:
            return None

        if len(isomers) == 1 or self.assign == "first":
            return isomers[0]

        if self.assign == "random":
            return isomers[np.random.randint(len(isomers))]

        if self.assign == "lowest":
            return self._select_lowest_energy(isomers)

        return isomers[0]

    def _remove_stereochemistry(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """Remove all stereochemistry from molecule (in-place)."""
        try:
            RemoveStereochemistry(mol)
            return mol
        except Exception as e:
            self.log(f"Failed to remove stereochemistry: {e}", "ERROR")
            return None

    # ==================== Energy Ranking ====================

    def _select_lowest_energy(self, isomers: List[Chem.Mol]) -> Chem.Mol:
        """Select stereoisomer with lowest MMFF94 (or UFF) energy."""
        scored = []

        for iso in isomers:
            energy = self._calculate_energy(iso)
            if energy is not None:
                scored.append((energy, iso))

        if scored:
            return min(scored, key=lambda x: x[0])[1]
        else:
            self.log("Energy calculation failed for all isomers", "WARNING")
            return isomers[0]

    def _calculate_energy(self, mol: Chem.Mol) -> Optional[float]:
        """Calculate conformer energy (MMFF94 or UFF). Returns None on failure."""
        try:
            mol_h = Chem.AddHs(mol)

            if AllChem.EmbedMolecule(mol_h, randomSeed=self.seed) != 0:
                return None

            if AllChem.MMFFHasAllMoleculeParams(mol_h):
                props = AllChem.MMFFGetMoleculeProperties(mol_h)
                ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
                ff.Minimize()
                return ff.CalcEnergy()
            else:
                AllChem.UFFOptimizeMolecule(mol_h)
                ff = AllChem.UFFGetMoleculeForceField(mol_h)
                return ff.CalcEnergy()

        except Exception:
            return None

    # ==================== Logging ====================

    def log(self, message: str, level: str = 'INFO') -> None:
        """
        Standard logging interface.

        Args:
            message: Log message
            level: Log level (INFO, WARNING, ERROR, DEBUG)
        """
        # Suppress logs during init phase if flag is set
        if self._in_init_phase and self._suppress_init_logs:
            return

        if self.logger:
            self.logger.log(level.lower(), self._nametag, message)
        elif getattr(self, 'verbose', True):
            print(f"{self._nametag} | {level} | {message}")
