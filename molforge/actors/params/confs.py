"""
Unified conformer generation parameters.

Supports multiple backends (RDKit, OpenEye, etc.) with shared and backend-specific parameters.
"""

from typing import Literal, Optional, Dict, Any
from dataclasses import dataclass, field

import os
import shutil
import multiprocessing as mp

from .base import BaseParams


@dataclass
class GenerateConfsParams(BaseParams):
    """
    Unified conformer generation parameters.

    Provides a single parameter class that supports multiple backends,
    with both shared parameters (work across all backends) and
    backend-specific parameters (only used by relevant backend).
    """

    # ==================== Backend Selection ====================

    backend: Literal['rdkit', 'openeye'] = 'rdkit'
    """Conformer generation backend: 'rdkit' or 'openeye'"""

    # ==================== Shared Parameters ====================
    # These parameters work across all backends

    max_confs: int = 200
    """Maximum number of conformers to generate per molecule"""

    rms_threshold: float = 0.5
    """RMS threshold for conformer pruning (Angstroms)"""

    SMILES_column: str = 'curated_smiles'
    """DataFrame column containing SMILES strings"""

    names_column: str = 'molecule_chembl_id'
    """DataFrame column containing molecule identifiers"""

    dropna: bool = True
    """Drop rows with missing SMILES before processing"""

    convert_to_rdkit: bool = True
    """Convert conformers to RDKit Mol objects when extracting (OpenEye only, RDKit is native)"""

    # ==================== RDKit-Specific Parameters ====================
    # These are only used when backend='rdkit'

    use_random_coords: bool = True
    """Use random coordinates for initial embedding (RDKit)"""

    random_seed: int = 42
    """Random seed for reproducibility (RDKit)"""

    num_threads: int = 0
    """Number of threads for generation, 0 = all cores (RDKit)"""

    use_uff: bool = True
    """Use UFF force field for conformer optimization (RDKit)"""

    max_iterations: int = 200
    """Maximum optimization iterations per conformer (RDKit)"""

    # ==================== OpenEye-Specific Parameters ====================
    # These are only used when backend='openeye'

    mode: Literal['classic', 'macrocycle', 'rocs', 'pose', 'dense', 'fastrocs'] = 'classic'
    """OMEGA generation mode (OpenEye)"""

    use_gpu: bool = True
    """Enable GPU acceleration if available (OpenEye)"""

    mpi_np: int = -1
    """Number of MPI processes, -1 = auto-detect (OpenEye)"""

    strict: bool = True
    """Only process molecules with complete (stereochemistry) information (OpenEye)"""

    flipper: bool = False
    """Enable unassigned stereo-isomer flipping (OpenEye)"""

    flipper_warts: bool = False
    """Add suffix to flipped stereoisomers names (OpenEye)"""

    flipper_maxcenters: int = 4
    """Maximum stereocenters to enumerate with flipper (OpenEye)"""

    oeomega_path: Optional[str] = None
    """Path to OMEGA executable (OpenEye), None = auto-detect"""

    # ==================== Advanced ====================

    backend_kwargs: Dict[str, Any] = field(default_factory=dict)
    """Additional backend-specific keyword arguments"""

    def _validate_params(self) -> None:
        """Validate parameter values."""
        # Shared validation
        if self.max_confs <= 0:
            raise ValueError("max_confs must be positive")
        if self.rms_threshold <= 0:
            raise ValueError("rms_threshold must be positive")

        # RDKit-specific validation
        if self.backend == 'rdkit':
            if self.random_seed < 0:
                raise ValueError("random_seed must be non-negative")
            if self.max_iterations <= 0:
                raise ValueError("max_iterations must be positive")

        # OpenEye-specific validation
        if self.backend == 'openeye':
            if self.flipper and self.flipper_maxcenters <= 0:
                raise ValueError("flipper_maxcenters must be positive")

    def _post_init_hook(self) -> None:
        """Post-initialization setup."""
        # OpenEye: Auto-detect OMEGA executable if needed
        if self.backend == 'openeye' and self.oeomega_path is None:
            self.oeomega_path = self._find_omega_executable()

        # Auto-detect CPU count for MPI if needed
        if self.backend == 'openeye' and self.mpi_np == -1:
            self.mpi_np = max(1, mp.cpu_count() - 1)

    def _find_omega_executable(self) -> str:
        """
        Locate OMEGA executable automatically.

        Returns:
            Path to OMEGA executable

        Raises:
            RuntimeError: If OMEGA not found
        """
        possible_paths = [
            os.path.expanduser("~/OpenEye/openeye/bin/oeomega"),
            "/usr/local/bin/oeomega",
            "oeomega",  # Assume in PATH
        ]

        for path in possible_paths:
            if os.path.exists(path) or shutil.which(path):
                return path

        raise RuntimeError(
            "OMEGA executable not found. "
            "Please install OpenEye tools or set 'oeomega_path' parameter."
        )
