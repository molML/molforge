"""
Abstract base class for conformer generation backends.

All conformer generation methods (RDKit, OpenEye, etc.) must implement this interface.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Iterator, TYPE_CHECKING
from rdkit import Chem
import pandas as pd

from ..registry import Backend

if TYPE_CHECKING:
    from ...actors.params.confs import GenerateConfsParams
    from ...actors.protocol import PipelineContext


class ConformerBackend(Backend):
    """
    Abstract base class for conformer generation backends.

    Provides a unified interface for different conformer generation methods,
    ensuring consistent behavior and making it easy to add new backends.
    """

    def __init__(self, params: 'GenerateConfsParams', logger=None, context: Optional['PipelineContext'] = None):
        """
        Initialize backend with parameters, logger, and context.

        Args:
            params: Conformer generation parameters
            logger: Logger instance for consistent logging (optional)
            context: Pipeline context for file paths and shared state (optional)
        """
        self.params = params
        self.logger = logger
        self.context = context

    @abstractmethod
    def generate_conformers(
        self,
        smiles_list: list[str],
        names_list: list[str]
    ) -> Dict[str, Dict[str, Any]]:
        """
        Generate conformers for SMILES.

        Args:
            smiles_list: List of SMILES strings
            names_list: List of molecule identifiers

        Returns:
            Dict mapping names to results. Each result contains:
                - 'success': bool - Whether generation succeeded
                - 'n_conformers': int - Number of conformers generated
                - 'status': str - Status message (empty for success, error description for failure)
        """
        pass

    @abstractmethod
    def get_endpoint(self) -> Any:
        """
        Return endpoint for downstream actors.

        The endpoint type depends on backend implementation:
        - File path (str) - For file-based backends (OpenEye)
        - Molecule generator (Iterator[Chem.Mol]) - For in-memory backends (RDKit)
        - MolBox object - For specialized storage

        Returns:
            Backend-specific endpoint
        """
        pass

    @abstractmethod
    def extract_molecules(self) -> Iterator[Chem.Mol]:
        """
        Extract all generated molecules in order.

        Yields:
            RDKit molecules with conformers attached
        """
        pass

    @abstractmethod
    def get_report_dataframe(self) -> pd.DataFrame:
        """
        Get full conformer generation report as DataFrame.

        Returns:
            DataFrame with columns matching OMEGA format:
                - Molecule: SMILES string
                - Title: Molecule name/identifier
                - Rotors: Rotatable bond count
                - Conformers: Number of conformers generated
                - ElapsedTime(s): Generation time
                - Status: Status message (empty string for success)
        """
        pass

    @abstractmethod
    def get_successful_names(self) -> list[str]:
        """
        Get list of successfully generated molecule names.

        Returns:
            List of molecule names that have conformers
        """
        pass

    @property
    def method(self) -> str:
        """
        Backend identifier for logging with consistent width.

        Returns:
            Formatted logging prefix (e.g., '[RDKIT-ETKDG]', '[   OEMEGA   ]')
        """
        width = 13
        backend_name = getattr(self.__class__, '__backend_name__', 'BACKEND')
        return f"[{backend_name.upper():^{width}s}]"

    def log(self, message: str, level: str = 'INFO'):
        """
        Helper for consistent logging across backends.

        Args:
            message: Log message
            level: Log level (INFO, WARNING, ERROR, DEBUG)
        """
        if self.logger:
            self.logger.log(level.lower(), self.method, message)

    def __repr__(self) -> str:
        """String representation."""
        backend_name = getattr(self.__class__, '__backend_name__', 'BACKEND')
        return f"{self.__class__.__name__}(backend={backend_name})"
