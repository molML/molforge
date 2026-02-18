"""
Abstract base class for ChEMBL data source backends.

All ChEMBL data retrieval methods (SQL, API, etc.) must implement this interface.
"""

from abc import ABC, abstractmethod
from typing import Optional, TYPE_CHECKING
import pandas as pd

from ..registry import Backend

if TYPE_CHECKING:
    from ...actors.params.source import ChEMBLSourceParams
    from ...configuration.context import PipelineContext


class ChEMBLSourceBackend(Backend):
    """
    Abstract base class for ChEMBL data source backends.

    Provides a unified interface for different ChEMBL data retrieval methods,
    ensuring consistent behavior and making it easy to add new backends.
    """

    def __init__(self, params: 'ChEMBLSourceParams', logger=None, context: Optional['PipelineContext'] = None):
        """
        Initialize backend with parameters, logger, and context.

        Args:
            params: ChEMBL source parameters
            logger: Logger instance for consistent logging (optional)
            context: Pipeline context for accessing input_id and shared state (optional)
        """
        self.params = params
        self.logger = logger
        self.context = context

    @abstractmethod
    def fetch_activities(self, target_chembl_id: str) -> pd.DataFrame:
        """
        Fetch activity data for a ChEMBL protein target.

        Args:
            target_chembl_id: ChEMBL identifier for the target protein

        Returns:
            DataFrame containing activity data with standardized columns

        Raises:
            ValueError: If target_chembl_id is invalid
            Exception: Backend-specific errors (network, database, etc.)
        """
        pass

    @property
    @abstractmethod
    def source_description(self) -> str:
        """
        Description of data source for metadata.

        Returns:
            Source description (e.g., 'ChEMBL SQL', 'ChEMBL API')
        """
        pass

    @property
    def method(self) -> str:
        """
        Backend identifier for logging with consistent width.

        Returns:
            Formatted logging prefix (e.g., '[     SQL     ]', '[     API     ]')
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
