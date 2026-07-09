"""
Unified ChEMBL data source actor.

Supports multiple backends (SQL, API) with a single consistent interface.
The backend can be selected via the ChEMBLSourceParams.backend parameter.

Architecture:
    - Unified actor interface (ChEMBLSource)
    - Pluggable backend system (SQL, API)
    - Consistent output format across backends
    - Automatic backend registration via registry

Usage:
    >>> from molforge.actors.source import ChEMBLSource
    >>> from molforge.actors.params.source import ChEMBLSourceParams
    >>>
    >>> # SQL backend (default)
    >>> params = ChEMBLSourceParams(backend='sql', db_path='/path/to/chembl.db')
    >>> actor = ChEMBLSource(params, logger)
    >>> output = actor(input_data)
    >>>
    >>> # API backend
    >>> params = ChEMBLSourceParams(backend='api', n=5000)
    >>> actor = ChEMBLSource(params, logger)
    >>> output = actor(input_data)
"""

import pandas as pd
from typing import List

from .base import BaseActor
from .protocol import ActorOutput
from ..backends.registry import BackendRegistry
from ..configuration.steps import Steps
# Import source package to trigger backend registration
import molforge.backends.source  # noqa: F401


class ChEMBLSource(BaseActor):
    """ChEMBL data source actor."""
    __step_name__ = Steps.SOURCE
    __initial__   = True   # initial step: expected to be the first step in the pipeline
    """
    Unified ChEMBL data source actor supporting multiple backends.

    Fetches activity data for protein targets from ChEMBL using either
    a local SQLite database or the web API.

    Backend options:
        - 'sql': Queries local ChEMBL database
        - 'api': Fetches from ChEMBL web API with rate limiting

    Workflow:
        1. Get target_chembl_id from context (input_id)
        2. Fetch activity data via backend
        3. Return DataFrame with activity records

    Returns DataFrame with standardized ChEMBL activity columns.
    """

    def __post_init__(self):
        """Initialize ChEMBL data source backend."""
        # Get backend class from unified registry
        # Use step name ('source') as namespace
        backend_class = BackendRegistry.get('source', self.backend)

        # Instantiate backend with params, logger, and context
        self.backend_instance = backend_class(self._params, logger=self.logger, context=self._context)

        self.log(f"Initialized ChEMBL source with {self.backend} backend")

    @property
    def required_columns(self) -> List[str]:
        """No input columns required - this is a data source."""
        return []

    @property
    def output_columns(self) -> List[str]:
        """Columns returned by ChEMBL data fetch."""
        return [
            'activity_id', 'canonical_smiles', 'standard_value', 'standard_type',
            'standard_units', 'standard_relation', 'pchembl_value', 'molecule_chembl_id',
            'target_chembl_id', 'assay_type', 'document_year'
        ]

    @property
    def required_context_keys(self) -> List[str]:
        """Requires input_id (ChEMBL target ID) in context."""
        return ['input_id']

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Fetch activity data for a ChEMBL protein target.

        Target ID is retrieved from pipeline context via 'input_id' key.

        Args:
            data: Empty DataFrame (data source actor)

        Returns:
            DataFrame containing activity data

        Raises:
            ValueError: If target_chembl_id not in context
            Exception: Backend-specific errors
        """
        # Get target ID from context (set by pipeline)
        target_chembl_id = self.get_context('input_id')
        if not target_chembl_id:
            raise ValueError("target_chembl_id not found in context. Set 'input_id' in pipeline.")

        # Fetch activities via backend
        activities_df = self.backend_instance.fetch_activities(target_chembl_id)

        return activities_df

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with backend metadata."""
        target_id = self.get_context('input_id')

        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'source': self.backend_instance.source_description,
                'backend': self.backend,
                'target_chembl_id': target_id,
                'n_activities': len(data),
                'search_all': self.search_all
            }
        )
