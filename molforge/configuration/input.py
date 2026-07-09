"""
Input parser for pipeline and standalone actor execution.

Separates input validation logic from pipeline execution.
"""

from typing import Union, Optional, List, Tuple, TYPE_CHECKING
from pathlib import Path
import pandas as pd

if TYPE_CHECKING:
    from ..actors.base import BaseActor


class InputParser:
    """
    Parse and validate inputs for actor execution.

    Handles three input patterns:
    1. ChEMBL ID string - for data source actors
    2. CSV file path - loads DataFrame
    3. DataFrame - direct data input

    Validates input matches first actor's requirements via required_context_keys.
    """

    @staticmethod
    def parse(
        actors: List['BaseActor'],
        input_data: Union[str, pd.DataFrame],
        input_id: Optional[str] = None
    ) -> Tuple[pd.DataFrame, str, str, str]:
        """
        Parse input for actor execution.

        Args:
            actors: List of actors (checks first for requirements)
            input_data: ChEMBL ID string, CSV file path, or DataFrame
            input_id: Optional identifier (required for DataFrame input)

        Returns:
            Tuple of (initial_data, input_id, input_type, input_source):
                - initial_data: DataFrame to pass to first actor (empty for data sources)
                - input_id: Identifier for this input
                - input_type: Type of input ('ChEMBL ID', 'CSV File', 'DataFrame')
                - input_source: Description of input source

        Raises:
            ValueError: If input doesn't match actor requirements
            FileNotFoundError: If CSV file not found
        """
        first_actor = actors[0]
        required_keys = first_actor.required_context_keys

        # Determine what first actor needs
        needs_target_id = 'input_id' in required_keys

        # Parse based on input type
        if isinstance(input_data, str):
            initial_data, metadata = InputParser._parse_string(input_data, input_id, needs_target_id)
        elif isinstance(input_data, pd.DataFrame):
            initial_data, metadata = InputParser._parse_dataframe(input_data, input_id, needs_target_id)
        else:
            raise ValueError(
                f"input_data must be str or DataFrame, got {type(input_data)}"
            )

        return initial_data, metadata['input_id'], metadata['input_type'], metadata['input_source']

    @staticmethod
    def _parse_string(
        input_data: str,
        input_id: Optional[str],
        needs_target_id: bool
    ) -> Tuple[pd.DataFrame, dict]:
        """Parse string input (ChEMBL ID or CSV path)."""

        # CSV file
        if input_data.endswith('.csv'):
            if needs_target_id:
                raise ValueError(
                    "First actor requires target ID (e.g., ChEMBL ID), but CSV file provided. "
                    "Use ChEMBL ID string instead."
                )
            return InputParser._handle_csv(input_data, input_id)

        # ChEMBL ID
        elif input_data.startswith('CHEMBL'):
            if not needs_target_id:
                raise ValueError(
                    "First actor expects DataFrame input, but ChEMBL ID provided. "
                    "Use CSV file or DataFrame instead."
                )
            return InputParser._handle_chembl_id(input_data, input_id)

        else:
            raise ValueError(
                f"String input must be ChEMBL ID (starts with 'CHEMBL') "
                f"or CSV path (ends with '.csv'), got: {input_data}"
            )

    @staticmethod
    def _parse_dataframe(
        df: pd.DataFrame,
        input_id: Optional[str],
        needs_target_id: bool
    ) -> Tuple[pd.DataFrame, dict]:
        """Parse DataFrame input."""

        if needs_target_id:
            raise ValueError(
                "First actor requires target ID (e.g., ChEMBL ID), but DataFrame provided. "
                "Use ChEMBL ID string instead."
            )

        if input_id is None:
            raise ValueError("input_id required when passing DataFrame")

        context_dict = {
            'input_id': input_id,
            'input_type': 'DataFrame',
            'input_source': f"DataFrame '{input_id}' ({df.shape[0]}x{df.shape[1]})"
        }

        return df, context_dict

    @staticmethod
    def _handle_csv(csv_path: str, input_id: Optional[str]) -> Tuple[pd.DataFrame, dict]:
        """Handle CSV file input."""
        path = Path(csv_path)

        if not path.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        if not path.suffix.lower() == '.csv':
            raise ValueError(f"File must have .csv extension: {csv_path}")

        try:
            df = pd.read_csv(path)
        except Exception as e:
            raise ValueError(f"Error reading CSV file: {e}")

        context_dict = {
            'input_id': input_id or path.stem,  # Use filename if no input_id provided
            'input_type': 'CSV File',
            'input_source': str(path.absolute())
        }

        return df, context_dict

    @staticmethod
    def _handle_chembl_id(chembl_id: str, input_id: Optional[str]) -> Tuple[pd.DataFrame, dict]:
        """Handle ChEMBL ID input."""
        context_dict = {
            'input_id': chembl_id,  # ChEMBL ID is the input_id
            'input_type': 'ChEMBL ID',
            'input_source': f"ChEMBL ID: {chembl_id}"
        }

        # Return empty DataFrame - data source actor fetches using input_id from context
        return pd.DataFrame(), context_dict
