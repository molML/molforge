"""
Unified conformer generation actor.

Supports multiple backends (RDKit, OpenEye) with a single consistent interface.
The backend can be selected via the GenerateConfsParams.backend parameter.

Architecture:
    - Unified actor interface (GenerateConfs)
    - Pluggable backend system (RDKit, OpenEye)
    - Consistent output format across backends
    - Automatic backend registration via registry

Usage:
    >>> from molforge.actors.confs_unified import GenerateConfs
    >>> from molforge.actors.params.confs_unified import GenerateConfsParams
    >>>
    >>> # RDKit backend (default)
    >>> params = GenerateConfsParams(backend='rdkit', max_confs=200)
    >>> actor = GenerateConfs(params, context)
    >>> output = actor(input_data)
    >>>
    >>> # OpenEye backend (requires license)
    >>> params = GenerateConfsParams(backend='openeye', mode='macrocycle')
    >>> actor = GenerateConfs(params, context)
    >>> output = actor(input_data)
"""

import pandas as pd

from .base import BaseActor
from .protocol import ActorOutput
from ..backends.registry import BackendRegistry
from ..configuration.steps import Steps
# Import confs package to trigger backend registration
import molforge.backends.confs  # noqa: F401


class GenerateConfs(BaseActor):
    """Conformer generation actor."""
    __step_name__ = Steps.CONFS
    """
    Unified conformer generation actor supporting multiple backends.

    Generates 3D conformers using either RDKit ETKDG or OpenEye OMEGA.

    Backend options:
        - 'rdkit': Uses RDKit ETKDGv3 with parallel processing
        - 'openeye': Uses OMEGA executable (requires OpenEye license)

    Workflow:
        1. Extract SMILES and names from DataFrame
        2. Generate conformers via backend
        3. Merge backend report into DataFrame
        4. Add conformer_success column

    Returns DataFrame with conformer statistics and success flags.
    """

    def __post_init__(self):
        """Initialize conformer generation backend."""
        # Get backend class from unified registry
        # Use step name ('confs') as namespace
        backend_class = BackendRegistry.get('confs', self.backend)

        # Instantiate backend with params, logger, and context
        self.backend_instance = backend_class(self._params, logger=self.logger, context=self._context)

        self.log(f"Initialized conformer generation with {self.backend} backend")

    @property
    def required_columns(self) -> list[str]:
        """Columns required in input DataFrame."""
        return [self.SMILES_column]

    @property
    def output_columns(self) -> list[str]:
        """
        Columns added to output DataFrame.

        Backend report columns (merged from generation, matches OMEGA format):
            - Molecule: SMILES string used for conformer generation
            - Title: Molecule identifier/name used during generation
            - Rotors: Number of rotatable bonds detected
            - n_conformers: Number of conformers successfully generated (renamed from Conformers)
            - ElapsedTime(s): Generation time in seconds
            - Status: Backend-specific status message (empty for success)

        Actor-added columns:
            - conformer_success: Boolean flag indicating successful generation
        """
        return ['Molecule', 'Title', 'Rotors', 'n_conformers', 'ElapsedTime(s)', 'Status', 'conformer_success']

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Generate conformers for all SMILES in DataFrame.

        Args:
            data: Input DataFrame with SMILES and names

        Returns:
            DataFrame with conformer statistics added
        """
        df = data.copy()

        # Title column is the name used during conformer generation (matches OMEGA format)
        title_col = self.output_columns[1]
        if title_col in df.columns:
            raise ValueError(
                f"DataFrame contains '{title_col}' column which will conflict with conformer report. "
                f"Please rename this column before conformer generation."
            )

        # Setup Title column for mapping results
        # Handle cases where names_column doesn't exist or has duplicates
        if (self.names_column not in df.columns) or \
           (len(df[self.names_column].unique()) != len(df)):
            self.log(
                f"Names column '{self.names_column}' not in DataFrame or contains duplicate names. "
                f"Using enumeration for naming.",
                level='WARNING'
            )
            df[title_col] = [f'{self.backend.upper()}_{i}' for i in range(len(df))]
        else:
            df[title_col] = df[self.names_column]

        # Extract SMILES and names as lists
        smiles_list = df[self.SMILES_column].tolist()
        names_list = df[title_col].tolist()

        self.log(f"Generating conformers for {len(df)} SMILES")

        # Generate conformers via backend
        self.backend_instance.generate_conformers(smiles_list, names_list)

        # Get backend report and successful names
        report_df = self.backend_instance.get_report_dataframe()
        success_names = self.backend_instance.get_successful_names()

        # Merge report into DataFrame on Title column
        df = pd.merge(
            df,
            report_df,
            on=title_col,
            how='left',
            suffixes=('', f' ({self.backend.upper()})')
        )

        # Mark conformer success based on presence in successful names
        df['conformer_success'] = df[title_col].isin(success_names)

        # Rename 'Conformers' column to 'n_conformers' for consistency
        if 'Conformers' in df.columns:
            df = df.rename(columns={'Conformers': 'n_conformers'})

        # Ensure n_conformers has a value (0 for missing)
        if 'n_conformers' not in df.columns:
            df['n_conformers'] = 0
        df['n_conformers'] = df['n_conformers'].fillna(0).astype(int)

        # Calculate statistics
        success_count = df['conformer_success'].sum()
        failure_count = len(df) - success_count
        total_conformers = df['n_conformers'].sum()
        avg_conformers = total_conformers / success_count if success_count > 0 else 0

        # Log summary
        self.log(
            f"Conformer generation complete: "
            f"{success_count} succeeded, {failure_count} failed, "
            f"avg {avg_conformers:.1f} conformers/molecule"
        )

        # Log failure summary if any (matches original implementation)
        if failure_count > 0 and 'Status' in df.columns:
            failure_df = df[~df['conformer_success']]
            if len(failure_df) > 0:
                failure_summary = failure_df.groupby('Status').size()
                self.log(
                    f"Failure summary:\n{failure_summary.to_frame('count').to_markdown()}",
                    level='WARNING'
                )

        # Drop failed conformer generations if requested
        if self.dropna:
            initial_count = len(df)
            df = df[df['conformer_success']]
            dropped = initial_count - len(df)
            if dropped > 0:
                self.log(f"Dropped {dropped} rows with failed conformer generation")

        return df

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with backend endpoint and metadata."""
        # Calculate statistics
        success_count = data['conformer_success'].sum()
        total_conformers = data['n_conformers'].sum()
        avg_conformers = total_conformers / success_count if success_count > 0 else 0

        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'backend': self.backend,
                'molecules_processed': len(data),
                'success_count': int(success_count),
                'total_conformers': int(total_conformers),
                'avg_conformers_per_molecule': float(avg_conformers)
            },
            endpoint=self.backend_instance.get_endpoint()
        )

    def extract_molecules(self):
        """
        Extract generated molecules from backend.

        Yields:
            RDKit molecules with conformers attached

        Example:
            >>> output = actor(input_data)
            >>> for mol in actor.extract_molecules():
            ...     print(f"{mol.GetProp('_Name')}: {mol.GetNumConformers()} conformers")
        """
        return self.backend_instance.extract_molecules()
