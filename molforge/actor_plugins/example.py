"""
Plugin Actor Example and Template

This file demonstrates the complete pattern for creating custom plugin actors
in molforge. Plugin actors extend the pipeline with domain-specific functionality
while maintaining compatibility with the core framework.

Plugin Structure:
    1. Parameter class: Defines configurable parameters with validation
    2. Actor class: Implements processing logic following BaseActor protocol
    3. Optional dependencies: Handle external libraries gracefully

Usage:
    Place this file in molforge/actor_plugins/ and it will be automatically
    discovered. Reference the actor in pipeline configurations using the
    __step_name__ identifier.

Use this module as a starting point when writing your own plugin actors.
"""

from dataclasses import dataclass
from typing import List, Optional
import pandas as pd

from molforge.actors.base import BaseActor
from molforge.actors.protocol import ActorOutput
from molforge.actors.params.base import BaseParams


# ============================================================================
# PARAMETER CLASS
# ============================================================================

@dataclass
class CalculatePropertiesParams(BaseParams):
    """
    Configuration parameters for molecular property calculation.

    Defines all configurable options for the property calculation actor,
    including validation logic that runs before execution begins. See the
    per-field docstrings below for details on each parameter.
    """

    # Actor-specific parameters
    properties: Optional[List[str]] = None
    """List of RDKit descriptor names to calculate. Each name must exist in
    the rdkit.Chem.Descriptors module. Defaults to a set of common descriptors
    when left as None."""

    smiles_column: str = 'curated_smiles'
    """DataFrame column containing the SMILES strings to process."""

    filter_invalid: bool = False
    """Whether to drop molecules whose descriptor calculation failed (rows with
    None values in the computed property columns)."""

    def __post_init__(self):
        """
        Initialize default values and trigger validation.

        Called automatically after dataclass initialization. Sets default
        property list if none provided, then validates all parameters.
        """
        if self.properties is None:
            # Default to common molecular descriptors
            self.properties = ['ExactMolWt', 'MolLogP', 'TPSA', 'NumHDonors', 'NumHAcceptors']

        # Call parent initialization (triggers _validate_params)
        super().__post_init__()

    def _validate_params(self) -> None:
        """
        Validate parameter values before actor initialization.

        This method is called during parameter initialization and should raise
        exceptions for invalid configurations. Early validation prevents runtime
        errors during pipeline execution.

        Raises:
            ImportError: If RDKit is not available.
            AttributeError: If specified property names are invalid.
            ValueError: If properties list is empty.
        """
        # Check RDKit availability
        try:
            from rdkit.Chem import Descriptors
        except ImportError:
            raise ImportError(
                "RDKit is required for property calculation. "
                "Install with: pip install rdkit"
            )

        # Validate properties list is not empty
        if not self.properties or len(self.properties) == 0:
            raise ValueError("Properties list cannot be empty")

        # Validate each property name exists in RDKit
        for prop_name in self.properties:
            if not hasattr(Descriptors, prop_name):
                raise AttributeError(
                    f"Property '{prop_name}' not found in rdkit.Chem.Descriptors. "
                    f"Available descriptors can be listed with dir(Descriptors)."
                )


# ============================================================================
# ACTOR CLASS
# ============================================================================

class CalculateProperties(BaseActor):
    """
    Calculate molecular descriptors for molecules in the pipeline.

    This actor computes RDKit molecular descriptors and adds them as new columns
    to the pipeline DataFrame. Each specified property is calculated for every
    molecule, enabling downstream filtering, analysis, or modeling.

    The actor demonstrates standard plugin patterns:
        - Parameter-driven configuration
        - Input validation using required_columns
        - Output specification via output_columns
        - Robust error handling for individual molecules
        - Integration with pipeline logging system

    Example:
        Configure in pipeline YAML:
            steps: [source, curate, properties]
            properties:
                properties: [NumHDonors, NumHAcceptors]
                filter_invalid: true

        Or programmatically:
            params = CalculatePropertiesParams(
                properties=['ExactMolWt', 'HeavyAtomCount'],
                smiles_column='smiles'
            )
            actor = CalculateProperties(params)
            result = actor.process(dataframe)
    """

    # Plugin metadata (required)
    __step_name__ = 'properties'  # Used in pipeline configuration
    __param_class__ = CalculatePropertiesParams

    # Optional: declare upstream steps this plugin depends on. The registry
    # validates that each listed step appears earlier in the pipeline.
    # __dependencies__ = ['curate']  # (uncomment if this plugin requires curated SMILES)

    @property
    def required_columns(self) -> List[str]:
        """
        Specify columns that must exist in input DataFrame.

        The pipeline validates these requirements before actor execution,
        providing clear error messages if columns are missing.

        Returns:
            List of required column names.
        """
        return [self.smiles_column]

    @property
    def output_columns(self) -> List[str]:
        """
        Specify columns that will be added to DataFrame.

        This documents actor behavior and enables validation in downstream
        actors that depend on these columns. Here the output columns are
        dynamic: it returns the configured descriptor names, so the reported
        columns follow whatever properties the user requested.

        Returns:
            List of property names that will be added as columns.
        """
        return self.properties

    def __post_init__(self):
        """
        Initialize actor-specific state after parameter initialization.

        This optional method is called automatically after BaseActor.__init__
        completes. Use it for actor-specific setup that requires access to
        parameters or logging.

        Note:
            Avoid heavy computation here. This method is called during pipeline
            construction, not execution. Resource-intensive initialization should
            occur in process() method.
        """
        self.log(
            f"Property calculator initialized for {len(self.properties)} descriptors: "
            f"{', '.join(self.properties)}"
        )

        # Store descriptor functions for performance (avoid repeated getattr calls)
        from rdkit.Chem import Descriptors
        self._descriptor_functions = {
            prop: getattr(Descriptors, prop) for prop in self.properties
        }

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate molecular properties for all molecules in DataFrame.

        This is the core processing method that implements actor logic. It is
        called by the pipeline for each batch of molecules, with input validation
        and error handling managed by the BaseActor wrapper.

        Processing strategy:
            1. Create copy of input DataFrame to avoid mutations
            2. Parse each SMILES into an RDKit molecule exactly ONCE
            3. Calculate every descriptor for that molecule
            4. Handle calculation failures gracefully (None values)
            5. Optionally filter molecules with failed calculations

        Args:
            data: Input DataFrame containing SMILES column and any additional data
                  from upstream actors.

        Returns:
            DataFrame with original columns plus new property columns. If
            filter_invalid is True, rows with calculation failures are removed.

        Note:
            Calculation failures (invalid SMILES, descriptor errors) result in
            None values rather than exceptions. This prevents single molecule
            failures from terminating the entire pipeline.
        """
        from rdkit import Chem

        # Work on copy to avoid modifying pipeline state
        df = data.copy()

        self.log(
            f"Calculating {len(self.properties)} properties for {len(df)} molecules"
        )

        # Track molecule-level statistics for logging
        success_count = 0
        failure_count = 0

        # One list of values per requested descriptor, built in parallel as we
        # walk the molecules. This keeps the output as one column per property.
        columns = {prop_name: [] for prop_name in self.properties}

        # Parse each molecule ONCE, then compute all descriptors for it.
        for idx, smiles in enumerate(df[self.smiles_column]):
            mol = Chem.MolFromSmiles(smiles) if smiles is not None else None

            if mol is None:
                # Invalid/unparseable SMILES: no descriptors for this molecule.
                for prop_name in self.properties:
                    columns[prop_name].append(None)
                failure_count += 1
                continue

            # Molecule parsed successfully — compute each descriptor.
            for prop_name in self.properties:
                descriptor_func = self._descriptor_functions[prop_name]
                try:
                    columns[prop_name].append(descriptor_func(mol))
                except Exception as e:
                    # Individual descriptor failure — store None, keep going.
                    self.log(
                        f"Descriptor '{prop_name}' failed for molecule {idx}: {e}",
                        level='DEBUG'
                    )
                    columns[prop_name].append(None)
            success_count += 1

        # Attach the computed columns to the DataFrame.
        for prop_name in self.properties:
            df[prop_name] = columns[prop_name]

        # Log molecule-level summary statistics.
        self.log(
            f"Parsed {success_count}/{len(df)} molecules "
            f"({failure_count} failures)"
        )

        # Optional: Filter molecules with calculation failures
        if self.filter_invalid:
            initial_count = len(df)
            df = df.dropna(subset=self.properties)
            removed_count = initial_count - len(df)

            if removed_count > 0:
                self.log(f"Filtered {removed_count} molecules with calculation failures")

        return df

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """
        Create actor output with additional metadata.

        This method is called by BaseActor after successful processing to wrap
        the result DataFrame in an ActorOutput object. Override this method to
        add custom metadata or specify output endpoints.

        Args:
            data: Processed DataFrame from process() method.

        Returns:
            ActorOutput containing processed data and metadata.
        """
        # Calculate summary statistics for metadata
        property_stats = {}
        for prop in self.properties:
            if prop in data.columns:
                # Include basic statistics for each property
                property_stats[prop] = {
                    'mean': float(data[prop].mean()),
                    'std': float(data[prop].std()),
                    'min': float(data[prop].min()),
                    'max': float(data[prop].max()),
                    'n_valid': int(data[prop].notna().sum())
                }

        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'n_molecules': len(data),
                'properties_calculated': self.properties,
                'property_statistics': property_stats
            }
        )
