"""
ChEMBL data curation module for bioactivity data processing.

This module provides comprehensive curation functionality for ChEMBL bioactivity data,
including standardization of activity values, handling of duplicates, detection of 
suspicious data entries, and validation of experimental conditions.
"""

from typing import Optional, List, Set, Dict, Any, Tuple
import pandas as pd
import numpy as np

from .params.chembl import ChEMBLCuratorParams
from .base import BaseActor
from .protocol import ActorOutput
from ..configuration.steps import Steps


class ChEMBLCurator(BaseActor):
    """
    ChEMBL bioactivity and ADMET data curator with comprehensive quality control.

    Supports two major data paradigms:

    **1. Potency Assays** (IC50, Ki, EC50, Kd):
       - Concentration-based measurements with units (nM, µM, etc.)
       - Converts to -log10 scale (pIC50, pKi, etc.)
       - Requires standard_units and standard_relation

    **2. ADMET Properties** (Log BB, Log D, Permeability):
       - Often unitless or already in log scale
       - Values used as-is (no conversion)
       - standard_units and standard_relation optional/not applicable

    Pipeline operations:
    - Validates experimental conditions and assay consistency
    - Standardizes activity measurements (conditional on data type)
    - Detects and filters suspicious data points (transcription errors)
    - Handles duplicate SMILES entries with statistical aggregation

    Examples:
        >>> # IC50 potency assay
        >>> params = ChEMBLCuratorParams(
        ...     standard_type='IC50',
        ...     standard_units='nM',
        ...     assay_type='B'
        ... )

        >>> # Log BB for blood-brain barrier penetration
        >>> params = ChEMBLCuratorParams(
        ...     standard_type='Log BB',
        ...     standard_units=None,  # Unitless
        ...     standard_relation=None,  # Not applicable
        ...     assay_type='A',  # ADME assay
        ...     target_organism='Mus musculus'
        ... )

    Note: This class does NOT perform molecular standardization (desalting,
    stereochemistry, neutralization, etc.). Use CurateMol for SMILES curation.
    """
    __step_name__ = Steps.CHEMBL

    def __post_init__(self):
        """Post-initialization setup for output label."""
        # Placeholder for output label - updated during standardization
        self.output_label = "[WARNING] Activity values have not been standardized yet."

    @property
    def required_columns(self) -> List[str]:
        """Required input columns from ChEMBL API/SQL."""
        return [
            self.SMILES_column,
            'standard_type', 'standard_units', 'standard_relation',
            'standard_flag', 'potential_duplicate', 'standard_value',
            'pchembl_value', 'data_validity_comment',
            'document_year', 'document_chembl_id',
            'assay_type', 'assay_description', 'target_organism', 'bao_format',
            'molecule_chembl_id'
        ]

    @property
    def output_columns(self) -> List[str]:
        """Output columns depend on standard_type (dynamic)."""
        # Return placeholder, actual output_label set during processing
        return [self.output_label]

    def __repr__(self) -> str:
        """Return string representation with configuration parameters."""
        return f"ChEMBLCurator({self._params})"

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Curate ChEMBL bioactivity data through comprehensive pipeline.

        Performs the following operations:
        1. Validates experimental conditions
        2. Filters by assay type and quality criteria
        3. Standardizes activity values to common units
        4. Detects and removes suspicious data pairs
        5. Aggregates duplicate SMILES entries

        Args:
            data: Raw ChEMBL activity DataFrame

        Returns:
            Curated DataFrame with standardized activities

        Raises:
            ValueError: If SMILES column not found or invalid standard_type
        """
        if self.SMILES_column not in data.columns:
            raise ValueError(
                f"'{self.SMILES_column}' not in dataframe. "
                f"Available columns: {list(data.columns)}"
            )

        # Analyze conditions if verbose mode enabled
        if self.verbose:
            self.conditions_analysis(data)

        # Apply curation pipeline
        curated = data.copy()
        curated = self.curate_activies(curated)
        curated = self.standardize_acitivities(curated)
        curated = self.filter_suspicious_pairs(curated)
        curated = self.handle_duplicate_smiles(curated)

        return curated

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with ChEMBL curation metadata."""
        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'output_label': self.output_label,
                'standard_type': self.standard_type,
                'n_curated': len(data)
            }
        )
    
    # ==================== Condition Analysis ====================
    
    def conditions_analysis(self, activities: pd.DataFrame) -> None:
        """
        Analyze experimental conditions to validate configuration optimality.
        
        Compares the configured curation parameters against the most common
        experimental conditions in the dataset to ensure optimal data retention.
        
        Args:
            activities: Input activity DataFrame
        """
        fields = ['standard_type', 'assay_type', 'target_organism', 'bao_format']
        
        # Apply basic validity filters
        valid_activities = self._get_valid_activities(activities)
        
        # Analyze both raw and valid data
        self._show_condition_analysis(activities, 'RAW', fields)
        self._show_condition_analysis(valid_activities, 'VALID', fields)
    
    # ==================== Primary Curation ====================
    
    def curate_activies(self, activities: pd.DataFrame) -> pd.DataFrame:
        """
        Apply primary curation filters based on experimental conditions.

        Filters activities based on:
        - Standard type (required)
        - Units and relation (conditional - only for potency assays)
        - Data quality flags
        - Assay metadata
        - Exclusion of mutants and allosteric modulators

        Args:
            activities: Input activity DataFrame

        Returns:
            Filtered DataFrame
        """
        # Build conditional mask for pChEMBL requirement
        pchembl_condition = self._get_pchembl_condition(activities)

        # Build filtering conditions
        conditions = [
            # Core type filtering
            (activities['standard_type'] == self.standard_type),
            (activities['standard_flag'] == 1),
            (activities['potential_duplicate'] == 0),  # Important for data leakage

            # Data validity
            pchembl_condition,
            (activities['standard_value'].notna()),
            (activities[self.SMILES_column].notna()),

            # Metadata validity
            (activities['document_year'].notna()),
            (activities['document_chembl_id'].notna()),
            (activities['data_validity_comment'].isna()),

            # Assay curation
            (activities['assay_type'] == self.assay_type),
            (activities['target_organism'] == self.target_organism),
            (activities['bao_format'] == self.bao_format),

            # Exclude mutants and allosteric modulators
            (~activities['assay_description'].str.contains(
                self.mutant_regex, case=False, na=False)),
            (~activities['assay_description'].str.contains(
                self.allosteric_regex, case=False, na=False))
        ]

        # Conditionally add units filter (only if specified)
        if self.standard_units is not None:
            conditions.append(activities['standard_units'] == self.standard_units)

        # Conditionally add relation filter (only if specified)
        if self.standard_relation is not None:
            conditions.append(activities['standard_relation'] == self.standard_relation)

        # Combine all conditions
        combined_mask = conditions[0]
        for condition in conditions[1:]:
            combined_mask = combined_mask & condition

        cleaned = activities[combined_mask].copy()

        self.log(f'Conditional curation: {len(cleaned)}/{len(activities)}.')
        return cleaned
    
    # ==================== Value Standardization ====================
    
    def standardize_acitivities(self, activities: pd.DataFrame) -> pd.DataFrame:
        """
        Convert activity values to standardized units (type-aware).

        Behavior determined by standard_type category:
        - Potency: Convert to -log10 scale (pIC50, pKi, etc.)
        - Kinetic/ADMET/Activity: Use values as-is (no transformation)

        Args:
            activities: DataFrame with raw activity values

        Returns:
            DataFrame with standardized activity values
        """
        activities = activities.copy()
        standard_values = activities['standard_value'].astype(float)
        category = self._params.get_category()

        # Type-aware standardization
        if self._params.should_apply_log_transform():
            # Potency types: convert to -log10 scale
            conversion_map = {
                "IC50": (f"pIC50", self.molarity_to_pXC50),
                "XC50": (f"pXC50", self.molarity_to_pXC50),
                "EC50": (f"pEC50", self.molarity_to_pXC50),
                "AC50": (f"pAC50", self.molarity_to_pXC50),
                "ED50": (f"pED50", self.molarity_to_pXC50),
                "Ki": ("pKi", self.molarity_to_pKi),
                "Kd": ("pKd", self.molarity_to_pKd),
            }

            if self.standard_type in conversion_map:
                label, converter = conversion_map[self.standard_type]
                self.output_label = label
                activities[self.output_label] = standard_values.apply(converter)
                self.log(f'Standardization: {self.standard_type} → {self.output_label} (log10 conversion).')
            else:
                # Fallback for unmapped potency types
                self.output_label = f"p{self.standard_type}"
                activities[self.output_label] = -1 * np.log10(standard_values / self.C)
                self.log(f'Standardization: {self.standard_type} → {self.output_label} (generic log10 conversion).')

        else:
            # Non-potency types: pass through (kinetic, ADMET, activity)
            self.output_label = self.standard_type
            activities[self.output_label] = standard_values

            value_range = (standard_values.min(), standard_values.max())
            self.log(
                f'Standardization: {self.standard_type} ({category}) - no conversion. '
                f'Range: [{value_range[0]:.2g}, {value_range[1]:.2g}]'
            )

        return activities
    
    def molarity_to_pXC50(self, x: float) -> float:
        """
        Convert xC50 concentration to pXC50.
        
        Args:
            x: Concentration in nM (ChEMBL default)
            
        Returns:
            pXC50 value (-log10 of molar concentration)
        """
        return -1 * np.log10(x / self.C)
    
    def molarity_to_pKi(self, x: float) -> float:
        """
        Convert Ki to pKi.
        
        Ki is the equilibrium dissociation constant for the inhibitor.
        
        Args:
            x: Ki in nM (ChEMBL default)
            
        Returns:
            pKi value (-log10 of molar Ki)
        """
        return -1 * np.log10(x / self.C)
    
    def molarity_to_pKd(self, x: float) -> float:
        """
        Convert Kd to pKd.
        
        Kd is the equilibrium dissociation constant.
        
        Args:
            x: Kd in nM (ChEMBL default)
            
        Returns:
            pKd value (-log10 of molar Kd)
        """
        return -1 * np.log10(x / self.C)
    
    # ==================== Suspicious Data Detection ====================
    
    def detect_suspicious_pairs(self, activities: pd.DataFrame) -> Set[str]:
        """
        Detect SMILES with suspicious activity patterns.
        
        Identifies potential data errors based on:
        - Identical values (potential duplicates)
        - 3.0 log unit differences (1000-fold, common transcription error)
        - 1.0 log unit differences for Ki/Kd (10-fold, common reporting error)
        
        Reference: Kramer et al. on common errors in ChEMBL data
        
        Args:
            activities: DataFrame with standardized activities
            
        Returns:
            Set of suspicious SMILES strings
        """
        grouped = activities.groupby(self.SMILES_column)
        suspicious_smiles = set()
        
        for smiles, group in grouped:
            if len(group) < 2:
                continue
            
            values = group[self.output_label].values.astype(float)
            
            if self._has_suspicious_pattern(values):
                suspicious_smiles.add(smiles)
        
        log_level = 'WARNING' if suspicious_smiles else 'INFO'
        self.log(
            f'{len(suspicious_smiles)} suspicious SMILES detected. '
            f'mistakes_only={self.mistakes_only}, error_margin={self.error_margin}',
            level=log_level
        )

        return suspicious_smiles
    
    def filter_suspicious_pairs(self, activities: pd.DataFrame) -> pd.DataFrame:
        """
        Remove entries with suspicious activity patterns (type-aware).

        Only applies to potency types where unit conversion errors are common.

        Args:
            activities: Input DataFrame

        Returns:
            Filtered DataFrame without suspicious entries (or unchanged)
        """
        # Only check for suspicious pairs in potency assays
        if not self._params.should_check_suspicious_pairs():
            self.log(f'Skipping suspicious pair detection for {self.standard_type} ({self._params.get_category()}).')
            return activities

        suspicious = self.detect_suspicious_pairs(activities)
        filtered = activities[~activities[self.SMILES_column].isin(suspicious)]
        self.log(f'Dropped {len(activities) - len(filtered)} entries from {len(suspicious)} suspicious SMILES.')
        return filtered
    
    # ==================== Duplicate Handling ====================

    def handle_duplicate_smiles(self, activities: pd.DataFrame,
                                SMILES_column: Optional[str] = None) -> pd.DataFrame:
        """
        Process and aggregate duplicate SMILES entries with comprehensive logging.

        Performs statistical aggregation:
        1. Groups by SMILES, calculates stats (mean, std, range)
        2. Warns about conflicting ChEMBL IDs with MARKDOWN TABLE
        3. Filters groups by std/range thresholds
        4. Selects earliest occurrence as representative
        5. Logs detailed aggregation summary

        Aggregates multiple measurements for the same SMILES if they are
        consistent (low standard deviation and range).

        Args:
            activities: Input DataFrame
            SMILES_column: Column name for SMILES (uses default if None)

        Returns:
            DataFrame with aggregated duplicates and statistics
        """
        if SMILES_column is None:
            SMILES_column = self.SMILES_column

        # Calculate group statistics
        grouped = self._calculate_group_statistics(activities, SMILES_column)

        # Warn about conflicting ChEMBL IDs with full markdown table
        self._check_conflicting_ids(grouped)

        # Filter groups by quality thresholds
        grouped['keep'] = (
            (grouped['std'] < self.std_threshold) &
            (abs(grouped['max'] - grouped['min']) < self.range_threshold)
        )

        # Select representative rows from kept groups
        processed_rows = self._process_kept_groups(activities, grouped, SMILES_column)

        # Log aggregation summary
        self._log_aggregation_summary(grouped, len(activities))

        return pd.DataFrame(processed_rows).reset_index(drop=True)

    # ==================== Helper Methods ====================
    
    def _apply_params(self, params: ChEMBLCuratorParams) -> None:
        """Apply configuration parameters to instance."""
        for key, value in params.to_dict().items():
            setattr(self, key, value)
    
    def _get_valid_activities(self, activities: pd.DataFrame) -> pd.DataFrame:
        """Get activities passing basic validity checks."""
        pchembl_condition = self._get_pchembl_condition(activities)
        
        return activities[
            pchembl_condition &
            (activities['standard_flag'] == 1) &
            (activities['potential_duplicate'] == 0) &
            (activities['standard_value'].notna()) &
            (activities[self.SMILES_column].notna()) &
            (activities['document_year'].notna()) &
            (activities['document_chembl_id'].notna()) &
            (activities['data_validity_comment'].isna())
        ]
    
    def _get_pchembl_condition(self, activities: pd.DataFrame) -> pd.Series:
        """Get condition for pChEMBL value requirement (type-aware)."""
        if self._params.requires_pchembl():
            return activities['pchembl_value'].notna()
        return True  # No pchembl requirement for kinetic, ADMET, activity types
    
    def _show_condition_analysis(self, df: pd.DataFrame, label: str, 
                                 fields: List[str]) -> None:
        """Display analysis of experimental conditions."""
        # Group by fields and rank by frequency
        grouped = df.groupby(fields, dropna=False).size()
        top = pd.DataFrame(
            grouped.sort_values(ascending=False).head(self.top_n),
            columns=['entries']
        ).reset_index()

        # Calculate percentages
        top['%'] = (round(top['entries'] / len(df), 2) * 100).astype(int)
        top['entries'] = top['entries'].astype(str) + f"/{len(df)}"
        
        # Check if configuration is optimal
        if not top.empty:
            self._check_configuration_optimality(top, label, df)

        # Translate BAO codes to readable labels (after optimality check to not affect config slicing)
        bao_to_label = {v: k for k, v in self._params._BAO_FORMAT_MAP.items()}
        top.insert(top.columns.get_loc('bao_format') + 1, 'assay_format',
                   top['bao_format'].map(bao_to_label).fillna('unmapped'))

        self.log(f"{label}\n{top.to_markdown()}")
    
    def _check_configuration_optimality(self, top: pd.DataFrame, label: str, df: pd.DataFrame) -> None:
        """Check if current configuration matches most common conditions."""
        top_config = top.iloc[0][:-2].to_dict()
        current_config = {
            k: getattr(self, k) 
            for k in top_config.keys()
        }
        
        if top_config != current_config:
            self.log(
                f'{label} | Configured curation is not optimal.',
                level='WARNING'
            )
            bao_to_label = {v: k for k, v in self._params._BAO_FORMAT_MAP.items()}

            lines = []
            for (c_k, c_v), (_, t_v) in zip(current_config.items(), top_config.items()):
                if c_v == t_v:
                    continue
                if c_k == 'bao_format':
                    c_label = bao_to_label.get(c_v, 'unmapped')
                    t_label = bao_to_label.get(t_v, 'unmapped')
                    lines.append(('assay_format', f"{c_label} ({c_v})", f"{t_label} ({t_v})"))
                else:
                    lines.append((c_k, repr(c_v), repr(t_v)))

            max_key = max(len(k) for k, _, _ in lines)
            max_cv = max(len(cv) for _, cv, _ in lines)

            recommended_str = '\n'.join(
                f"  {k:{max_key}}:\t{cv:>{max_cv}} -> {tv}" for k, cv, tv in lines
            )
            
            self.log(
                f"Recommended configuration for better data retention:\n"
                f"{recommended_str}",
                level='WARNING'
            )

        else:
            target_id = df.target_chembl_id.unique()[0] if len(df.target_chembl_id.unique()) == 1 else "multiple targets"
            self.log(
                f'{label} | Configured curation is optimal for {target_id}.',
                level='INFO'
            )
    
    def _has_suspicious_pattern(self, values: np.ndarray) -> bool:
        """Check if values contain suspicious patterns."""
        for i in range(len(values)):
            for j in range(i + 1, len(values)):
                diff = abs(values[i] - values[j])
                
                # Check for identical or 3-log unit difference
                if ((not self.mistakes_only and diff <= self.error_margin) or 
                    abs(diff - 3.0) < self.error_margin):
                    return True
                
                # Additional check for Ki/Kd: 10-fold errors
                if (self.standard_type in ['Ki', 'Kd'] and 
                    abs(diff - 1.0) < self.error_margin):
                    return True
        
        return False
    
    def _calculate_group_statistics(self, activities: pd.DataFrame, 
                                   smiles_column: str) -> pd.DataFrame:
        """Calculate statistics for grouped SMILES."""
        grouped = activities.groupby(smiles_column).agg({
            self.output_label: ['mean', 'min', 'max', 'std', 'count'],
            'molecule_chembl_id': list,
            'document_year': list,
        }).fillna(0)
        
        # Flatten column names
        grouped.columns = [
            self.output_label, 'min', 'max', 'std', 'count',
            'molecule_chembl_ids', 'document_years'
        ]
        
        # Add earliest year information
        grouped['earliest_year'] = grouped['document_years'].apply(min).astype(int)
        grouped['earliest_idx'] = [
            years.index(earliest) 
            for years, earliest in zip(
                grouped['document_years'], 
                grouped['earliest_year']
            )
        ]
        
        return grouped
    
    def _check_conflicting_ids(self, grouped: pd.DataFrame) -> None:
        """Check for multiple ChEMBL IDs for same SMILES."""
        conflicting_ids = grouped['molecule_chembl_ids'].apply(
            lambda x: len(set(x)) > 1
        )
        
        if conflicting_ids.any():
            self.log(
                f"Multiple ChEMBL IDs exist for (canonical) SMILES: \n"
                f"{grouped[conflicting_ids].to_markdown()}",
                level='WARNING'
            )
            self.log(
                "Taking first entry(s) for the names of curated data.",
                level='WARNING'
            )
    
    def _process_kept_groups(self, activities: pd.DataFrame, grouped: pd.DataFrame,
                            smiles_column: str) -> List[pd.Series]:
        """Process groups that pass quality thresholds."""
        processed_rows = []
        
        for smiles, stats in grouped[grouped['keep']].iterrows():
            # Get earliest occurrence of this SMILES
            matching_rows = activities[activities[smiles_column] == smiles]
            row = matching_rows.iloc[stats['earliest_idx']].copy()
            
            # Add aggregated statistics
            for stat in [self.output_label, 'min', 'max', 'std', 'count']:
                row[stat] = stats[stat]
            
            processed_rows.append(row)
        
        return processed_rows
    
    def _log_aggregation_summary(self, grouped: pd.DataFrame, original_count: int) -> None:
        """Log summary statistics about duplicate aggregation."""
        smiles_to_remove = len(grouped[~grouped['keep']])
        smiles_to_keep = len(grouped[grouped['keep']])
        smiles_aggregated = sum(grouped[grouped['keep']]['count'] > 1)

        self.log(
            f"Removed {smiles_to_remove} (canonical) SMILES entries with "
            f"std > {self.std_threshold} or range > {self.range_threshold}."
        )
        self.log(
            f"Aggregated {smiles_aggregated} duplicate SMILES entries. "
            f"Final: {smiles_to_keep}/{original_count} rows."
        )