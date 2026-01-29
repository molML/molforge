"""Configuration parameters for distribution-based molecular curation."""

from typing import Optional, List, Dict
from dataclasses import dataclass, field
import warnings

from .base import BaseParams


@dataclass
class PropertyThreshold(BaseParams):
    """
    Threshold configuration for a molecular property.
    
    Supports three types of thresholds (priority: absolute > statistical > quantile):
    - Absolute: Fixed min/max values
    - Statistical: Mean ± n*std
    - Quantile: Percentile-based bounds
    """
    
    # Absolute thresholds (highest priority)
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    
    # Statistical thresholds (mean ± n*std)
    statistical_lower: Optional[float] = None  # e.g., -2.0 for mean - 2*std
    statistical_upper: Optional[float] = None  # e.g., 2.0 for mean + 2*std
    
    # Quantile thresholds (percentiles)
    quantile_lower: Optional[float] = None  # e.g., 0.025 for 2.5th percentile
    quantile_upper: Optional[float] = None  # e.g., 0.975 for 97.5th percentile
    
    def _validate_params(self) -> None:
        """Validate threshold configuration."""
        
        # Validate absolute bounds
        if self.min_value is not None and self.max_value is not None:
            if self.min_value >= self.max_value:
                raise ValueError(
                    f"min_value ({self.min_value}) must be less than "
                    f"max_value ({self.max_value})"
                )
        
        # Validate statistical bounds
        if self.statistical_lower is not None and self.statistical_upper is not None:
            if self.statistical_lower >= self.statistical_upper:
                raise ValueError(
                    f"statistical_lower ({self.statistical_lower}) must be less than "
                    f"statistical_upper ({self.statistical_upper})"
                )
        
        # Validate quantile bounds
        if self.quantile_lower is not None:
            if not 0 <= self.quantile_lower < 1:
                raise ValueError(
                    f"quantile_lower must be in [0, 1), got {self.quantile_lower}"
                )
        
        if self.quantile_upper is not None:
            if not 0 < self.quantile_upper <= 1:
                raise ValueError(
                    f"quantile_upper must be in (0, 1], got {self.quantile_upper}"
                )
        
        if self.quantile_lower is not None and self.quantile_upper is not None:
            if self.quantile_lower >= self.quantile_upper:
                raise ValueError(
                    f"quantile_lower ({self.quantile_lower}) must be less than "
                    f"quantile_upper ({self.quantile_upper})"
                )
        
        # Warn about potentially problematic statistical bounds
        if self.statistical_lower is not None and self.statistical_lower > 0:
            warnings.warn(
                f"statistical_lower is positive ({self.statistical_lower}), meaning the "
                f"lower bound is ABOVE the mean. This may create an overly restrictive range.",
                UserWarning
            )
        
        if self.statistical_upper is not None and self.statistical_upper < 0:
            warnings.warn(
                f"statistical_upper is negative ({self.statistical_upper}), meaning the "
                f"upper bound is BELOW the mean. This may create an overly restrictive range.",
                UserWarning
            )
        
        # Warn about conflicting threshold types
        threshold_types = sum([
            self.min_value is not None or self.max_value is not None,
            self.statistical_lower is not None or self.statistical_upper is not None,
            self.quantile_lower is not None or self.quantile_upper is not None
        ])
        
        if threshold_types > 1:
            warnings.warn(
                f"Multiple threshold types detected. Priority: absolute > statistical > quantile. "
                f"Only the highest priority threshold will be applied.",
                UserWarning
            )


@dataclass
class CurateDistributionParams(BaseParams):
    """Configuration parameters for distribution-based molecular curation."""
    
    # Column names
    SMILES_column: str = 'curated_smiles'
    tokens_column: str = 'tokens'
    
    # Property-specific thresholds
    thresholds: Dict[str, PropertyThreshold] = field(default_factory=dict)
    
    # Global threshold rules (applied to all properties without explicit thresholds)
    global_statistical_threshold: Optional[float] = None  # e.g., 2.0 for ±2σ
    global_quantile_threshold: Optional[float] = None     # e.g., 0.025 for 2.5%-97.5%
    
    # Property computation control
    properties: Optional[List[str]] | str = None  # 'all', list of properties, or None (auto)
    
    # Valid molecular properties
    _VALID_PROPERTIES: List[str] = field(default_factory=lambda: [
        "num_atoms", "num_rings", "size_largest_ring", "num_tokens",
        "tokens_atom_ratio", "c_atom_ratio", "longest_aliph_carbon",
        "molecular_weight", "logp", "tpsa", "num_rotatable_bonds",
        "num_h_donors", "num_h_acceptors", "fsp3", "num_aromatic_rings",
        "num_stereocenters", "num_heteroatoms", "heteroatom_ratio"
    ])
    
    # Token curation
    curate_tokens: bool = True
    filter_unknown_tokens: bool = True
    token_frequency_threshold: Optional[float] = None  # Percentage (0-100)
    
    # Processing options
    compute_properties: bool = True
    dropna: bool = True
    
    # Visualization
    plot_distributions: bool = True
    perform_pca: bool = False
    
    def _validate_params(self) -> None:
        """Validate parameter configuration."""
        
        # Validate token curation settings
        if not self.curate_tokens:
            if self.token_frequency_threshold is not None:
                raise ValueError(
                    "token_frequency_threshold requires curate_tokens=True"
                )
            if self.filter_unknown_tokens:
                raise ValueError(
                    "filter_unknown_tokens requires curate_tokens=True"
                )
        
        if self.token_frequency_threshold is not None:
            if not 0 < self.token_frequency_threshold < 100:
                raise ValueError(
                    f"token_frequency_threshold must be in (0, 100), "
                    f"got {self.token_frequency_threshold}"
                )
        
        # Validate properties parameter
        if self.properties is not None and self.properties != 'all':
            if not isinstance(self.properties, list):
                raise ValueError(
                    f"properties must be 'all', a list, or None, got {type(self.properties)}"
                )
            
            invalid_props = [p for p in self.properties if p not in self._VALID_PROPERTIES]
            if invalid_props:
                raise ValueError(
                    f"Invalid properties: {invalid_props}. "
                    f"Valid: {self._VALID_PROPERTIES}"
                )
        
        # Validate global thresholds (mutually exclusive)
        if self.global_quantile_threshold is not None and self.global_statistical_threshold is not None:
            raise ValueError(
                "Cannot set both global_quantile_threshold and global_statistical_threshold. "
                "Please choose one."
            )
        
        if self.global_statistical_threshold is not None:
            if self.global_statistical_threshold <= 0:
                raise ValueError(
                    f"global_statistical_threshold must be positive, "
                    f"got {self.global_statistical_threshold}"
                )
        
        if self.global_quantile_threshold is not None:
            if not 0 < self.global_quantile_threshold < 0.5:
                raise ValueError(
                    f"global_quantile_threshold must be in (0, 0.5), "
                    f"got {self.global_quantile_threshold}"
                )
    
    def _post_init_hook(self):
        # Apply global thresholds (if set) - quantile takes precedence over statistical
        if self.global_quantile_threshold is not None:
            if self.properties is None:
                self.properties = 'all'
                
            for name in self._VALID_PROPERTIES:
                if name not in self.thresholds or self.thresholds[name] is None:
                    self.thresholds[name] = PropertyThreshold(
                        quantile_lower=self.global_quantile_threshold,
                        quantile_upper=1.0 - self.global_quantile_threshold
                    )
        elif self.global_statistical_threshold is not None:
            if self.properties is None:
                self.properties = 'all'
                
            for name in self._VALID_PROPERTIES:
                if name not in self.thresholds or self.thresholds[name] is None:
                    self.thresholds[name] = PropertyThreshold(
                        statistical_lower=-self.global_statistical_threshold,
                        statistical_upper=self.global_statistical_threshold
                    )
        
        # Process user-defined thresholds
        for name in self._VALID_PROPERTIES:
            if name not in self.thresholds or self.thresholds[name] is None:
                self.thresholds[name] = PropertyThreshold()
                
            elif isinstance(self.thresholds[name], dict):
                try:
                    self.thresholds[name] = PropertyThreshold(**self.thresholds[name])
                except TypeError as e:
                    raise ValueError(f"Invalid threshold configuration for '{name}': {e}")
                    
            elif not isinstance(self.thresholds[name], PropertyThreshold):
                raise TypeError(f"property threshold '{name}' must be either None, dict or PropertyThreshold, got {type(self.thresholds[name])}")