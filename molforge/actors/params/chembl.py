from typing import List, Dict, Literal, Optional, NamedTuple
from dataclasses import dataclass, field

from .base import BaseParams


# ==================== Type Registry (Internal) ====================

class TypeSpec(NamedTuple):
    """
    Specification defining curation behavior for a ChEMBL standard_type.

    Attributes:
        category: Data category (potency, kinetic, admet, activity)
        requires_pchembl: Whether pchembl_value is required for quality control
        requires_units: Whether standard_units must be specified
        apply_log_transform: Whether to apply -log10 transformation during standardization
        apply_suspicious_pairs: Whether to check for unit conversion errors
        default_units: Default units if not specified by user
        description: Human-readable description of the measurement type
    """
    category: str
    requires_pchembl: bool
    requires_units: bool
    apply_log_transform: bool
    apply_suspicious_pairs: bool
    default_units: str | None
    description: str


# Module-level registry mapping standard_type strings to curation behavior
_TYPE_REGISTRY: Dict[str, TypeSpec] = {
    # Potency: Concentration-based affinity/inhibition measurements
    'IC50': TypeSpec('potency', True, True, True, True, 'nM', 'Half-maximal inhibitory concentration'),
    'EC50': TypeSpec('potency', True, True, True, True, 'nM', 'Half-maximal effective concentration'),
    'AC50': TypeSpec('potency', True, True, True, True, 'nM', 'Half-maximal activity concentration'),
    'Ki': TypeSpec('potency', True, True, True, True, 'nM', 'Inhibition constant (equilibrium dissociation)'),
    'Kd': TypeSpec('potency', True, True, True, True, 'nM', 'Dissociation constant'),
    'XC50': TypeSpec('potency', True, True, True, True, 'nM', 'Generic half-maximal concentration'),

    # Kinetic: Rate constants for binding kinetics
    'kon': TypeSpec('kinetic', False, True, False, False, 'M-1.s-1', 'Association rate constant'),
    'k_off': TypeSpec('kinetic', False, True, False, False, 's-1', 'Dissociation rate constant'),
    'koff': TypeSpec('kinetic', False, True, False, False, 's-1', 'Dissociation rate constant (alternate)'),
    'Kon': TypeSpec('kinetic', False, True, False, False, 'M-1.s-1', 'Association rate constant (capitalized)'),
    'Koff': TypeSpec('kinetic', False, True, False, False, 's-1', 'Dissociation rate constant (capitalized)'),

    # ADMET: Pharmacokinetic and physicochemical properties
    'Log BB': TypeSpec('admet', False, False, False, False, None, 'Blood-brain barrier penetration (log ratio)'),
    'Log D': TypeSpec('admet', False, False, False, False, None, 'Distribution coefficient at pH 7.4'),
    'Log P': TypeSpec('admet', False, False, False, False, None, 'Partition coefficient (octanol/water)'),
    'Clearance': TypeSpec('admet', False, True, False, False, 'mL.min-1.kg-1', 'Metabolic clearance rate'),
    'T1/2': TypeSpec('admet', False, True, False, False, 'h', 'Plasma half-life'),
    'Solubility': TypeSpec('admet', False, True, False, False, 'uM', 'Aqueous solubility'),
    'Permeability': TypeSpec('admet', False, True, False, False, 'cm.s-1', 'Membrane permeability'),
    'Caco-2': TypeSpec('admet', False, True, False, False, 'cm.s-1', 'Caco-2 cell permeability'),
    'MDCK': TypeSpec('admet', False, True, False, False, 'cm.s-1', 'MDCK cell permeability'),

    # Activity: Percent-based and relative activity measures
    'Inhibition': TypeSpec('activity', False, True, False, False, '%', 'Percent inhibition'),
    'Percent Effect': TypeSpec('activity', False, True, False, False, '%', 'Percent effect'),
    'Activity': TypeSpec('activity', False, False, False, False, None, 'Generic activity measure'),
    'Potency': TypeSpec('activity', False, False, False, False, None, 'Generic potency measure'),
}


def _get_type_spec(standard_type: str) -> TypeSpec:
    """Get type specification for a standard_type."""
    if standard_type not in _TYPE_REGISTRY:
        available = sorted(_TYPE_REGISTRY.keys())
        raise ValueError(
            f"Unsupported standard_type '{standard_type}'. "
            f"Supported types ({len(available)}): {', '.join(available)}"
        )
    return _TYPE_REGISTRY[standard_type]


def _list_types_by_category(category: str) -> list[str]:
    """Get standard_types for a specific category."""
    return sorted([
        stype for stype, spec in _TYPE_REGISTRY.items()
        if spec.category == category
    ])

@dataclass
class ChEMBLCuratorParams(BaseParams):
    """
    Configuration for ChEMBL assay curation with automatic type-driven behavior.

    Curation logic adapts based on standard_type category:
    - Potency (IC50, Ki): pchembl required, log transform, suspicious pair detection
    - Kinetic (kon, k_off): no pchembl, no transforms, different units
    - ADMET (Log BB, Clearance): no pchembl, pass-through values
    - Activity (% Inhibition): no pchembl, no transforms

    Use ChEMBLCuratorParams.list_supported_types() to see all supported types.
    """

    # ==================== Core Parameters ====================

    standard_type: str = 'IC50'
    """ChEMBL standard_type to curate (drives type-specific curation behavior)."""
    standard_units: Optional[str] = None
    """Target standard_units; auto-set from the type registry when None."""
    standard_relation: Optional[str] = '='
    """Activity relation to keep (e.g. '=', '>', '<')."""
    assay_type: Literal['B', 'F', 'A', 'T'] = 'B'
    """ChEMBL assay type: 'B' (binding), 'F' (functional), 'A' (ADME), or 'T' (toxicity)."""
    target_organism: str = 'Homo sapiens'
    """Target organism to filter assays by."""
    assay_format: Literal['general', 'protein', 'cell', 'organism', 'tissue', 'microsome'] = 'protein'
    """Assay format; mapped to a BAO format code in _post_init_hook."""

    top_n: int = 5
    """Number of top experimental-condition groupings shown in the logged RAW/VALID analysis tables (reporting only; does not affect curation)."""

    mutant_regex: str = 'mutant|mutation|variant'
    """Regex used to detect and exclude mutant/variant assays."""
    allosteric_regex: str = 'allosteric'
    """Regex used to detect and exclude allosteric assays."""

    C: float = 1e9
    """Molarity conversion constant, auto-set from standard_units in _post_init_hook. Only change for non-standard unit conversion."""

    mistakes_only: bool = True
    """If True, flag a SMILES as suspicious only when its repeated standardized values differ by ~3 log units (or ~1 for Ki/Kd); if False, also flag SMILES whose repeated values are (near-)identical. Flagged SMILES are dropped by filter_suspicious_pairs."""
    error_margin: float = 0.0001
    """Numerical tolerance used when comparing values for mistake detection."""

    SMILES_column: str = 'canonical_smiles'
    """DataFrame column containing SMILES strings."""
    std_threshold: float = 0.5
    """Maximum per-SMILES standard deviation of standardized values for a duplicate group to be kept during aggregation (applied to all standard_types)."""
    range_threshold: float = 0.5
    """Maximum per-SMILES value range (max - min of standardized values) for a duplicate group to be kept during aggregation (applied to all standard_types)."""

    bao_format: Optional[str] = None
    """Derived BioAssay Ontology format code; overwritten in _post_init_hook from assay_format (do not set manually)."""

    # ==================== Private Fields (excluded from serialization) ====================

    _type_spec: Optional[TypeSpec] = field(default=None, init=False, repr=False)

    _ALLOWED_VALUES: Dict[str, List[str | None]] = field(default_factory=lambda: {
        'standard_units':    ['nM', 'uM', 'mM', 'M', 'pM', None],
        'standard_relation': ['=', '>', '<', '>=', '<=', '~', None],
        'assay_type':        ['B', 'F', 'A', 'T'],
        'assay_format':      ['general', 'protein', 'cell', 'organism', 'tissue', 'microsome'],
    }, init=False, repr=False)
    
    _BAO_FORMAT_MAP: Dict[str, str] = field(default_factory=lambda: {
        # BioAssay Ontology (BAO) format codes for assay system classification
        # See: http://www.bioassayontology.org/

        # Potency-focused formats
        'general': 'BAO_0000019',       # assay format (unspecified)
        'protein': 'BAO_0000357',       # single protein format (biochemical, 51% of binding assays)
        'cell': 'BAO_0000219',          # cell-based format (67% of functional assays)

        # ADMET-focused formats
        'organism': 'BAO_0000218',      # organism-based format (in vivo, 93% of ADMET assays)
        'tissue': 'BAO_0000221',        # tissue-based format (ex vivo preparations)
        'microsome': 'BAO_0000251',     # microsome format (metabolic stability studies)
        'subcellular': 'BAO_0000220',   # subcellular format (organelles, fractions)

        # Specialized formats
        'membrane': 'BAO_0000249',      # cell membrane format (membrane preparations)
        'nucleic_acid': 'BAO_0000225',  # nucleic acid format (DNA/RNA targeting)
    }, init=False, repr=False)


    _UNIT_CONVERSION: Dict[str, float] = field(default_factory=lambda: {
        'pM': 0.001,
        'nM': 1.0,
        'uM': 1000.0,
        'mM': 1e6,
        'M': 1e9,
    }, init=False, repr=False)
    
    def _validate_params(self) -> None:
        """Validate parameters against allowed values and internal consistency."""
        # Validate allowed values (skip standard_units - type-specific)
        for param, allowed in self._ALLOWED_VALUES.items():
            if param == 'standard_units':
                continue  # Units are type-specific, validated by registry
            if hasattr(self, param):
                value = getattr(self, param)
                self._validate_policy(param, value, allowed)

        # Type-specific validation (requires _type_spec from _post_init_hook)
        if self._type_spec:
            if self._type_spec.requires_units and self.standard_units is None:
                raise ValueError(
                    f"Type '{self.standard_type}' ({self._type_spec.category}) requires standard_units. "
                    f"Default: {self._type_spec.default_units}"
                )

    def _post_init_hook(self) -> None:
        """Post-initialization hook for type lookup and auto-configuration."""
        # Lookup type specification
        try:
            self._type_spec = _get_type_spec(self.standard_type)
        except ValueError as e:
            raise ValueError(
                f"{e}\nTo see supported types, use: "
                f"ChEMBLCuratorParams.list_supported_types()"
            )

        # Auto-set units from registry if not provided
        if self.standard_units is None:
            self.standard_units = self._type_spec.default_units

        # Map BAO format
        self.bao_format = self._BAO_FORMAT_MAP[self.assay_format]

        # Set conversion constant (only for potency types with unit conversion)
        if self.standard_units in self._UNIT_CONVERSION:
            self.C = self._UNIT_CONVERSION[self.standard_units] * 1e9

    # ==================== Type Query Methods ====================

    def get_category(self) -> str:
        """Get data category (potency, kinetic, admet, activity)."""
        return self._type_spec.category if self._type_spec else 'unknown'

    def requires_pchembl(self) -> bool:
        """Whether this type requires pchembl_value for quality control."""
        return self._type_spec.requires_pchembl if self._type_spec else False

    def should_apply_log_transform(self) -> bool:
        """Whether to apply -log10 transformation during standardization."""
        return self._type_spec.apply_log_transform if self._type_spec else False

    def should_check_suspicious_pairs(self) -> bool:
        """Whether to check for unit conversion errors (3-log, 1-log differences)."""
        return self._type_spec.apply_suspicious_pairs if self._type_spec else False

    @classmethod
    def list_supported_types(cls) -> list[str]:
        """Get all supported standard_types."""
        return sorted(_TYPE_REGISTRY.keys())

    @classmethod
    def get_potency_types(cls) -> list[str]:
        """Get all potency-category types."""
        return _list_types_by_category('potency')

    @classmethod
    def get_kinetic_types(cls) -> list[str]:
        """Get all kinetic-category types."""
        return _list_types_by_category('kinetic')

    @classmethod
    def get_admet_types(cls) -> list[str]:
        """Get all ADMET-category types."""
        return _list_types_by_category('admet')


# ==================== Example Configurations ====================

def get_ic50_config():
    """
    Standard IC50 curation configuration for potency assays.

    Example for binding assays with human targets.
    """
    return ChEMBLCuratorParams(
        standard_type='IC50',
        standard_units='nM',
        standard_relation='=',
        assay_type='B',
        target_organism='Homo sapiens',
        assay_format='protein'
    )

def get_ki_config():
    """
    Standard Ki curation configuration for equilibrium binding constants.

    Uses tighter thresholds than IC50 due to higher measurement precision.
    """
    return ChEMBLCuratorParams(
        standard_type='Ki',
        standard_units='nM',
        standard_relation='=',
        assay_type='B',  # Binding assays for Ki
        target_organism='Homo sapiens',
        assay_format='protein',
        std_threshold=0.3,  # Tighter threshold for Ki
        range_threshold=0.3  # Tighter threshold for Ki
    )

def get_logbb_config():
    """
    Log BB (blood-brain barrier penetration) configuration for ADMET assays.

    Typical for mouse BBB studies. Log BB values are unitless log ratios
    (log(C_brain/C_blood)), typically ranging from -3 to +2.
    """
    return ChEMBLCuratorParams(
        standard_type='Log BB',
        standard_units=None,  # Unitless (log ratio)
        standard_relation=None,  # Not applicable for measurements
        assay_type='A',  # ADME assay
        target_organism='Mus musculus',  # Often measured in mice
        assay_format='organism',  # In vivo organism-based
        std_threshold=0.3,  # Tighter for ADMET
        range_threshold=0.3
    )

def get_logd_config():
    """
    Log D (distribution coefficient at pH 7.4) configuration.

    Log D values are unitless partition coefficients, typically -5 to +5.
    """
    return ChEMBLCuratorParams(
        standard_type='Log D',
        standard_units=None,  # Unitless
        standard_relation=None,
        assay_type='A',  # ADME assay
        target_organism='Homo sapiens',
        assay_format='organism',  # In vivo organism-based
        std_threshold=0.3,
        range_threshold=0.3
    )