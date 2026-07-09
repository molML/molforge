"""
Clean parameter system supporting core actors (typed) and plugins (dynamic).
"""

from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field, fields
import warnings

from ..actors.params import *
from .registry import ActorRegistry
from .serialization import ConfigCodec


@dataclass
class PipeParams(BaseParams):
    """
    Configuration parameters for MolForge Pipeline.
    
    Supports both core actors (with explicit type hints) and plugin actors
    (with dynamic parameter handling).
    """
    
    # Pipeline configuration
    steps: List[str] = None
    """Ordered list of actor step names to execute (defaults to the standard 5-step pipeline when None)."""

    # Pipeline behavior
    return_none_on_fail: bool = False
    """On step failure, return None in place of the step output; with False, the
    partial step output is returned. Step errors are caught internally and
    surfaced through the return value."""
    output_root: str = "MOLFORGE_OUT"
    """Root directory for pipeline outputs."""
    write_output: bool = True
    """Write the final pipeline output DataFrame to disk (per-step checkpoints are
    controlled separately by write_checkpoints)."""
    write_checkpoints: bool = False
    """Write intermediate checkpoints between steps."""

    # Logging
    console_log_level: str = 'INFO'
    """Logging level for console output."""
    file_log_level: str = 'DEBUG'
    """Logging level for the log file."""

     # Parameters
    override_actor_params: bool = True
    """When True, copy each pipeline-level field value onto every actor's params
    (and any nested param dataclass) for any field whose name matches — e.g. the
    shared `verbose` flag. Matching top-level fields are always overwritten,
    regardless of the actor's own value."""

    # Core actors (explicit, type-hinted)
    source_params: Optional[ChEMBLSourceParams] = None
    """Parameters for the source actor (auto-created if None)."""
    chembl_params: Optional[ChEMBLCuratorParams] =  None
    """Parameters for the chembl actor (auto-created if None)."""
    curate_params: Optional[CurateMolParams]     =  None
    """Parameters for the curate actor (auto-created if None)."""
    tokens_params: Optional[TokenizeDataParams]  =  None
    """Parameters for the tokens actor (auto-created if None)."""
    distributions_params: Optional[CurateDistributionParams]  =  None
    """Parameters for the distributions actor (auto-created if None)."""
    confs_params:  Optional[GenerateConfsParams] =  None
    """Parameters for the confs actor (auto-created if None)."""

    # Plugin actors (dynamic, no type hints)
    plugin_params: Dict[str, Any] = field(default_factory=dict)
    """Per-step parameters for plugin actors, keyed by step name."""
    
    # Default pipeline
    _DEFAULT_STEPS = [
        'source',  # Unified data source (backend='sql' or 'api')
        'chembl',
        'curate',
        'tokens',
        'distributions',
    ]
    
    def __post_init__(self):
        """Initialize and validate parameters after dataclass creation."""
        if self.steps is None:
            self.steps = self._DEFAULT_STEPS.copy()
        
        super().__post_init__()
    
    def _validate_params(self) -> None:
        """Handle parameter initialization, overrides, and pipeline validation."""
        # Process all actors in the pipeline
        for step_name in self.steps:
            actor_info = ActorRegistry.get_actor_info(step_name)
            if not actor_info:
                continue

            if actor_info.get('is_plugin', False):
                self._initialize_plugin_params(step_name, actor_info)
            else:
                self._initialize_core_params(step_name, actor_info)

        # Validate pipeline configuration
        is_valid, error_message = ActorRegistry.validate_steps(self.steps)
        if not is_valid:
            raise ValueError(f"Invalid pipeline configuration: {error_message}")

        self._warn_vocab_antipattern()
        self._warn_step_position()

    def _warn_step_position(self) -> None:
        """
        Warn when a position-constrained step is misplaced.

        An actor may set ``__initial__ = True`` (expected to be the first step)
        or ``__terminal__ = True`` (expected to be the last step).
        """
        last = len(self.steps) - 1
        for i, step in enumerate(self.steps):
            info = ActorRegistry.get_actor_info(step)
            if not info:
                continue
            actor = info['class']
            if getattr(actor, '__initial__', False) and i != 0:
                warnings.warn(
                    f"Step '{step}' is an initial step and is expected to be the first "
                    f"step in the pipeline, but is preceded by {self.steps[:i]}.",
                    UserWarning,
                )
            if getattr(actor, '__terminal__', False) and i != last:
                warnings.warn(
                    f"Step '{step}' is a terminal step and is expected to be the last "
                    f"step in the pipeline, but is followed by {self.steps[i + 1:]}.",
                    UserWarning,
                )

    def _warn_vocab_antipattern(self) -> None:
        """
        Warn when a provided vocabulary is neither enforced nor extended.

        With a fixed vocabulary (``vocab_file`` set, ``dynamically_update_vocab=False``)
        and ``filter_unknown_tokens=False``, out-of-vocabulary tokens are neither
        filtered nor added, so the curated vocabulary may diverge from the base
        vocabulary.
        """
        if 'tokens' not in self.steps or 'distributions' not in self.steps:
            return
        tp, dp = self.tokens_params, self.distributions_params
        if tp is None or dp is None:
            return
        if (getattr(tp, 'vocab_file', None)
                and not getattr(tp, 'dynamically_update_vocab', True)
                and not getattr(dp, 'filter_unknown_tokens', True)):
            warnings.warn(
                "A fixed vocabulary is provided (vocab_file set, "
                "dynamically_update_vocab=False) while filter_unknown_tokens=False: "
                "out-of-vocabulary tokens are neither filtered nor added, so the "
                "curated vocabulary may diverge from the base vocabulary. Enable "
                "filter_unknown_tokens to enforce the vocabulary, or "
                "dynamically_update_vocab to extend it.",
                UserWarning,
            )

    def _initialize_core_params(self, step_name: str, actor_info: Dict) -> None:
        """Initialize parameters for core actors."""
        param_attr = ActorRegistry.get_param_attr(step_name)
        param_class = actor_info['param_class']

        # Create default params if not provided
        if getattr(self, param_attr) is None:
            setattr(self, param_attr, param_class())

        # Convert dict to param class, keeping only declared fields
        if isinstance(getattr(self, param_attr), dict):
            setattr(self, param_attr, ConfigCodec.build_param_object(param_class, getattr(self, param_attr)))

        # Override with pipe-level params
        if self.override_actor_params:
            self._override_params(getattr(self, param_attr))
    
    def _initialize_plugin_params(self, step_name: str, actor_info: Dict) -> None:
        """Initialize parameters for plugin actors."""
        param_class = actor_info['param_class']
        
        # Get or create plugin params
        if step_name not in self.plugin_params:
            self.plugin_params[step_name] = param_class()
        # Convert dict to param class, keeping only declared fields
        elif isinstance(self.plugin_params[step_name], dict):
            self.plugin_params[step_name] = ConfigCodec.build_param_object(param_class, self.plugin_params[step_name])
        
        # Override with pipe-level params
        if self.override_actor_params:
            self._override_params(self.plugin_params[step_name])
                
    def _override_params(self, actor_params) -> None:
        """Override actor parameters with pipe-level parameters."""
        pipe_fields = {f.name for f in fields(self)}
        actor_fields = {f.name for f in fields(actor_params)}
        
        # find common fields and override
        common_fields = pipe_fields & actor_fields
        for field_name in common_fields:
            setattr(actor_params, field_name, getattr(self, field_name))
    
        # Apply the same overrides to any nested params dataclass.
        for actor_field in fields(actor_params):
            nested_param = getattr(actor_params, actor_field.name)

            if hasattr(nested_param, '__dataclass_fields__'):
                self._override_params(nested_param)  # recursive

    def get_actor_params(self, step_name: str) -> Optional[Any]:
        """
        Get parameters for any actor (core or plugin).
        
        Args:
            step_name: Name of the pipeline step
            
        Returns:
            Parameter object for the actor, or None if not found
        """
        actor_info = ActorRegistry.get_actor_info(step_name)
        if not actor_info:
            return None
                
        if actor_info.get('is_plugin', False):
            return self.plugin_params.get(step_name)
        else:
            return getattr(self, ActorRegistry.get_param_attr(step_name), None)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary, including plugin params."""
        result = super().to_dict() if hasattr(super(), 'to_dict') else {}
        
        # Add all fields
        for f in fields(self):
            value = getattr(self, f.name)
            if hasattr(value, 'to_dict'):
                result[f.name] = value.to_dict()
            elif isinstance(value, dict):
                # Handle plugin_params dict
                result[f.name] = {
                    k: v.to_dict() if hasattr(v, 'to_dict') else v
                    for k, v in value.items()
                }
            else:
                result[f.name] = value
        
        return result
    
    def to_config(self) -> Dict[str, Any]:
        """Serialize to the canonical, flat configuration dict (see ConfigCodec)."""
        return ConfigCodec.to_config(self)

    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'PipeParams':
        """Build a params object from a canonical configuration dict (see ConfigCodec)."""
        return ConfigCodec.from_config(config, cls)

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipeParams':
        """Build a params object from a configuration dict (see from_config)."""
        return ConfigCodec.from_config(config_dict, cls)

    def to_json(self, filepath: str) -> None:
        """Save the configuration to a JSON file (canonical flat format)."""
        ConfigCodec.to_json(self, filepath)

    def to_yaml(self, filepath: str) -> None:
        """Save the configuration to a YAML file (canonical flat format)."""
        ConfigCodec.to_yaml(self, filepath)

    @classmethod
    def from_json(cls, config_path: str) -> 'PipeParams':
        """Load a configuration from a JSON file."""
        return ConfigCodec.from_json(config_path, cls)

    @classmethod
    def from_yaml(cls, config_path: str) -> 'PipeParams':
        """Load a configuration from a YAML file."""
        return ConfigCodec.from_yaml(config_path, cls)

    def _get_config_hash(self, length: int = 8) -> str:
        """Hash the configuration, bound to the MolForge version (see ConfigCodec)."""
        return ConfigCodec.config_hash(self, length)
    
# Clean user-facing API
class ForgeParams(PipeParams):
    """
    User-friendly alias for PipeParams.
    
    Use this when creating configurations for MolForge.
    """
    pass