"""
Clean parameter system supporting core actors (typed) and plugins (dynamic).
"""

from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field, fields
import hashlib
import json
import copy

from ..actors.params import *
from .registry import ActorRegistry

@dataclass
class PipeParams(BaseParams):
    """
    Configuration parameters for MolForge Pipeline.
    
    Supports both core actors (with explicit type hints) and plugin actors
    (with dynamic parameter handling).
    """
    
    # Pipeline configuration
    steps: List[str] = None

    # Pipeline behavior
    return_none_on_fail: bool = False
    output_root: str = "MOLFORGE_OUT"
    write_output: bool = True
    write_checkpoints: bool = False
    
    # Logging
    console_log_level: str = 'INFO'
    file_log_level: str = 'DEBUG'
    
     # Parameters
    override_actor_params: bool = True  # use pipe params instead of actor params when both exist

    # Core actors (explicit, type-hinted)
    source_params: Optional[ChEMBLSourceParams] = None
    chembl_params: Optional[ChEMBLCuratorParams] =  None
    curate_params: Optional[CurateMolParams]     =  None
    tokens_params: Optional[TokenizeDataParams]  =  None
    distributions_params: Optional[CurateDistributionParams]  =  None
    confs_params:  Optional[GenerateConfsParams] =  None
    
    # Plugin actors (dynamic, no type hints)
    plugin_params: Dict[str, Any] = field(default_factory=dict)
    
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
    
    def _initialize_core_params(self, step_name: str, actor_info: Dict) -> None:
        """Initialize parameters for core actors."""
        param_attr = ActorRegistry.get_param_attr(step_name)
        param_class = actor_info['param_class']

        # Create default params if not provided
        if getattr(self, param_attr) is None:
            setattr(self, param_attr, param_class())

        # Convert dict to param class
        if isinstance(getattr(self, param_attr), dict):
            setattr(self, param_attr, param_class(**getattr(self, param_attr)))

        # Override with pipe-level params
        if self.override_actor_params:
            self._override_params(getattr(self, param_attr))
    
    def _initialize_plugin_params(self, step_name: str, actor_info: Dict) -> None:
        """Initialize parameters for plugin actors."""
        param_class = actor_info['param_class']
        
        # Get or create plugin params
        if step_name not in self.plugin_params:
            self.plugin_params[step_name] = param_class()
        # Convert dict to param class
        elif isinstance(self.plugin_params[step_name], dict):
            self.plugin_params[step_name] = param_class(**self.plugin_params[step_name])
        
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
    
        # also check nested params (mostly relevant for stereo_params in CurateMolParams)
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
            param_attr = actor_info['param_attr']
            return getattr(self, param_attr, None)
    
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
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'PipeParams':
        """
        Create PipeParams from dictionary.
        
        Handles both core and plugin actor parameters with recursive nested dataclass conversion.
        """
        # remove config-specific entries that we don't want to leak into a new instance.
        # TODO: dynamically fix this so that we can log extra stuff without breaking the read method.
        for to_remove in ['log_time',  'input_type', 'input_id', 'input_source', 'run_id']:
            config_dict.pop(to_remove, None)
        
        # convert dictionaries to dataclasses
        for step_name in ActorRegistry.get_all_step_names():
            actor_info = ActorRegistry.get_actor_info(step_name)
            if not actor_info:
                continue
                
            param_attr = actor_info['param_attr']
            param_class = actor_info['param_class']
            
            if param_attr in config_dict and config_dict[param_attr] is not None:
                # retrieve dataclass fields
                valid_fields = {f.name for f in fields(param_class)}
                # ignore class-init level attributes                
                filtered_data = {k: v for k, v in config_dict[param_attr].items() if k in valid_fields}
                # assign back to config
                config_dict[param_attr] = param_class(**filtered_data)
        
        return cls(**config_dict)
    
    @classmethod
    def from_json(cls, config_path: str) -> 'PipeParams':
        """Load configuration from JSON file."""
        with open(config_path, 'r') as f:
            config_dict = json.load(f)
        
        return cls.from_dict(config_dict)
    
    def _get_params_dict(self) -> dict:
        """Convert params to nested dictionary."""
        def to_dict(obj):
            if hasattr(obj, 'to_dict'):
                result = obj.to_dict()
            elif isinstance(obj, dict):
                result = copy.deepcopy(obj)
            else:
                return obj
            
            # rcursively convert nested objects.
            return {k: to_dict(v) for k, v in result.items()}
        
        return to_dict(self)
    
    def to_json(self, filepath: str) -> None:
        """Save configuration to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self._get_params_dict(), f, indent=4, default=str)

    def _get_config_hash(self, length: int = 8) -> str:
        """Generate unique hash from configuration."""
        # ignore unimportant params for hashing.
        exclude = {'output_root', 'write_output', 'write_checkpoints',
                   'save_config', 'verbose', 'console_log_level', 'file_log_level'}
        
        def filter_dict(d):
            if not isinstance(d, dict):
                return d
            # filter out recursively, consider e.g., verbose exists in actors too.
            return {k: filter_dict(v) for k, v in d.items() if k not in exclude}
        
        config = filter_dict(self._get_params_dict())
        config_str = json.dumps(config, sort_keys=True, default=str)
        return hashlib.sha256(config_str.encode()).hexdigest()[:length].upper()
    
# Clean user-facing API
class ForgeParams(PipeParams):
    """
    User-friendly alias for PipeParams.
    
    Use this when creating configurations for MolForge.
    """
    pass