"""
Centralized (de)serialization and hashing for MolForge configurations.

``ConfigCodec`` is the single source for reading and writing pipeline
configurations. The pipeline, CLI, and public API all borrow it, so every entry
point shares one canonical format, one run-metadata policy, and one
version-bound hash.

Canonical format: a flat, step-keyed mapping. Pipeline-level settings are
top-level keys; each configured step's parameters live under a block keyed by the
step name, with an identical shape for core actors and plugins. Private
attributes (keys starting with ``_``) are never written, and the MolForge version
is recorded and folded into the hash.
"""

from __future__ import annotations

import hashlib
import json
import types
import warnings
from dataclasses import fields, is_dataclass
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Dict, Type, Union, get_args, get_origin, get_type_hints

from .registry import ActorRegistry


class ConfigCodec:
    """Reads, writes, and hashes pipeline configurations in the canonical format."""

    # Run-specific metadata recorded in a written config at run time. Ignored when
    # a params object is built from a config, and excluded from the hash.
    RUN_METADATA_KEYS = frozenset({
        'run_id', 'log_time', 'input_id', 'input_type', 'input_source',
    })

    # Excluded from the configuration hash: run metadata plus output and logging
    # settings, which do not change the scientific configuration.
    NON_CONFIG_KEYS = RUN_METADATA_KEYS | frozenset({
        'output_root', 'write_output', 'write_checkpoints', 'save_config',
        'verbose', 'console_log_level', 'file_log_level',
    })

    # ------------------------------------------------------------------ encode
    @classmethod
    def to_config(cls, params) -> Dict[str, Any]:
        """Serialize a params object to the canonical flat, step-keyed dict."""
        from .._version import __version__ as molforge_version

        core_attrs = {ActorRegistry.get_param_attr(s) for s in ActorRegistry.get_core_steps()}
        config: Dict[str, Any] = {'version': molforge_version}

        # Pipeline-level fields (exclude the per-actor param containers).
        for f in fields(params):
            if f.name.startswith('_') or f.name in core_attrs or f.name == 'plugin_params':
                continue
            config[f.name] = getattr(params, f.name)

        # Per-step blocks, keyed by step name (core and plugin alike).
        for step in params.steps:
            params_obj = params.get_actor_params(step)
            if params_obj is not None:
                config[step] = cls._to_plain(params_obj)

        return config

    # ------------------------------------------------------------------ decode
    @classmethod
    def from_config(cls, config: Dict[str, Any], params_cls: Type) -> Any:
        """
        Build a ``params_cls`` instance from a canonical configuration dict.

        Per-step blocks are routed to the matching core field or the
        ``plugin_params`` map via the registry, keeping only declared fields so
        attributes derived in ``__post_init__`` are re-derived on construction.
        The version and run-specific metadata are ignored; any remaining
        unrecognized key is dropped with a warning. Both the flat, step-keyed
        form and the field-keyed form (``source_params``, ``plugin_params``) are
        accepted.
        """
        config = dict(config or {})
        config.pop('version', None)
        for key in cls.RUN_METADATA_KEYS:
            config.pop(key, None)
        cls._expand_output(config)

        pipe_field_names = {f.name for f in fields(params_cls)}
        kwargs: Dict[str, Any] = {}
        plugin_params: Dict[str, Any] = {}
        unknown: list = []

        for key, value in config.items():
            info = ActorRegistry.get_actor_info(key)
            if info is not None:
                if isinstance(value, dict):
                    param_obj = cls.build_param_object(info['param_class'], value)
                    if info.get('is_plugin', False):
                        plugin_params[key] = param_obj
                    else:
                        kwargs[ActorRegistry.get_param_attr(key)] = param_obj
            elif key in pipe_field_names:
                kwargs[key] = value
            else:
                unknown.append(key)

        if unknown:
            cls._warn_unknown_keys(unknown, pipe_field_names)
        if plugin_params:
            kwargs['plugin_params'] = plugin_params

        return params_cls(**kwargs)

    # ---------------------------------------------------------------- file I/O
    @classmethod
    def to_json(cls, params, filepath: str) -> None:
        """Write the configuration to a JSON file."""
        with open(filepath, 'w') as f:
            json.dump(cls.to_config(params), f, indent=4, default=str)

    @classmethod
    def to_yaml(cls, params, filepath: str) -> None:
        """Write the configuration to a YAML file."""
        import yaml
        with open(filepath, 'w') as f:
            yaml.safe_dump(cls.to_config(params), f, sort_keys=False, default_flow_style=False)

    @classmethod
    def from_json(cls, config_path: str, params_cls: Type) -> Any:
        """Load a configuration from a JSON file."""
        with open(config_path, 'r') as f:
            data = json.load(f)
        return cls.from_config(cls._require_mapping(data, config_path), params_cls)

    @classmethod
    def from_yaml(cls, config_path: str, params_cls: Type) -> Any:
        """Load a configuration from a YAML file."""
        import yaml
        with open(config_path, 'r') as f:
            data = yaml.safe_load(f)
        return cls.from_config(cls._require_mapping(data, config_path), params_cls)

    @classmethod
    def from_file(cls, config_path: str, params_cls: Type) -> Any:
        """Load a configuration from a JSON or YAML file, chosen by extension."""
        suffix = Path(config_path).suffix.lower()
        if suffix in ('.yaml', '.yml'):
            return cls.from_yaml(config_path, params_cls)
        if suffix == '.json':
            return cls.from_json(config_path, params_cls)
        raise ValueError(f"Unsupported config format: {suffix}")

    # -------------------------------------------------------------------- hash
    @classmethod
    def config_hash(cls, params, length: int = 8) -> str:
        """
        Hash the configuration, bound to the MolForge version.

        Output, logging, and run-specific metadata are excluded so the hash
        identifies the configuration itself; the version is recorded in the
        canonical dict, so a hash match implies a version match.
        """
        def strip(d):
            if not isinstance(d, dict):
                return d
            return {k: strip(v) for k, v in d.items() if k not in cls.NON_CONFIG_KEYS}

        payload = strip(cls.to_config(params))
        digest = json.dumps(payload, sort_keys=True, default=str)
        return hashlib.sha256(digest.encode()).hexdigest()[:length].upper()

    # ----------------------------------------------------------------- helpers
    @classmethod
    def _to_plain(cls, obj: Any) -> Any:
        """Recursively convert params objects and containers to plain data.

        Uses each params object's ``to_dict`` (which drops private ``_`` keys), so
        private and run-time attributes never reach the serialized form.
        """
        if hasattr(obj, 'to_dict'):
            return {k: cls._to_plain(v) for k, v in obj.to_dict().items()}
        if isinstance(obj, dict):
            return {k: cls._to_plain(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [cls._to_plain(v) for v in obj]
        return obj

    @classmethod
    def build_param_object(cls, param_class, data: Dict[str, Any]):
        """Construct a params dataclass from a config block.

        Constructor fields are populated from the block, and a field whose type is
        itself a params dataclass is rebuilt recursively at any depth. Computed
        ``init=False`` fields are left for the class to derive. A key that matches
        no field is an unknown parameter and raises, so typos surface at load time.
        """
        hints = cls._type_hints(param_class)
        flds = fields(param_class)
        field_names = {f.name for f in flds}
        init_names = {f.name for f in flds if f.init}

        unknown = set(data) - field_names
        if unknown:
            raise TypeError(
                f"Unknown parameter(s) {sorted(unknown)} for {param_class.__name__}. "
                f"Valid parameters: {sorted(init_names)}."
            )

        kwargs: Dict[str, Any] = {}
        for name in init_names & set(data):
            value = data[name]
            nested = cls._dataclass_type(hints.get(name))
            if nested is not None and isinstance(value, dict):
                value = cls.build_param_object(nested, value)
            kwargs[name] = value
        return param_class(**kwargs)

    @staticmethod
    def _type_hints(param_class) -> Dict[str, Any]:
        """Resolve a class's type hints, tolerating unresolvable annotations."""
        try:
            return get_type_hints(param_class)
        except Exception:
            return {}

    @staticmethod
    def _dataclass_type(hint):
        """Return the params dataclass behind a hint, unwrapping ``Optional``.

        Handles both ``Optional[X]``/``Union[X, None]`` and the PEP 604 ``X | None``
        syntax, so a plugin may annotate a nested-params field either way.
        """
        if is_dataclass(hint):
            return hint
        if get_origin(hint) in (Union, types.UnionType):
            for arg in get_args(hint):
                if is_dataclass(arg):
                    return arg
        return None

    @staticmethod
    def _expand_output(config: Dict[str, Any]) -> None:
        """Expand the grouped ``output: {dir, checkpoints}`` shorthand in place.

        Explicit flat keys take precedence, so a config may use either form.
        """
        output = config.pop('output', None)
        if isinstance(output, dict):
            if 'dir' in output:
                config.setdefault('output_root', output['dir'])
            if 'checkpoints' in output:
                config.setdefault('write_checkpoints', output['checkpoints'])

    @staticmethod
    def _warn_unknown_keys(unknown, pipe_field_names) -> None:
        """Warn that unrecognized keys were ignored, suggesting close matches."""
        known = set(pipe_field_names) | set(ActorRegistry.get_all_step_names())
        parts = []
        for key in unknown:
            match = get_close_matches(str(key), known, n=1)
            parts.append(f"'{key}' (did you mean '{match[0]}'?)" if match else f"'{key}'")
        warnings.warn(
            f"Ignoring unrecognized configuration key(s): {', '.join(parts)}.",
            stacklevel=2,
        )

    @staticmethod
    def _require_mapping(data, source: str):
        """Reject an empty or non-mapping configuration with a clear message."""
        if not data:
            raise ValueError(
                f"Configuration '{source}' is empty. Provide a configuration with "
                f"at least one pipeline step."
            )
        if not isinstance(data, dict):
            raise ValueError(
                f"Configuration '{source}' must be a mapping at the top level."
            )
        return data
