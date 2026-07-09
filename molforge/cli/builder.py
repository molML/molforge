"""
Shared parameter-assembly point for MolForge.

``build_params`` is the SINGLE place where a ``ForgeParams`` object is assembled
from any combination of:
    - a loaded config file (YAML/JSON),
    - an explicit list of steps,
    - a list of ``key=value`` / ``step.key=value`` overrides,
    - an output root.

Both the CLI (``run`` / wizard) and a future GUI funnel through this function so
that validation, type coercion and step routing behave identically everywhere.

The chosen steps are placed into the params dict before construction, so
``ForgeParams.__post_init__`` initializes and validates the params of every
selected step. ``--steps`` and ``--set`` overrides are layered on top of a loaded
``--config``. Assembly finishes through :class:`ConfigCodec`, so step blocks are
routed and filtered exactly as every other entry point.

Depends only on pyyaml (an existing ``cli`` extra) for file loading.
"""

from __future__ import annotations

import json
from dataclasses import fields
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ..configuration.params import ForgeParams
from ..configuration.registry import ActorRegistry
from ..configuration.serialization import ConfigCodec
from ..configuration import introspect


__all__ = ["build_params", "parse_override"]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_params(
    steps: Optional[Union[str, List[str]]] = None,
    config_path: Optional[str] = None,
    overrides: Optional[List[str]] = None,
    output_root: Optional[str] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> ForgeParams:
    """
    Assemble a validated ``ForgeParams`` from all supported sources.

    Args:
        steps: Comma-separated string or list of step names. Overrides any
            ``steps`` loaded from ``config_path``.
        config_path: Optional YAML/JSON configuration file to seed the params.
        overrides: List of ``"key=value"`` or ``"step.key=value"`` strings.
            Values are coerced to the target field's type using introspection
            metadata and applied on top of the config/steps.
        output_root: Optional output root directory (top-level override).
        extra: Optional dict merged into the params before construction (lowest
            precedence after config, useful for programmatic callers).

    Returns:
        A fully constructed and validated ``ForgeParams`` instance.
    """
    merged: Dict[str, Any] = {}

    if config_path:
        merged.update(_read_config_file(config_path))

    if extra:
        merged.update(extra)

    if steps is not None:
        merged["steps"] = _normalize_steps(steps)

    if output_root is not None:
        merged["output_root"] = output_root

    if overrides:
        for override in overrides:
            _apply_override(merged, override)

    # Assemble through the shared codec so the output shorthand, step routing, and
    # field filtering behave exactly as every other entry point.
    return ConfigCodec.from_config(merged, ForgeParams)


def parse_override(override: str) -> tuple[Optional[str], str, str]:
    """
    Split an override string into ``(step, field, raw_value)``.

    ``step`` is ``None`` for a bare ``key=value`` (top-level) override.

    Raises:
        ValueError: If the string is not of the form ``key=value``.
    """
    if "=" not in override:
        raise ValueError(
            f"Invalid override '{override}'. Expected KEY=VALUE or step.KEY=VALUE."
        )
    key, _, raw_value = override.partition("=")
    key = key.strip()
    raw_value = raw_value.strip()
    if not key:
        raise ValueError(f"Invalid override '{override}': empty key.")

    if "." in key:
        step_name, field_name = key.split(".", 1)
        return step_name.strip(), field_name.strip(), raw_value
    return None, key, raw_value


# ---------------------------------------------------------------------------
# Config loading / normalization
# ---------------------------------------------------------------------------

def _read_config_file(config_path: str) -> Dict[str, Any]:
    """Read a YAML/JSON config file into a raw dict (no normalization)."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    suffix = path.suffix.lower()
    if suffix in (".yaml", ".yml"):
        import yaml  # local import so builder is importable without a YAML load

        with open(path, "r") as f:
            raw = yaml.safe_load(f)
    elif suffix == ".json":
        with open(path, "r") as f:
            raw = json.load(f)
    else:
        raise ValueError(f"Unsupported config format: {path.suffix}")

    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError("Configuration file must contain a mapping at the top level.")

    return raw


def _normalize_steps(steps: Union[str, List[str]]) -> List[str]:
    """Turn a comma-separated string or iterable into a clean step list."""
    if isinstance(steps, str):
        return [s.strip() for s in steps.split(",") if s.strip()]
    return [str(s).strip() for s in steps if str(s).strip()]


# ---------------------------------------------------------------------------
# Override application + type coercion
# ---------------------------------------------------------------------------

def _apply_override(merged: Dict[str, Any], override: str) -> None:
    """Parse a single override string and route it into ``merged``."""
    step_name, field_name, raw_value = parse_override(override)

    if step_name is None:
        _apply_toplevel_override(merged, field_name, raw_value)
    else:
        _apply_step_override(merged, step_name, field_name, raw_value)


def _apply_toplevel_override(merged: Dict[str, Any], field_name: str, raw_value: str) -> None:
    """Coerce and set a bare ``key=value`` onto the top-level ForgeParams field."""
    forge_fields = {f.name: f for f in fields(ForgeParams)}
    if field_name not in forge_fields:
        raise ValueError(
            f"Unknown parameter '{field_name}'. "
            f"Use 'step.{field_name}=...' to target a specific step."
        )
    meta = introspect._field_metadata(forge_fields[field_name])
    merged[field_name] = _coerce(raw_value, meta)


def _apply_step_override(
    merged: Dict[str, Any], step_name: str, field_name: str, raw_value: str
) -> None:
    """Coerce and route a ``step.key=value`` override into the step's params."""
    param_class = ActorRegistry.get_param_class(step_name)
    if param_class is None:
        available = ", ".join(ActorRegistry.get_all_step_names())
        raise ValueError(f"Unknown step '{step_name}' in override. Available: {available}")

    value = _coerce_field(param_class, field_name, raw_value)

    # Route into the flat, step-keyed block; the codec maps it to the core field
    # or the plugin_params map.
    block = merged.get(step_name)
    if block is None:
        block = {}
        merged[step_name] = block
    if not isinstance(block, dict):
        raise ValueError(
            f"Cannot apply override to '{step_name}.{field_name}': "
            f"existing config for '{step_name}' is not a mapping."
        )
    block[field_name] = value


def _coerce_field(param_class, field_name: str, raw_value: str) -> Any:
    """Coerce ``raw_value`` to the type of ``field_name`` on ``param_class``."""
    for meta in introspect.describe_param_class(param_class):
        if meta["name"] == field_name:
            return _coerce(raw_value, meta)
    valid = ", ".join(m["name"] for m in introspect.describe_param_class(param_class))
    raise ValueError(
        f"Unknown parameter '{field_name}' for {param_class.__name__}. Valid: {valid}"
    )


def _coerce(raw_value: str, meta: Dict[str, Any]) -> Any:
    """Coerce a raw string to the target type described by ``meta``."""
    # Optional fields accept explicit null-ish values.
    if meta.get("is_optional") and raw_value.lower() in ("none", "null"):
        return None

    choices = meta.get("choices")
    if choices is not None:
        return _match_choice(raw_value, choices)

    if meta.get("is_bool"):
        return _to_bool(raw_value)

    return _coerce_by_typestr(raw_value, meta.get("type", "str"))


def _coerce_by_typestr(raw_value: str, type_str: str) -> Any:
    """Coerce using a best-effort interpretation of the type name string."""
    base = type_str.strip()

    if base == "bool":
        return _to_bool(raw_value)
    if base == "int":
        try:
            return int(raw_value)
        except ValueError:
            raise ValueError(f"Expected an integer, got '{raw_value}'.")
    if base == "float":
        try:
            return float(raw_value)
        except ValueError:
            raise ValueError(f"Expected a number, got '{raw_value}'.")
    if base.startswith(("List", "list", "Sequence", "Tuple", "tuple")):
        return _parse_list(raw_value)
    if base.startswith(("Dict", "dict", "Mapping")):
        return _parse_json(raw_value)

    # str, Optional[str], unknown -> keep as string
    return raw_value


def _parse_list(raw_value: str) -> Any:
    """Parse a list from JSON or a comma-separated string."""
    stripped = raw_value.strip()
    if stripped.startswith("["):
        return _parse_json(stripped)
    if stripped == "":
        return []
    return [item.strip() for item in stripped.split(",") if item.strip()]


def _parse_json(raw_value: str) -> Any:
    try:
        return json.loads(raw_value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Expected JSON value, got '{raw_value}': {exc}")


def _to_bool(raw_value: str) -> bool:
    lowered = raw_value.strip().lower()
    if lowered in ("true", "1", "yes", "y", "on"):
        return True
    if lowered in ("false", "0", "no", "n", "off"):
        return False
    raise ValueError(f"Expected a boolean (true/false), got '{raw_value}'.")


def _match_choice(raw_value: str, choices: List[Any]) -> Any:
    """Validate ``raw_value`` against a Literal's allowed choices."""
    lowered = raw_value.strip().lower()

    for choice in choices:
        if isinstance(choice, bool):
            if lowered in ("true", "false", "1", "0", "yes", "no") and _to_bool(raw_value) == choice:
                return choice
        elif choice is None:
            if lowered in ("none", "null"):
                return None
        elif isinstance(choice, str):
            if choice == raw_value:
                return choice
        else:
            if str(choice) == raw_value:
                return choice

    valid = ", ".join(repr(c) for c in choices)
    raise ValueError(f"Invalid value '{raw_value}'. Allowed choices: {valid}.")
