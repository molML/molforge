"""
Shared introspection layer for MolForge parameter dataclasses.

This module is pure Python (no CLI dependencies) so that BOTH the CLI and a
future GUI can discover, describe and validate the fields of any actor's
parameter dataclass through ONE shared, introspectable core.

It provides:
    - ``attribute_docs``: extract PEP 257 attribute docstrings at runtime.
    - ``describe_param_class``: structured metadata for each dataclass field.
    - ``describe_step``: full description of a pipeline step (params + backends).
    - ``list_steps``: the available core and plugin steps.

Nothing here imports questionary, argparse or any other UI machinery.
"""

from __future__ import annotations

import ast
import inspect
import importlib
import textwrap
import typing
from dataclasses import MISSING, fields, is_dataclass
from typing import Any, Dict, List, Optional, Type


__all__ = [
    "attribute_docs",
    "describe_param_class",
    "describe_step",
    "list_steps",
]


# ---------------------------------------------------------------------------
# Attribute docstrings (PEP 257 style: a string literal directly under a field)
# ---------------------------------------------------------------------------

def attribute_docs(cls: Type) -> Dict[str, str]:
    """
    Extract PEP 257 attribute docstrings for a dataclass.

    Parses the source of ``cls`` (and every class in its MRO) with ``ast`` and
    collects the string literal that appears *immediately* under each
    field/assignment. This makes per-field help available at runtime even
    though Python does not otherwise retain attribute docstrings.

    Inheritance is handled by walking the MRO from base to derived, so a
    subclass docstring overrides a base one and inherited fields (e.g.
    ``BaseParams.verbose``) are still picked up.

    Robust against unavailable source (returns ``{}`` / partial results).

    Args:
        cls: A dataclass (or any class) to introspect.

    Returns:
        Mapping of ``field_name -> docstring``.
    """
    docs: Dict[str, str] = {}

    try:
        mro = inspect.getmro(cls)
    except Exception:
        return docs

    # base -> derived so more-derived classes override inherited docs
    for klass in reversed(mro):
        if klass is object:
            continue

        try:
            source = inspect.getsource(klass)
        except (OSError, TypeError):
            continue

        try:
            tree = ast.parse(textwrap.dedent(source))
        except SyntaxError:
            continue

        class_def = next(
            (node for node in tree.body if isinstance(node, ast.ClassDef)),
            None,
        )
        if class_def is None:
            continue

        body = class_def.body
        for index, node in enumerate(body):
            name = _assignment_target_name(node)
            if name is None:
                continue

            if index + 1 >= len(body):
                continue

            following = body[index + 1]
            if (
                isinstance(following, ast.Expr)
                and isinstance(following.value, ast.Constant)
                and isinstance(following.value.value, str)
            ):
                docs[name] = inspect.cleandoc(following.value.value)

    return docs


def _assignment_target_name(node: ast.stmt) -> Optional[str]:
    """Return the single target name for an assignment node, else ``None``."""
    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
        return node.target.id
    if (
        isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
    ):
        return node.targets[0].id
    return None


# ---------------------------------------------------------------------------
# Type analysis
# ---------------------------------------------------------------------------

def _type_to_str(tp: Any) -> str:
    """Produce a compact, human-readable type name."""
    if tp is type(None):
        return "None"
    if isinstance(tp, str):
        return tp
    if isinstance(tp, type):
        return tp.__name__
    return str(tp).replace("typing.", "")


def _analyze_type(tp: Any) -> Dict[str, Any]:
    """
    Analyze a type annotation into introspection metadata.

    Returns a dict with keys: ``type`` (str), ``choices`` (list|None),
    ``is_bool`` (bool), ``is_optional`` (bool).
    """
    is_optional = False

    # Forward reference / string annotation: nothing more we can infer.
    if isinstance(tp, str):
        return {"type": tp, "choices": None, "is_bool": False, "is_optional": False}

    origin = typing.get_origin(tp)

    # Unwrap Optional / Union[..., None]
    if origin is typing.Union:
        args = typing.get_args(tp)
        none_type = type(None)
        non_none = [a for a in args if a is not none_type]
        if len(non_none) < len(args):
            is_optional = True
        if len(non_none) == 1:
            tp = non_none[0]
            origin = typing.get_origin(tp)
        else:
            # e.g. Union[List[str], str, None] -> no single choice/bool
            return {
                "type": _type_to_str(tp),
                "choices": None,
                "is_bool": False,
                "is_optional": is_optional,
            }

    if origin is typing.Literal:
        return {
            "type": _type_to_str(tp),
            "choices": list(typing.get_args(tp)),
            "is_bool": False,
            "is_optional": is_optional,
        }

    return {
        "type": _type_to_str(tp),
        "choices": None,
        "is_bool": tp is bool,
        "is_optional": is_optional,
    }


def _field_default(f) -> Any:
    """Resolve the concrete default for a dataclass field (or ``None``)."""
    if f.default is not MISSING:
        return f.default
    if f.default_factory is not MISSING:  # type: ignore[comparison-overlap]
        try:
            return f.default_factory()
        except Exception:
            return None
    return None


def _field_metadata(f, doc: str = "") -> Dict[str, Any]:
    """Build the metadata dict for a single dataclass field."""
    info = _analyze_type(f.type)
    return {
        "name": f.name,
        "type": info["type"],
        "default": _field_default(f),
        "choices": info["choices"],
        "is_bool": info["is_bool"],
        "is_optional": info["is_optional"],
        "doc": doc,
    }


def describe_param_class(cls: Type) -> List[Dict[str, Any]]:
    """
    Describe every public field of a dataclass.

    Args:
        cls: A dataclass parameter class.

    Returns:
        A list of dicts, one per field (skipping names starting with ``_``),
        each with keys: ``name``, ``type``, ``default``, ``choices``,
        ``is_bool``, ``is_optional``, ``doc``.
    """
    if cls is None or not is_dataclass(cls):
        return []

    docs = attribute_docs(cls)
    result: List[Dict[str, Any]] = []
    for f in fields(cls):
        if f.name.startswith("_"):
            continue
        result.append(_field_metadata(f, docs.get(f.name, "")))
    return result


# ---------------------------------------------------------------------------
# Step / pipeline description
# ---------------------------------------------------------------------------

def _step_backends(step_name: str) -> List[str]:
    """Return the list of backend names registered for a step (may be empty)."""
    from ..backends.registry import BackendRegistry

    # Backends self-register on import of molforge.backends.<step>; trigger it.
    try:
        importlib.import_module(f"molforge.backends.{step_name}")
    except ModuleNotFoundError:
        pass
    except Exception:
        pass

    return BackendRegistry.list_backends(step_name).get(step_name, [])


def describe_step(step_name: str) -> Dict[str, Any]:
    """
    Fully describe a pipeline step.

    Args:
        step_name: Registered step name (core or plugin).

    Returns:
        Dict with keys: ``step``, ``is_plugin``, ``param_class``, ``fields``
        (from ``describe_param_class``), and ``backends`` (possibly empty).

    Raises:
        ValueError: If the step is not registered.
    """
    from .registry import ActorRegistry

    info = ActorRegistry.get_actor_info(step_name)
    if info is None:
        available = ", ".join(ActorRegistry.get_all_step_names())
        raise ValueError(f"Unknown step '{step_name}'. Available: {available}")

    param_class = ActorRegistry.get_param_class(step_name)

    return {
        "step": step_name,
        "is_plugin": bool(info.get("is_plugin", False)),
        "param_class": param_class,
        "fields": describe_param_class(param_class) if param_class else [],
        "backends": _step_backends(step_name),
    }


def list_steps() -> Dict[str, List[str]]:
    """
    List the available core and plugin steps.

    Returns:
        ``{'core': [...], 'plugin': [...]}``.
    """
    from .registry import ActorRegistry

    return {
        "core": ActorRegistry.get_core_steps(),
        "plugin": ActorRegistry.get_plugin_steps(),
    }
