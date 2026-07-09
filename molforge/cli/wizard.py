"""
Interactive configuration wizard for MolForge.

Guides a user through picking pipeline steps, tuning pipeline-level and per-step
parameters, then funnels everything through ``builder.build_params`` (the same
assembly point the non-interactive CLI uses). The interactive layer stays thin:
prompt helpers only gather strings and hand ``step.field=value`` / bare
``key=value`` overrides to ``build_params``, so this flow can later back a
persistent UI. ``questionary`` is imported lazily so the rest of the CLI works
without the optional ``cli`` extra installed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..configuration.params import ForgeParams
from ..configuration import introspect
from . import builder


__all__ = ["run_wizard", "questionary_available"]


_MISSING_MESSAGE = (
    "The interactive wizard requires 'questionary'.\n"
    "Install it with:  pip install questionary\n"
    "(or: pip install 'molforge[cli]')"
)

# Top-level ForgeParams fields offered under the pipeline-options gate.
_PIPELINE_FIELDS = [
    "output_root",
    "write_checkpoints",
    "console_log_level",
    "file_log_level",
    "return_none_on_fail",
]


def questionary_available() -> bool:
    """Return True if questionary can be imported."""
    try:
        import questionary  # noqa: F401
        return True
    except ImportError:
        return False


def _require_questionary():
    """Import questionary or raise a friendly RuntimeError."""
    try:
        import questionary
        return questionary
    except ImportError as exc:
        raise RuntimeError(_MISSING_MESSAGE) from exc


def run_wizard(preselected_steps: Optional[List[str]] = None) -> Optional[ForgeParams]:
    """
    Interactively assemble a ``ForgeParams``.

    Args:
        preselected_steps: Steps to check by default (defaults to the standard
            pipeline, ``ForgeParams._DEFAULT_STEPS``).

    Returns:
        A validated ``ForgeParams`` instance, or ``None`` if the user cancelled
        or no steps were selected. Ctrl+C at any prompt cancels the wizard.
    """
    questionary = _require_questionary()

    try:
        steps = _prompt_steps(questionary, preselected_steps)
        if not steps:
            print("No steps selected; nothing to configure.")
            return None

        overrides: List[str] = []
        overrides.extend(_prompt_pipeline_options(questionary))
        for step in steps:
            overrides.extend(_prompt_step(questionary, step))

        try:
            params = builder.build_params(steps=steps, overrides=overrides)
        except Exception as exc:
            print(f"Failed to build configuration: {exc}")
            return None

        _finalize(questionary, params, steps)
        return params
    except KeyboardInterrupt:
        print("\nWizard cancelled.")
        return None


# ---------------------------------------------------------------------------
# Step selection
# ---------------------------------------------------------------------------

def _prompt_steps(questionary, preselected_steps: Optional[List[str]]) -> List[str]:
    """Prompt the user to pick pipeline steps via a checkbox list."""
    available = introspect.list_steps()
    default_steps = set(preselected_steps or ForgeParams._DEFAULT_STEPS)

    choices = []
    for step in available["core"]:
        choices.append(questionary.Choice(step, checked=step in default_steps))
    for step in available["plugin"]:
        choices.append(
            questionary.Choice(f"{step} (plugin)", value=step, checked=step in default_steps)
        )

    answer = questionary.checkbox(
        "Select pipeline steps (space to toggle, enter to confirm):",
        choices=choices,
    ).unsafe_ask()

    return answer or []


# ---------------------------------------------------------------------------
# Parameter prompting
# ---------------------------------------------------------------------------

def _prompt_pipeline_options(questionary) -> List[str]:
    """
    Gate and prompt the top-level pipeline options.

    Returns bare ``key=value`` overrides, which ``build_params`` routes to the
    top-level ForgeParams fields.
    """
    metas = {m["name"]: m for m in introspect.describe_param_class(ForgeParams)}
    selected = [metas[name] for name in _PIPELINE_FIELDS if name in metas]
    if not selected:
        return []

    print(f"\nPipeline options: {_default_summary(selected)}")
    if not questionary.confirm("Configure pipeline options?", default=False).unsafe_ask():
        return []

    overrides: List[str] = []
    for meta in selected:
        override = _field_override(questionary, "", meta)
        if override is not None:
            overrides.append(override)
    return overrides


def _prompt_step(questionary, step: str) -> List[str]:
    """Gate and prompt a step's fields; return a list of override strings."""
    try:
        info = introspect.describe_step(step)
    except ValueError as exc:
        print(f"  (skipping unknown step '{step}': {exc})")
        return []

    fields_meta = info["fields"]
    if not fields_meta:
        print(f"\nStep '{step}': no configurable parameters.")
        return []

    print(f"\nStep '{step}' defaults: {_default_summary(fields_meta)}")
    if not questionary.confirm(f"Configure '{step}'?", default=False).unsafe_ask():
        return []

    overrides: List[str] = []
    for meta in fields_meta:
        override = _field_override(questionary, f"{step}.", meta)
        if override is not None:
            overrides.append(override)
    return overrides


def _default_summary(fields_meta: List[Dict[str, Any]], limit: int = 6) -> str:
    """Render a short ``name=default, ...`` summary of a field list.

    The internal ``verbose`` flag is omitted and the list is capped at ``limit``
    fields, with a count of any remaining fields.
    """
    fields = [m for m in fields_meta if m["name"] != "verbose"]
    shown = ", ".join(f"{m['name']}={m.get('default')}" for m in fields[:limit])
    remaining = len(fields) - limit
    return f"{shown} (+{remaining} more)" if remaining > 0 else shown


def _field_override(questionary, prefix: str, meta: Dict[str, Any]) -> Optional[str]:
    """
    Prompt for a single field and build an override string.

    ``prefix`` is ``"<step>."`` for a step field or ``""`` for a top-level
    field. Returns ``"<prefix><name>=<value>"`` when the answer differs from the
    default, or ``None`` to keep the default.
    """
    name = meta["name"]
    doc = meta.get("doc") or ""
    default = meta.get("default")

    # Literal choices -> select
    if meta.get("choices") is not None:
        choices = [str(c) for c in meta["choices"]]
        default_str = str(default) if default is not None else None
        answer = questionary.select(
            f"{name}:",
            choices=choices,
            default=default_str if default_str in choices else None,
            instruction=doc,
        ).unsafe_ask()
        if answer is None or answer == default_str:
            return None
        return f"{prefix}{name}={answer}"

    # Booleans -> confirm
    if meta.get("is_bool"):
        answer = questionary.confirm(
            f"{name}? ({doc})" if doc else f"{name}?",
            default=bool(default),
        ).unsafe_ask()
        if answer is None or answer == bool(default):
            return None
        return f"{prefix}{name}={'true' if answer else 'false'}"

    # Complex containers -> optional JSON text (blank keeps default)
    if _is_complex(default, meta):
        answer = questionary.text(
            f"{name} (JSON, blank = default):",
            instruction=doc,
        ).unsafe_ask()
        if answer is None or answer.strip() == "":
            return None
        try:
            json.loads(answer)
        except json.JSONDecodeError:
            print(f"  (ignoring invalid JSON for '{name}')")
            return None
        return f"{prefix}{name}={answer.strip()}"

    # Everything else -> free text with the default shown.
    default_str = "" if default is None else str(default)
    answer = questionary.text(
        f"{name}:",
        default=default_str,
        instruction=doc,
    ).unsafe_ask()
    if answer is None:
        return None
    answer = answer.strip()
    if answer == default_str:
        return None
    if answer == "" and default is None:
        return None
    return f"{prefix}{name}={answer}"


def _is_complex(default: Any, meta: Dict[str, Any]) -> bool:
    """Heuristic: dict/list-valued fields are edited as JSON, not plain text."""
    if isinstance(default, (dict, list)):
        return True
    type_str = meta.get("type", "")
    return type_str.startswith(("Dict", "dict", "Mapping", "List", "list"))


# ---------------------------------------------------------------------------
# Finalize / session continuation
# ---------------------------------------------------------------------------

def _finalize(questionary, params: ForgeParams, steps: List[str]) -> None:
    """Offer to save the config and/or run it, so the session has an outcome."""
    saved = False
    if questionary.confirm("Save to a YAML file?", default=False).unsafe_ask():
        saved = _save_yaml(questionary, params)

    ran = False
    if questionary.confirm("Run now?", default=False).unsafe_ask():
        _run_now(questionary, params, steps)
        ran = True

    if not saved and not ran:
        print(f"\nAssembled pipeline: {' -> '.join(steps)}")


def _save_yaml(questionary, params: ForgeParams) -> bool:
    """Write the assembled config to a YAML file; return True on success."""
    path = questionary.text("Save path:", default="molforge_config.yaml").unsafe_ask()
    if not path:
        return False

    try:
        import yaml
    except ImportError:
        print("PyYAML is required to save configs (pip install pyyaml).")
        return False

    try:
        with open(path, "w") as f:
            yaml.safe_dump(params.to_dict(), f, sort_keys=False, default_flow_style=False)
    except Exception as exc:
        print(f"Failed to save configuration: {exc}")
        return False

    print(f"Configuration written to {Path(path).resolve()}")
    return True


def _run_now(questionary, params: ForgeParams, steps: List[str]) -> None:
    """Obtain an input, run the pipeline, and report the output directory."""
    if "source" in steps:
        input_data = questionary.text(
            "ChEMBL target ID (e.g. CHEMBL234):"
        ).unsafe_ask()
    else:
        input_data = questionary.text("Input file path (CSV):").unsafe_ask()

    if not input_data or not input_data.strip():
        print("No input provided; skipping run.")
        return
    input_data = input_data.strip()

    from ..forge import MolForge

    try:
        forge = MolForge(params)
        _, success, failed_step = forge.forge(input_data, return_info=True)
    except Exception as exc:
        print(f"Run failed: {exc}")
        return

    output_dir = getattr(forge._pipe, "output_dir", None) or params.output_root
    if success:
        print(f"Run complete. Outputs written to {output_dir}")
    else:
        print(f"Run failed at step '{failed_step}'. Partial outputs in {output_dir}")
