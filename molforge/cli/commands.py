"""
Command implementations for MolForge CLI.

Handles run, init, info, and validate commands.
"""

import sys
from pathlib import Path
from typing import List, Optional
import argparse

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

from ..forge import MolForge
from ..configuration.params import PipeParams
from ..configuration.registry import ActorRegistry
from ..configuration import introspect
from ..backends.registry import BackendRegistry
from . import display
from .builder import build_params


def run_pipeline(args: argparse.Namespace) -> int:
    """
    Execute MolForge pipeline.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    if not args.quiet:
        display.print_simple_banner()

    if args.config and not HAS_YAML and str(args.config).lower().endswith(('.yaml', '.yml')):
        display.print_error("PyYAML is required for YAML config files. Install with: pip install pyyaml")
        return 1

    # Assemble parameters through the shared builder (config + steps + --set),
    # so overrides are applied on top of any --config and every selected step's
    # params get initialized and validated.
    try:
        overrides = getattr(args, 'set', None)
        params = build_params(
            steps=args.steps,
            config_path=args.config,
            overrides=overrides,
            output_root=args.output,
        )
    except Exception as e:
        display.print_error(f"Configuration error: {e}")
        return 1

    # Create forge instance (output_root already baked into params)
    try:
        if not args.quiet:
            display.print_info("Initializing pipeline...")
        forge = MolForge(params)
        if not args.quiet:
            display.print_success("Configuration validated")
    except Exception as e:
        display.print_error(f"Initialization failed: {e}")
        return 1

    # Process inputs (no output_dir passed - already set at construction)
    if args.batch:
        return run_batch(forge, args.inputs, None, args.quiet)
    else:
        return run_single(forge, args.inputs[0], None, args.quiet)


def run_single(forge: MolForge, input_data: str, output_dir: Optional[str] = None, quiet: bool = False) -> int:
    """
    Run single pipeline execution.

    Args:
        forge: MolForge instance
        input_data: Input data (ChEMBL ID, CSV file, etc.)
        output_dir: Ignored (output_root set at MolForge construction)
        quiet: Suppress progress messages

    Returns:
        Exit code
    """
    if not quiet:
        print(f"\nProcessing: {input_data}")

    try:
        # Execute pipeline
        df, success, failed_step = forge.forge(
            input_data,
            return_info=True
        )
        
        if success:
            if not quiet:
                display.print_success("Pipeline completed successfully")
                print(f"Output: {forge._pipe.output_dir}")
                print(f"Rows: {len(df) if df is not None else 0}")
            return 0
        else:
            display.print_error(f"Pipeline failed at step: {failed_step}")
            return 1
            
    except Exception as e:
        display.print_error(f"Execution failed: {e}")
        if not quiet:
            import traceback
            traceback.print_exc()
        return 1


def run_batch(forge: MolForge, inputs: List[str], output_dir: Optional[str] = None, quiet: bool = False) -> int:
    """
    Run batch pipeline executions.

    Args:
        forge: MolForge instance
        inputs: List of inputs to process
        output_dir: Ignored (output_root set at MolForge construction)
        quiet: Suppress progress messages

    Returns:
        Exit code (1 if any failed, 0 if all succeeded)
    """
    if not quiet:
        print(f"\nBatch processing {len(inputs)} inputs...")

    results = {}
    failed_count = 0

    for i, input_data in enumerate(inputs, 1):
        if not quiet:
            print(f"\n[{i}/{len(inputs)}] Processing: {input_data}")

        try:
            df, success, failed_step = forge.forge(
                input_data,
                return_info=True
            )
            
            results[input_data] = success
            if success:
                if not quiet:
                    display.print_success(f"Completed: {input_data}")
            else:
                failed_count += 1
                display.print_error(f"Failed at {failed_step}: {input_data}")
                
        except Exception as e:
            failed_count += 1
            results[input_data] = False
            display.print_error(f"Error processing {input_data}: {e}")
    
    # Summary
    if not quiet:
        print("\n" + "="*60)
        print(f"Batch Summary: {len(inputs) - failed_count}/{len(inputs)} succeeded")
        if failed_count > 0:
            print(f"Failed inputs:")
            for inp, success in results.items():
                if not success:
                    print(f"  - {inp}")
    
    return 1 if failed_count > 0 else 0


def init_config(args: argparse.Namespace) -> int:
    """
    Generate a configuration.

    Launches the interactive wizard when requested with ``--interactive`` (or
    by default when no ``--template`` is given and an interactive terminal with
    questionary is available). Otherwise copies a static template.

    Args:
        args: Parsed command-line arguments

    Returns:
        Exit code
    """
    template = getattr(args, 'template', None)
    interactive = getattr(args, 'interactive', False)

    # Decide whether to run the wizard.
    if interactive or (template is None and _interactive_capable()):
        result = _launch_wizard(args)
        if result is not None:
            return result
        # questionary unavailable and not explicitly requested -> fall through

    if not HAS_YAML:
        display.print_error("PyYAML is required for config generation. Install with: pip install pyyaml")
        return 1

    template_name = template or 'basic'
    output_path = Path(args.output)
    
    # Get template content
    try:
        template_content = get_template(template_name)
    except ValueError as e:
        display.print_error(str(e))
        return 1
    
    # Write to file
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(template_content)
        
        display.print_success(f"Configuration template created: {output_path}")
        display.print_info(f"Template type: {template_name}")
        display.print_info(f"Edit the file and run: molforge run <input> --config {output_path}")
        return 0
        
    except Exception as e:
        display.print_error(f"Failed to write config file: {e}")
        return 1


def _interactive_capable() -> bool:
    """Return True if the wizard can be launched (questionary + TTY)."""
    from . import wizard
    return wizard.questionary_available() and bool(getattr(sys.stdout, 'isatty', lambda: False)())


def _launch_wizard(args: argparse.Namespace) -> Optional[int]:
    """
    Run the interactive configuration wizard.

    Returns an exit code, or ``None`` to signal the caller to fall back to the
    static template flow (only when questionary is missing and the wizard was
    not explicitly requested).
    """
    from . import wizard

    if not wizard.questionary_available():
        if getattr(args, 'interactive', False):
            display.print_error(
                "Interactive mode requires 'questionary'. Install with: pip install questionary"
            )
            return 1
        return None  # fall back to template copy

    try:
        params = wizard.run_wizard()
    except Exception as e:
        display.print_error(f"Wizard failed: {e}")
        return 1

    if params is None:
        display.print_warning("Configuration wizard cancelled.")
        return 1

    display.print_success("Configuration ready.")
    display.print_info(f"Pipeline steps: {' -> '.join(params.steps)}")
    return 0


def show_info(args: argparse.Namespace) -> int:
    """
    Display system information.

    Args:
        args: Parsed command-line arguments

    Returns:
        Exit code
    """
    if getattr(args, 'actors', False):
        print_actors_info()
        return 0

    if getattr(args, 'backends', False):
        print_backends_info(getattr(args, 'actor', None))
        return 0

    if getattr(args, 'examples', False):
        print_examples()
        return 0

    if getattr(args, 'version', False):
        print_version()
        return 0

    # A bare step name (e.g. `molforge info confs`) describes that step.
    target = getattr(args, 'target', None)
    if target and target not in ('version', 'actors', 'backends', 'examples'):
        return print_step_info(target)

    # Default: show version
    print_version()
    return 0


def print_step_info(step_name: str) -> int:
    """Print a step's parameter fields (name/type/default/choices/doc) and backends."""
    try:
        info = introspect.describe_step(step_name)
    except ValueError as e:
        display.print_error(str(e))
        return 1

    kind = 'plugin' if info['is_plugin'] else 'core'
    display.print_section(f"Step: {step_name} ({kind})")

    if info['backends']:
        print(f"  Backends: {', '.join(info['backends'])}")

    fields_meta = info['fields']
    print(f"\n  Parameters ({len(fields_meta)}):")
    for meta in fields_meta:
        header = f"    {meta['name']} : {meta['type']}"
        if meta['choices'] is not None:
            header += f"  choices={meta['choices']}"
        print(header)
        print(f"      default: {meta['default']!r}")
        if meta['doc']:
            print(f"      {meta['doc']}")
    return 0


def validate_config(args: argparse.Namespace) -> int:
    """
    Validate configuration file.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Exit code
    """
    if not HAS_YAML:
        display.print_error("PyYAML is required. Install with: pip install pyyaml")
        return 1
    
    config_path = Path(args.config)
    
    if not config_path.exists():
        display.print_error(f"Config file not found: {config_path}")
        return 1
    
    try:
        params = load_config(str(config_path))
        display.print_success(f"Configuration is valid: {config_path}")
        
        # Show summary
        if not args.quiet:
            print(f"\nPipeline steps ({len(params.steps)}):")
            for i, step in enumerate(params.steps, 1):
                print(f"  {i}. {step}")
        
        return 0
        
    except Exception as e:
        display.print_error(f"Invalid configuration: {e}")
        return 1


# Helper functions

def load_config(config_path: str) -> PipeParams:
    """
    Load a YAML/JSON configuration file into a validated ForgeParams.

    Delegates to ``build_params`` so that step blocks are mapped to their
    parameter classes dynamically via the registry (core *and* plugin steps),
    rather than through a hardcoded actor->param table.
    """
    return build_params(config_path=config_path)


def get_template(template_name: str) -> str:
    """Get configuration template by name."""
    # Load templates from yaml files
    template_dir = Path(__file__).parent / 'templates'
    template_file = template_dir / f"{template_name}.yaml"
    
    if not template_file.exists():
        available_templates = [f.stem for f in template_dir.glob('*.yaml')]
        available = ', '.join(sorted(available_templates))
        raise ValueError(f"Unknown template '{template_name}'. Available: {available}")
    
    with open(template_file, 'r') as f:
        return f.read()


def print_version():
    """Print version information."""
    from .. import __version__
    display.print_simple_banner()
    print(f"Version: {__version__}")
    print("Python API: Available via 'from molforge import MolForge'")
    print("Documentation: https://github.com/molML/molforge")


def print_actors_info():
    """Print information about available actors."""
    display.print_section("Available Actors")
    
    for step in ActorRegistry.get_all_step_names():
        info = ActorRegistry.get_actor_info(step)
        if info:
            class_name = info['class'].__name__
            is_plugin = info.get('is_plugin', False)
            plugin_marker = " [plugin]" if is_plugin else ""
            
            display.print_table_row(step, f"{class_name}{plugin_marker}", width1=20)
    
    print("\nRun 'molforge info backends' to see backend options")


def print_backends_info(actor_filter: Optional[str] = None):
    """Print information about available backends."""
    display.print_section("Available Backends")
    
    backends = BackendRegistry.list_backends(actor_filter)
    
    if not backends:
        print("  No backends registered")
        return
    
    for actor_type, backend_list in backends.items():
        print(f"\n  {actor_type}:")
        for backend in backend_list:
            print(f"    • {backend}")


def print_examples():
    """Print usage examples."""
    display.print_section("Usage Examples")
    
    examples = [
        ("Run pipeline with ChEMBL target", "molforge run CHEMBL234"),
        ("Run with config file", "molforge run CHEMBL234 --config pipeline.yaml"),
        ("Run from CSV file", "molforge run compounds.csv --output ./results"),
        ("Batch processing", "molforge run CHEMBL234 CHEMBL279 --batch"),
        ("Generate config template", "molforge init --template basic --output config.yaml"),
        ("Show available actors", "molforge info actors"),
        ("Show backends", "molforge info backends"),
        ("Validate config", "molforge validate config.yaml"),
    ]
    
    for desc, cmd in examples:
        print(f"\n  {desc}:")
        print(f"  $ {cmd}")
    
    print()

