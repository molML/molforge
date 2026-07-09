"""
Main entry point for MolForge CLI.

Provides command-line interface for running pipelines, generating configs,
and viewing system information.
"""

import sys
import argparse
from pathlib import Path
from typing import List, Optional

from . import display
from . import commands


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for CLI."""
    parser = argparse.ArgumentParser(
        prog='molforge',
        description='MolForge - configurable molecular data processing, curation, and conformer generation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  molforge run CHEMBL234
  molforge run CHEMBL234 --steps source,chembl,curate,confs --set confs.backend=openeye
  molforge run compounds.csv --config pipeline.yaml
  molforge init                       # interactive configuration wizard
  molforge info confs                 # show a step's parameters
  molforge validate config.yaml

Documentation: https://github.com/molML/molforge
        """
    )
    
    # Global options
    parser.add_argument(
        '--version',
        action='store_true',
        help='Show version and exit'
    )
    
    parser.add_argument(
        '--diagram',
        action='store_true',
        help='Display architecture diagram and exit'
    )
    
    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Run command
    run_parser = subparsers.add_parser(
        'run',
        help='Execute MolForge pipeline',
        description='Run pipeline on one or more inputs'
    )
    run_parser.add_argument(
        'inputs',
        nargs='+',
        help='Input data (ChEMBL IDs, CSV files, etc.)'
    )
    run_parser.add_argument(
        '--config', '-c',
        help='Configuration file (YAML)'
    )
    run_parser.add_argument(
        '--output', '-o',
        help='Output directory'
    )
    run_parser.add_argument(
        '--steps',
        help='Comma-separated list of steps (overrides config; default: use config or standard pipeline)'
    )
    run_parser.add_argument(
        '--set',
        action='append',
        metavar='KEY=VALUE',
        help=(
            'Override a parameter (repeatable). Use "key=value" for a top-level '
            'ForgeParams field or "step.key=value" to target a step, '
            'e.g. --set confs.max_confs=50 --set chembl.standard_type=Ki'
        )
    )
    run_parser.add_argument(
        '--batch', '-b',
        action='store_true',
        help='Process all inputs in batch mode'
    )
    run_parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress messages'
    )
    
    # Init command
    init_parser = subparsers.add_parser(
        'init',
        help='Generate configuration template',
        description='Create a new configuration file from template'
    )
    
    # Dynamically discover available templates
    template_dir = Path(__file__).parent / 'templates'
    available_templates = sorted([f.stem for f in template_dir.glob('*.yaml')])
    
    init_parser.add_argument(
        '--template', '-t',
        choices=available_templates,
        default=None,
        help=(
            f'Template type (available: {", ".join(available_templates)}). '
            'If omitted, the interactive wizard is launched when possible, '
            'otherwise the "basic" template is used.'
        )
    )
    init_parser.add_argument(
        '--interactive', '-i',
        action='store_true',
        help='Launch the interactive configuration wizard (requires questionary)'
    )
    init_parser.add_argument(
        '--output', '-o',
        default='molforge_config.yaml',
        help='Output file path (default: molforge_config.yaml)'
    )
    
    # Info command
    info_parser = subparsers.add_parser(
        'info',
        help='Display system information',
        description='Show information about actors, backends, or examples'
    )
    info_parser.add_argument(
        'target',
        nargs='?',
        default='version',
        metavar='TARGET',
        help=(
            "What to show: 'actors', 'backends', 'examples', 'version', "
            "or a step name (e.g. 'confs') to describe its parameters and backends "
            "(default: version)"
        )
    )
    info_parser.add_argument(
        '--actor',
        help='Filter backends by actor name'
    )
    
    # Validate command
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate configuration file',
        description='Check configuration file for errors'
    )
    validate_parser.add_argument(
        'config',
        help='Configuration file to validate'
    )
    validate_parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress detailed output'
    )
    
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """
    Main entry point for CLI.
    
    Args:
        argv: Command-line arguments (defaults to sys.argv)
        
    Returns:
        Exit code (0 for success, 1 for failure)
    """
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Handle global flags
    if args.version:
        from .. import __version__
        display.print_simple_banner()
        print(f"Version: {__version__}")
        return 0
    
    if args.diagram:
        display.print_diagram()
        return 0
    
    # No command specified
    if not args.command:
        parser.print_help()
        return 0
    
    # Dispatch to command handlers
    try:
        if args.command == 'run':
            return commands.run_pipeline(args)
        elif args.command == 'init':
            return commands.init_config(args)
        elif args.command == 'info':
            # Convert positional target to flags for show_info
            args.version = (args.target == 'version')
            args.actors = (args.target == 'actors')
            args.backends = (args.target == 'backends')
            args.examples = (args.target == 'examples')
            return commands.show_info(args)
        elif args.command == 'validate':
            return commands.validate_config(args)
        else:
            parser.print_help()
            return 1
    
    except KeyboardInterrupt:
        display.print_warning("\nOperation cancelled by user")
        return 1
    except Exception as e:
        display.print_error(f"Unexpected error: {e}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
