# MolForge CLI

Command-line interface for MolForge.

## Installation

Install MolForge with CLI support:

```bash
pip install -e ".[cli]"
```

This installs the `molforge` command with optional PyYAML support for configuration files. The interactive `init` wizard additionally requires `questionary`.

## Quick Start

```bash
# Show help
molforge --help
python -m molforge.cli --help      # equivalent

# Display architecture diagram
molforge --diagram

# Run pipeline with a ChEMBL target
molforge run CHEMBL234

# Generate a configuration template
molforge init --template basic --output config.yaml

# Run with a configuration file
molforge run CHEMBL234 --config config.yaml

# Batch processing
molforge run CHEMBL234 CHEMBL279 CHEMBL280 --batch
```

## Shared Parameter Core

The CLI (and any future UI) assembles a `ForgeParams` object through a single shared function, `molforge.cli.builder.build_params`. It merges, in order, a loaded config file, an explicit `--steps` list, `--set` overrides, and an output root — so validation, type coercion, and step routing behave identically across every entry point.

## Commands

### `molforge run`

Execute the MolForge pipeline on one or more inputs.

**Usage:**
```bash
molforge run <input> [options]
```

**Arguments:**
- `inputs`: one or more ChEMBL IDs, CSV files, or other data sources.

**Options:**
- `--config, -c`: configuration file (YAML/JSON).
- `--output, -o`: output directory.
- `--steps`: comma-separated list of steps (overrides the config's `steps`).
- `--set KEY=VALUE`: repeatable parameter override. Use `key=value` for a top-level `ForgeParams` field, or dotted `step.key=value` to target a step (e.g. `--set confs.max_confs=50 --set chembl.standard_type=Ki`). Values are coerced to the target field's type.
- `--batch, -b`: process all inputs in batch mode.
- `--quiet, -q`: suppress progress messages.

`--steps` and `--set` overrides are applied on top of any `--config`.

**Examples:**
```bash
# Single ChEMBL target
molforge run CHEMBL234

# Override steps and parameters inline
molforge run CHEMBL234 --steps source,chembl,curate,confs --set confs.backend=openeye

# From CSV file
molforge run compounds.csv --output ./results

# With configuration, then tweak
molforge run CHEMBL2971 --config pipeline.yaml --set distributions.global_statistical_threshold=3.0

# Batch processing
molforge run CHEMBL234 CHEMBL279 --batch
```

### `molforge init`

Generate a configuration file. With `--interactive` (or by default in an interactive terminal when no `--template` is given and `questionary` is installed), launches an interactive wizard. Otherwise copies a static template.

**Usage:**
```bash
molforge init [options]
```

**Options:**
- `--interactive, -i`: launch the interactive configuration wizard (requires `questionary`).
- `--template, -t`: template name (see below).
- `--output, -o`: output file path (default: `molforge_config.yaml`).

**Templates:**
- **api**: minimal pipeline using the ChEMBL API backend.
- **basic**: source → chembl → curate (minimal pipeline).
- **conformers**: conformer generation pipeline.
- **distributions**: pipeline with distribution-based property filtering.
- **full**: all core pipeline steps.

**Examples:**
```bash
# Interactive wizard
molforge init --interactive

# Basic template
molforge init --template basic

# Full pipeline template
molforge init --template full --output my_config.yaml
```

### `molforge info`

Display system information.

**Usage:**
```bash
molforge info [target]
```

**Targets:**
- `version`: show MolForge version (default).
- `actors`: list available actors (core and plugin).
- `backends`: show backend options.
- `examples`: display usage examples.
- `<step>`: a bare step name (e.g. `confs`) describes that step's parameters — name, type, default, choices, and docstring — via live introspection.

**Options:**
- `--actor`: filter backends by actor name (with the `backends` target).

**Examples:**
```bash
# Show version
molforge info
molforge info version

# List actors
molforge info actors

# Describe a step's parameters
molforge info confs

# Show backends
molforge info backends
molforge info backends --actor confs

# Usage examples
molforge info examples
```

### `molforge validate`

Validate a configuration file.

**Usage:**
```bash
molforge validate <config>
```

**Arguments:**
- `config`: configuration file to validate.

**Options:**
- `--quiet, -q`: suppress detailed output.

**Examples:**
```bash
molforge validate config.yaml
molforge validate config.yaml --quiet
```

## Configuration Files

Configuration files use YAML (or JSON) to define pipeline parameters. Each step is a nested block keyed by step name; an `output:` block maps `dir` → `output_root` and `checkpoints` → `write_checkpoints`.

### Basic Configuration

```yaml
# Pipeline steps
steps:
  - source
  - chembl
  - curate

# Actor configurations (keyed by step name)
source:
  backend: sql
  version: 36

chembl:
  standard_type: IC50
  standard_units: nM
  assay_type: B
  assay_format: protein

curate:
  mol_steps: [desalt, neutralize, sanitize]
  stereo_policy: assign

# Output settings
output:
  dir: ./molforge_output
  checkpoints: false
```

### Full Configuration

Run `molforge init --template full` for a complete example with all actors.

## Global Options

Available for all commands:

- `--version`: show version and exit.
- `--diagram`: display the architecture diagram and exit.
- `--help, -h`: show the help message.

## Python API

The CLI complements the Python API. Both interfaces are available:

```python
from molforge import MolForge, ForgeParams

params = ForgeParams(steps=['source', 'chembl', 'curate'])
forge = MolForge(params)
df = forge.forge('CHEMBL234')
```

```bash
# CLI equivalent
molforge run CHEMBL234 --steps source,chembl,curate
```

## Dependencies

- **Core**: argparse (built-in)
- **Optional**: `pyyaml` (config files), `questionary` (interactive `init` wizard)

Install with CLI support:
```bash
pip install -e ".[cli]"
```

## Architecture

```
molforge/cli/
├── __init__.py       # Package initialization
├── main.py           # Entry point & argument parsing
├── commands.py       # Command implementations (run, init, info, validate)
├── builder.py        # Shared build_params core (config + steps + --set)
├── wizard.py         # Interactive configuration wizard
├── display.py        # Display utilities & formatting
└── templates/        # Configuration templates
    ├── api.yaml
    ├── basic.yaml
    ├── conformers.yaml
    ├── distributions.yaml
    └── full.yaml
```

## Troubleshooting

**Command not found:**
```bash
pip install -e ".[cli]"
```

**YAML errors:**
```bash
pip install pyyaml
```

**Interactive wizard unavailable:**
```bash
pip install questionary
```

## See Also

- [Main README](../../README.md): MolForge overview and full parameter reference.
