# Command-Line Interface

[← Back to README](../README.md) · [Documentation index](../README.md#documentation)

MolForge ships a CLI. Invoke it either as the installed `molforge` command or as a module:

```bash
molforge --help
python -m molforge.cli --help    # equivalent
```

Both the CLI and any future UI assemble a `ForgeParams` through a single shared core, `molforge.cli.builder.build_params`, so validation, type coercion, and step routing behave identically everywhere. See [Pipeline Configuration](configuration.md) for the parameters this builds, and the [Python API](python-api.md) for the equivalent programmatic entry point.

## `molforge run`

```bash
molforge run <input> [--config FILE] [--steps a,b,c] [--set step.param=value ...] [--output DIR] [--batch]
```

- `input`: one or more ChEMBL IDs, CSV files, or other data sources.
- `--config, -c`: YAML/JSON configuration file to seed the pipeline.
- `--steps`: comma-separated step list (overrides the config's `steps`).
- `--set`: repeatable parameter override. Use `key=value` for a top-level `ForgeParams` field, or dotted `step.param=value` to target a step. Values are coerced to the target field's type.
- `--output, -o`: output root directory.
- `--batch, -b`: process every input in batch mode.
- `--quiet, -q`: suppress progress messages.

`--steps` and `--set` overrides are applied **on top of** any `--config`.

```bash
# Single ChEMBL target with the standard pipeline
molforge run CHEMBL234

# Override steps and a couple of parameters inline
molforge run CHEMBL234 --steps source,chembl,curate,confs \
    --set confs.backend=openeye --set chembl.standard_type=Ki --set chembl.standard_units=nM

# Seed from a config file, then tweak
molforge run CHEMBL2971 --config demo/jak2_ic50.yaml --set confs.max_confs=50

# Batch processing
molforge run CHEMBL234 CHEMBL279 CHEMBL280 --batch
```

## `molforge init`

```bash
molforge init [--interactive] [--template NAME] [--output FILE]
```

Generates a configuration file. With `--interactive` (or by default in an interactive terminal when no `--template` is given), it launches a `questionary`-based wizard. Otherwise it copies a static template.

Available templates: **api**, **basic**, **conformers**, **distributions**, **full**.

```bash
molforge init                                       # interactive wizard
molforge init --template full --output my_config.yaml
molforge init --template conformers
```

## `molforge info`

```bash
molforge info [actors|backends|<step>|examples|version]
```

- `actors`: list all registered actors (core and plugin).
- `backends`: show backend options (filter with `--actor <step>`).
- `<step>` (e.g. `confs`): describe that step's parameters — name, type, default, choices, and docstring — via live introspection.
- `examples`: print usage examples.
- `version` (default): show the MolForge version.

```bash
molforge info actors
molforge info confs                # describe the confs step's parameters
molforge info backends --actor confs
molforge --diagram                 # print the architecture diagram
```

## `molforge validate`

```bash
molforge validate <config>
```

Loads and validates a configuration file, reporting the resolved pipeline steps.

```bash
molforge validate demo/jak2_ic50.yaml
```

See the [example run](example.md) for a complete config file and the output a `molforge run` produces.
