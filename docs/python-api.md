# MolForge Python API

[← Back to README](../README.md) · [Documentation index](../README.md#documentation)

The public API is:

| Method | Description |
|--------|-------------|
| `MolForge(params)` | Construct a pipeline from a `ForgeParams` object, dict, or config file |
| `forge(input_data, ...)` | Forge a single input (ChEMBL ID, CSV path, or DataFrame) into a DataFrame |
| `forge_to_csv(input_data, output_path, ...)` | Forge a single input and write the result to a CSV |
| `batch_forge(input_list, output_dir)` | Forge multiple inputs, writing one CSV per input |

Pipeline behavior and per-actor parameters are configured through `ForgeParams` — see [Pipeline Configuration](configuration.md) and the per-actor pages linked from it.

## Quick Start

```python
from molforge import MolForge, ForgeParams

# Configure pipeline
params = ForgeParams(
    steps=['source', 'chembl', 'curate'],
    source_params={'backend': 'sql'},
    chembl_params={'standard_type': 'IC50', 'standard_units': 'nM'},
    curate_params={'mol_steps': ['desalt', 'neutralize', 'sanitize']}
)

# Run pipeline
forge = MolForge(params)
df = forge.forge("CHEMBL234")
```

## Complete Pipeline Example

```python
from molforge import MolForge, ForgeParams

params = ForgeParams(
    steps=['source', 'chembl', 'curate', 'tokens', 'distributions', 'confs'],

    # Data retrieval
    source_params={'backend': 'sql', 'version': 36},

    # ChEMBL curation
    chembl_params={
        'standard_type': 'IC50',
        'standard_units': 'nM',
        'assay_type': 'B',
        'assay_format': 'protein'
    },

    # Molecule standardization
    curate_params={
        'mol_steps': ['desalt', 'neutralize', 'sanitize', 'handleStereo'],
        'stereo_policy': 'assign',
        'dropna': True
    },

    # Tokenization
    tokens_params={'dynamically_update_vocab': True},

    # Property filtering
    distributions_params={'global_statistical_threshold': 2.0, 'plot_distributions': True},

    # Conformer generation
    confs_params={'backend': 'rdkit', 'max_confs': 200, 'rms_threshold': 0.5}
)

forge = MolForge(params)
df = forge.forge("CHEMBL234")
print(f"Processed {len(df)} molecules")
```

## Loading from a Config File

`MolForge` can be constructed directly from a YAML or JSON configuration file via the `config_path` argument, instead of passing a `ForgeParams` object. Additional keyword arguments override individual config values.

```python
from molforge import MolForge

# From a config file
forge = MolForge(config_path="pipeline_config.yaml")
df = forge.forge("CHEMBL2971")
```

The same config files are used by the [CLI](cli.md) (`molforge run --config ...`). See the [example run](example.md) for a complete config file (`demo/jak2_ic50.yaml`).

## Related

- [Pipeline Configuration](configuration.md)
- [Command-Line Interface](cli.md)
- [Plugin Development](plugins.md)
- [Example run](example.md)
