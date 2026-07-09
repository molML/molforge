<div align="center">
  <img src="logo.png" alt="MolForge Logo" width="400"/>

  # MolForge

  **A configurable pipeline for molecular data processing, curation, and conformer generation**

  [![Version](https://img.shields.io/badge/version-1.0.0-blue)](https://github.com/molML/molforge)
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

</div>

---

## Overview

MolForge is a configurable pipeline for molecular data processing, curation, and conformer generation. The framework processes molecular data through modular steps including ChEMBL data retrieval, molecule standardization, tokenization, property-based filtering, and conformer generation.

---

## Core Features

- **Data Retrieval**: Fetch bioactivity data from ChEMBL (SQL or API backends)
- **Molecule Curation**: Standardize molecules through configurable curation steps
- **SMILES Tokenization**: Convert SMILES strings to token sequences with vocabulary management
- **Property-Based Filtering**: Distribution-based molecular property curation
- **Conformer Generation**: Generate 3D conformer ensembles (RDKit or OpenEye)
- **Plugin System**: Extend functionality with custom actors (ships with scaffold, split, and property demo plugins)
- **Command-Line Interface**: Run, configure, introspect, and validate pipelines from the shell

---

## Installation

### Standard Installation

```bash
git clone https://github.com/molML/molforge.git
cd molforge
conda env create -f environment.yaml
conda activate molforge
```

### With OpenEye Toolkit

OpenEye Toolkit requires a commercial license. Academic licenses are available from [OpenEye Academic Licensing](https://www.eyesopen.com/academic-licensing).

```bash
# Create environment
conda create -n molforge python=3.12
conda activate molforge

# Install OpenEye first
conda install -c openeye openeye-toolkits

# Install remaining dependencies
conda env update -f environment.yaml --prune

# Set license path
export OE_LICENSE=/path/to/oe_license.txt
```

---

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

The public API is:

| Method | Description |
|--------|-------------|
| `MolForge(params)` | Construct a pipeline from a `ForgeParams` object, dict, or config file |
| `forge(input_data, ...)` | Forge a single input (ChEMBL ID, CSV path, or DataFrame) into a DataFrame |
| `forge_to_csv(input_data, output_path, ...)` | Forge a single input and write the result to a CSV |
| `batch_forge(input_list, output_dir)` | Forge multiple inputs, writing one CSV per input |

See the [Python API](docs/python-api.md) and [Command-Line Interface](docs/cli.md) pages for the full interface.

---

## Documentation

- [Pipeline Configuration](docs/configuration.md) — `ForgeParams`/`PipeParams` pipeline-level parameters, base parameters, and how actor params are set.
- **Actor reference** — one page per core actor:
  - [ChEMBLSource (`source`)](docs/actors/source.md) — retrieve bioactivity data from ChEMBL (SQL/API).
  - [ChEMBLCurator (`chembl`)](docs/actors/chembl.md) — type-driven ChEMBL bioactivity curation.
  - [CurateMol (`curate`)](docs/actors/curate.md) — molecule standardization steps.
  - [TokenizeData (`tokens`)](docs/actors/tokens.md) — SMILES tokenization and vocabulary.
  - [CurateDistribution (`distributions`)](docs/actors/distributions.md) — property-distribution and token filtering (includes vocabulary handling).
  - [GenerateConfs (`confs`)](docs/actors/confs.md) — 3D conformer generation (RDKit/OpenEye).
- [Using Actors Standalone](docs/standalone-actors.md) — run a single actor outside the pipeline, and which actors need a pipeline context.
- [Plugin Development](docs/plugins.md) — build your own plugin actor; the included `properties`, `scaffold`, and `split` plugins.
- [Command-Line Interface](docs/cli.md) — `run`, `init`, `info`, `validate`, `--set`, `--steps`, `--diagram`.
- [Python API](docs/python-api.md) — `MolForge`, `forge`/`forge_to_csv`/`batch_forge`, config-file loading.
- [Example Run](docs/example.md) — end-to-end JAK2 (CHEMBL2971) IC50 pipeline with outputs and logs.

---

## Contributing

Contributions are welcome through [GitHub Issues](https://github.com/molML/molforge/issues).

## License

MIT License
