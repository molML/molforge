<div align="center">
  <img src="logo.png" alt="MolForge Logo" width="400"/>

  # MolForge

  **A configurable pipeline for molecular data processing, curation, and conformer generation**

  [![Status](https://img.shields.io/badge/status-public%20beta-orange)](https://github.com/molML/molforge)
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

</div>

---

## Overview

MolForge processes molecular data through modular, configurable steps including ChEMBL data retrieval, molecule standardization, tokenization, and conformer generation. It's designed for cheminformatics workflows and machine learning applications.

**This is a beta release.** We welcome bug reports, feature requests, and feedback through [GitHub Issues](https://github.com/molML/molforge/issues).

## Installation

### Option 1: Using Conda (Recommended)

The recommended way to install MolForge is using conda, which handles all dependencies including RDKit:

```bash
# Clone the repository
git clone https://github.com/molML/molforge.git
cd molforge

# Create and activate the conda environment
conda env create -f environment.yaml
conda activate molforge
```

**Verify installation:**
```bash
python -c "import molforge; print(f'MolForge v{molforge.__version__} installed successfully')"
```

### Option 2: Using pip

If you prefer pip, make sure you have RDKit installed first:

```bash
# Install RDKit (if not already installed)
conda install -c conda-forge rdkit

# Clone and install MolForge
git clone https://github.com/molML/molforge.git
cd molforge
pip install -e .

# Optional: Install CLI dependencies
pip install -e ".[cli]"
```

**Verify installation:**
```bash
python -c "import molforge; print(f'MolForge v{molforge.__version__} installed successfully')"
```

### Optional Dependencies

- **OpenEye Toolkit** (for OpenEye conformer backend, requires commercial license):
  ```bash
  conda install -c openeye openeye-toolkits
  ```

- **MolBox** (for MolBox integration):
  ```bash
  pip install molbox
  ```

**Requirements:** Python ≥3.8, pandas, rdkit, numpy

## Quick Start

```python
from molforge import MolForge, ForgeParams, CurateMolParams

# Configure pipeline
params = ForgeParams(
    steps=['source', 'chembl', 'curate', 'tokens'],
    curate_params=CurateMolParams(
        mol_steps=['desalt', 'neutralize', 'sanitize']
    )
)

# Run pipeline
forge = MolForge(params)
df = forge.forge("CHEMBL234")  # Process ChEMBL target
```

## Pipeline Steps

**Core Actors:**
- `source` - Fetch bioactivity data from ChEMBL
- `chembl` - Curate ChEMBL bioactivity data
- `curate` - Standardize molecules (desalting, neutralization, stereochemistry)
- `tokens` - Tokenize SMILES strings
- `distributions` - Distribution-based molecular filtering
- `confs` - Generate 3D conformer ensembles

## Example: Molecule Curation

```python
from molforge import ForgeParams, CurateMolParams

params = ForgeParams(
    steps=['curate'],
    curate_params=CurateMolParams(
        SMILES_column='SMILES',
        mol_steps=['desalt', 'neutralize', 'sanitize', 'handleStereo'],
        stereo_policy='assign',
        dropna=True
    )
)

forge = MolForge(params)
df = forge.forge("molecules.csv")
```

## Example: Conformer Generation

```python
from molforge import ForgeParams, GenerateConfsParams

params = ForgeParams(
    steps=['curate', 'confs'],
    confs_params=GenerateConfsParams(
        backend='rdkit',
        max_confs=200,
        rms_threshold=0.5
    )
)

forge = MolForge(params)
df = forge.forge("molecules.csv")
```

## Command Line Interface

```bash
# Run pipeline from config file
molforge run --config config.yaml --input CHEMBL234

# Get information about available actors
molforge info actors

# Show architecture diagram
molforge info architecture
```

## Contributing & Feedback

This is a **public beta release**. We actively welcome:

- 🐛 **Bug reports** - Found something broken? Let us know!
- 💡 **Feature requests** - Have ideas for improvements?
- 📝 **Feedback** - How can we make MolForge better?
- 🔧 **Pull requests** - Contributions are welcome!

Please open an issue on our [GitHub Issues page](https://github.com/molML/molforge/issues).

## Documentation

For detailed documentation, configuration examples, and API reference, visit our [documentation](https://github.com/molML/molforge/wiki) (coming soon).

## License

MIT License
