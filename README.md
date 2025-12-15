<div align="center">
  <img src="logo.png" alt="MolForge Logo" width="400"/>

  # MolForge

  **A configurable pipeline for molecular data processing, curation, and conformer generation**

  [![Status](https://img.shields.io/badge/status-public%20beta-orange)](https://github.com/molML/molforge)
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

</div>

---

## Overview

MolForge is a configurable pipeline for molecular data processing, curation, and conformer generation. The framework processes molecular data through modular steps including ChEMBL data retrieval, molecule standardization, tokenization, and conformer generation.

---

## Core Features

- **Data Retrieval**: Fetch bioactivity data from ChEMBL (SQL or API backends)
- **Molecule Curation**: Standardize molecules through configurable curation steps
- **SMILES Tokenization**: Convert SMILES strings to token sequences with vocabulary management
- **Property-Based Filtering**: Distribution-based molecular property curation
- **Conformer Generation**: Generate 3D conformer ensembles (RDKit or OpenEye)
- **Plugin System**: Extend functionality with custom actors

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

OpenEye Toolkit requires a commercial license. Academic licenses available from [OpenEye Academic Licensing](https://www.eyesopen.com/academic-licensing).

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
    chembl_params={'standard_type': 'IC50'},
    curate_params={'mol_steps': ['desalt', 'neutralize', 'sanitize']}
)

# Run pipeline
forge = MolForge(params)
df = forge.forge("CHEMBL234")
```

---

## Base Parameters

All actor parameter classes inherit from `BaseParams`, which provides common functionality.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `verbose` | `True` | `bool` | Enable verbose logging output |

---

## Pipeline Configuration

The `ForgeParams` class configures the overall pipeline behavior and individual actor parameters.

### Pipeline Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `steps` | `['source', 'chembl', 'curate', 'tokens', 'distributions']` | `List[str]` | List of actor steps to execute in order |
| `return_none_on_fail` | `False` | `bool` | Return None instead of raising exception on failure |
| `output_root` | `"MOLFORGE_OUT"` | `str` | Root directory for output files |
| `write_output` | `True` | `bool` | Write pipeline output to disk |
| `write_checkpoints` | `False` | `bool` | Save intermediate results at each step |
| `console_log_level` | `'INFO'` | `str` | Console logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR') |
| `file_log_level` | `'DEBUG'` | `str` | File logging level |
| `override_actor_params` | `True` | `bool` | Use pipeline-level params to override actor params |

### Actor Parameters

Each actor has a corresponding parameter class. Specify actor parameters using the `{actor}_params` attribute:

```python
ForgeParams(
    steps=['source', 'chembl'],
    source_params=ChEMBLSourceParams(backend='sql'),
    chembl_params={'standard_type': 'IC50'}  # Can also use dict
)
```

---

## Actor Reference

### 1. ChEMBLSource (`source`)

Retrieves bioactivity data from ChEMBL database. Supports SQL (local database) and API (web service) backends. SQL backend is significantly faster than API once downloaded (<1s vs. >30s per protein, depending on the target).

#### Shared Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `backend` | `'sql'` | `'sql'` \| `'api'` | Data source backend |
| `n` | `1000` | `int` | Number of activity entries to retrieve per query (> 0) |
| `search_all` | `True` | `bool` | Fetch all entries instead of limiting to n |

#### SQL Backend Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `db_path` | `None` | `str` \| `None` | Path to local ChEMBL SQLite database |
| `version` | `'latest'` | `str` \| `int` | ChEMBL version (e.g., 36) or 'latest' |
| `auto_download` | `True` | `bool` | Automatically download database if not found |
| `download_dir` | `"./data/chembl"` | `str` | Directory to store downloaded database |

#### API Backend Parameters

API backend currently uses only shared parameters.

#### Example

```python
# SQL backend
source_params = ChEMBLSourceParams(
    backend='sql',
    version=36,
    auto_download=True
)

# API backend
source_params = ChEMBLSourceParams(
    backend='api',
    n=5000,
    search_all=False
)
```

---

### 2. ChEMBLCurator (`chembl`)

Curates ChEMBL bioactivity data with automatic type-driven behavior. Curation logic adapts based on `standard_type` category (potency, kinetic, ADMET, activity).

#### Core Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `standard_type` | `'IC50'` | `str` | Measurement type (see supported types below) |
| `standard_units` | `None` | `str` \| `None` | Measurement units (auto-set from type if None) |
| `standard_relation` | `'='` | `'='` \| `'>'` \| `'<'` \| `'>='` \| `'<='` \| `'~'` \| `None` | Relation operator for measurements |
| `assay_type` | `'B'` | `'B'` \| `'F'` \| `'A'` \| `'T'` | Assay type (Binding, Functional, ADMET, Toxicity) |
| `target_organism` | `'Homo sapiens'` | `str` | Target organism for filtering |
| `assay_format` | `'protein'` | `str` | Assay system format (see options below), maps to BAO format and is mostly for readability |

#### Quality Control Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `top_n` | `5` | `int` | Number of top-ranked assays to consider per target |
| `mistakes_only` | `True` | `bool` | Only flag suspicious unit conversion errors |
| `error_margin` | `0.0001` | `float` | Tolerance for suspicious pair detection |
| `std_threshold` | `0.5` | `float` | Standard deviation threshold for filtering |
| `range_threshold` | `0.5` | `float` | Range threshold for filtering |

#### Advanced Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `mutant_regex` | `'mutant\|mutation\|variant'` | `str` | Regex pattern to identify mutant assays |
| `allosteric_regex` | `'allosteric'` | `str` | Regex pattern to identify allosteric assays |
| `C` | `1e9` | `float` | Conversion constant for molarity (auto-set from units). Only change if you desire non-standard unit conversion |
| `SMILES_column` | `'canonical_smiles'` | `str` | Column name for SMILES strings |
| `bao_format` | Auto-set | `str` | BioAssay Ontology format code (derived from assay_format) |

#### Supported Standard Types

**Potency** (IC50, EC50, AC50, Ki, Kd, XC50)
**Kinetic** (kon, k_off, koff, Kon, Koff)
**ADMET** (Log BB, Log D, Log P, Clearance, T1/2, Solubility, Permeability, Caco-2, MDCK)
**Activity** (Inhibition, Percent Effect, Activity, Potency)

Use `ChEMBLCuratorParams.list_supported_types()` to see all supported types.

#### Assay Format Options

- `'general'` - Unspecified format
- `'protein'` - Single protein format (biochemical)
- `'cell'` - Cell-based format
- `'organism'` - Organism-based format (in vivo)
- `'tissue'` - Tissue-based format (ex vivo)
- `'microsome'` - Microsome format (metabolic stability)

#### Example

```python
# IC50 potency assay
chembl_params = ChEMBLCuratorParams(
    standard_type='IC50',
    standard_units='nM',
    assay_type='B',
    assay_format='protein'
)

# ADMET assay
chembl_params = ChEMBLCuratorParams(
    standard_type='Log BB',
    assay_type='A',
    assay_format='organism',
    target_organism='Mus musculus'
)
```

---

### 3. CurateMol (`curate`)

Standardizes molecules through configurable curation steps including desalting, neutralization, stereochemistry handling, and SMILES canonicalization.

#### Core Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `mol_steps` | See below | `List[str]` | Ordered list of molecule curation steps |
| `smiles_steps` | `['canonical', 'kekulize']` | `List[str]` | SMILES processing steps |
| `SMILES_column` | `'canonical_smiles'` | `str` | Column name for SMILES strings |
| `dropna` | `True` | `bool` | Remove molecules with missing SMILES |
| `check_duplicates` | `True` | `bool` | Check for duplicate molecules |
| `duplicates_policy` | `'first'` | `'first'` \| `'last'` \| `False` | How to handle duplicates |

Default `mol_steps`: `['desalt', 'removeIsotope', 'removeHs', 'tautomers', 'neutralize', 'sanitize', 'handleStereo', 'computeProps']`

#### Step-Specific Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `neutralize_policy` | `'keep'` | `'keep'` \| `'remove'` | Keep or remove molecules after neutralization failure |
| `desalt_policy` | `'keep'` | `'keep'` \| `'remove'` | Keep or remove molecules after desalting |
| `brute_force_desalt` | `False` | `bool` | Use aggressive SMILES-based desalting for problematic molecules |
| `max_tautomers` | `512` | `int` | Maximum tautomers to enumerate |
| `step_timeout` | `60` | `int` | Timeout in seconds for individual curation steps |

#### Stereochemistry Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `stereo_policy` | `'assign'` | `'keep'` \| `'remove'` \| `'assign'` \| `'enumerate'` | Stereochemistry handling strategy |
| `assign_policy` | `'random'` | `'first'` \| `'random'` \| `'lowest'` | Selection strategy for assignment |
| `max_isomers` | `32` | `int` | Maximum stereoisomers to enumerate |
| `try_embedding` | `False` | `bool` | Use 3D embedding for stereocenter assignment |
| `only_unassigned` | `True` | `bool` | Only process unassigned stereocenters |
| `only_unique` | `True` | `bool` | Remove duplicate stereoisomers |
| `random_seed` | `42` | `int` | Random seed for reproducible assignment |

#### Available Curation Steps

**Molecule Steps** (`mol_steps`):
- `'desalt'` - Remove counterions and salts
- `'removeIsotope'` - Remove isotope labels
- `'removeHs'` - Remove explicit hydrogens
- `'tautomers'` - Canonicalize tautomeric form
- `'neutralize'` - Neutralize charged species
- `'sanitize'` - Apply RDKit sanitization
- `'handleStereo'` - Process stereochemistry
- `'computeProps'` - Compute molecular properties

**SMILES Steps** (`smiles_steps`):
- `'canonical'` - Generate canonical SMILES
- `'kekulize'` - Convert to Kekulé form

#### Example

```python
# Standard curation
curate_params = CurateMolParams(
    mol_steps=['desalt', 'neutralize', 'sanitize'],
    stereo_policy='assign',
    dropna=True
)

# Enumerate stereoisomers
curate_params = CurateMolParams(
    mol_steps=['desalt', 'neutralize', 'sanitize', 'handleStereo'],
    stereo_policy='enumerate',
    max_isomers=16
)
```

---

### 4. TokenizeData (`tokens`)

Tokenizes SMILES strings into token sequences for machine learning applications. Supports dynamic vocabulary updating which may be useful for building a single vocabulary across several forge inputs (e.g., an array of protein targets).

#### Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `vocab_file` | `None` | `str` \| `None` | Path to existing vocabulary file |
| `SMILES_column` | `'curated_smiles'` | `str` | Column name for SMILES strings |
| `dynamically_update_vocab` | `True` | `bool` | Update vocabulary with new tokens from input |

#### Example

```python
# Dynamic vocabulary
tokens_params = TokenizeDataParams(
    SMILES_column='curated_smiles',
    dynamically_update_vocab=True
)

# Fixed vocabulary
tokens_params = TokenizeDataParams(
    vocab_file='vocab.json',
    dynamically_update_vocab=False
)
```

---

### 5. CurateDistribution (`distributions`)

Filters molecules based on molecular property distributions using statistical or quantile thresholds. A global threshold can be set for all specified properties for easy-of-use. Beware of the fact that distribution thresholds apply to set set of molecules at that stage of the pipeline, these are typically non-normal, low sample size (±1k mols). For large scale curation, it is best to calibrate thresholds on a larger chemical space.

#### Core Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `SMILES_column` | `'curated_smiles'` | `str` | Column name for SMILES strings |
| `tokens_column` | `'tokens'` | `str` | Column name for token sequences |
| `properties` | `None` | `List[str]` \| `'all'` \| `None` | Properties to compute (None = auto-detect) |
| `compute_properties` | `True` | `bool` | Compute molecular properties |
| `dropna` | `True` | `bool` | Remove molecules with missing values |

#### Threshold Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `thresholds` | `{}` | `Dict[str, PropertyThreshold]` | Property-specific threshold configurations |
| `global_statistical_threshold` | `None` | `float` \| `None` | Global statistical threshold (e.g., 2.0 for ±2σ) |
| `global_quantile_threshold` | `None` | `float` \| `None` | Global quantile threshold (e.g., 0.025 for 2.5%-97.5%) |

#### Token Curation Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `curate_tokens` | `True` | `bool` | Enable token-based filtering |
| `filter_unknown_tokens` | `True` | `bool` | Remove molecules with unknown tokens |
| `token_frequency_threshold` | `None` | `float` \| `None` | Minimum token frequency percentage (0-100) |

#### Visualization Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `plot_distributions` | `True` | `bool` | Generate distribution plots |
| `perform_pca` | `False` | `bool` | Perform PCA analysis on properties |

#### Available Properties

`num_atoms`, `num_rings`, `size_largest_ring`, `num_tokens`, `tokens_atom_ratio`, `c_atom_ratio`, `longest_aliph_c_chain`, `molecular_weight`, `logp`, `tpsa`, `num_rotatable_bonds`, `num_h_donors`, `num_h_acceptors`, `fsp3`, `num_aromatic_rings`, `num_stereocenters`, `num_heteroatoms`, `heteroatom_ratio`

#### PropertyThreshold Configuration

`PropertyThreshold` supports three types of thresholds (priority: absolute > statistical > quantile):

| Parameter | Type | Description |
|-----------|------|-------------|
| `min_value` | `float` \| `None` | Absolute minimum value |
| `max_value` | `float` \| `None` | Absolute maximum value |
| `statistical_lower` | `float` \| `None` | Lower bound as mean - n×std (e.g., -2.0) |
| `statistical_upper` | `float` \| `None` | Upper bound as mean + n×std (e.g., 2.0) |
| `quantile_lower` | `float` \| `None` | Lower percentile (0-1, e.g., 0.025) |
| `quantile_upper` | `float` \| `None` | Upper percentile (0-1, e.g., 0.975) |

#### Example

```python
from molforge import PropertyThreshold

# Global statistical threshold
distributions_params = CurateDistributionParams(
    global_statistical_threshold=2.0,  # ±2σ
    plot_distributions=True
)

# Property-specific thresholds
distributions_params = CurateDistributionParams(
    thresholds={
        'molecular_weight': PropertyThreshold(min_value=200, max_value=500),
        'logp': PropertyThreshold(statistical_lower=-2.0, statistical_upper=2.0),
        'num_atoms': PropertyThreshold(quantile_lower=0.05, quantile_upper=0.95)
    }
)
```

---

### 6. GenerateConfs (`confs`)

Generates 3D conformer ensembles for molecules. Supports RDKit and OpenEye backends with backend-specific optimization parameters.

#### Shared Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `backend` | `'rdkit'` | `'rdkit'` \| `'openeye'` | Conformer generation backend |
| `max_confs` | `200` | `int` | Maximum conformers to generate per molecule (> 0) |
| `rms_threshold` | `0.5` | `float` | RMS threshold for conformer pruning in Angstroms (> 0) |
| `SMILES_column` | `'curated_smiles'` | `str` | Column name for SMILES strings |
| `names_column` | `'molecule_chembl_id'` | `str` | Column name for molecule identifiers |
| `dropna` | `True` | `bool` | do not attempt generation for invalid SMILES or molecules |

#### RDKit Backend Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `use_random_coords` | `True` | `bool` | Use random coordinates for initial embedding |
| `random_seed` | `42` | `int` | Random seed for reproducibility (≥ 0) |
| `num_threads` | `0` | `int` | Number of threads (0 = all available cores) |
| `use_uff` | `True` | `bool` | Use UFF force field for optimization |
| `max_iterations` | `200` | `int` | Maximum optimization iterations per conformer (> 0) |

#### OpenEye Backend Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `mode` | `'classic'` | `str` | OMEGA generation mode (see options below) |
| `use_gpu` | `True` | `bool` | Enable GPU acceleration if available |
| `mpi_np` | `-1` | `int` | Number of MPI processes (-1 = auto-detect) |
| `strict` | `True` | `bool` | Use strict generation mode, fails SMILES with unassigned stereocenters |
| `flipper` | `False` | `bool` | Enable stereocenter flipping for unassigned stereocenters |
| `flipper_warts` | `False` | `bool` | Add suffix to flipped stereoisomers names (may cause issues for downstream tasks) |
| `flipper_maxcenters` | `4` | `int` | Maximum stereocenters to enumerate with flipper (> 0) |
| `oeomega_path` | `None` | `str` \| `None` | Path to OMEGA executable (None = auto-detect) |

#### OpenEye OMEGA Modes

- `'classic'` - Standard conformer generation, use this if unsure
- `'macrocycle'` - Optimized for macrocyclic molecules
- `'rocs'` - Optimized for ROCS overlay applications
- `'pose'` - Optimized for docking pose generation
- `'dense'` - Dense conformer sampling
- `'fastrocs'` - Fast ROCS-ready generation

#### Example

```python
# RDKit backend
confs_params = GenerateConfsParams(
    backend='rdkit',
    max_confs=200,
    rms_threshold=0.5,
    use_uff=True
)

# OpenEye backend
confs_params = GenerateConfsParams(
    backend='openeye',
    max_confs=500,
    mode='classic',
    use_gpu=True
)

# OpenEye for macrocycles
confs_params = GenerateConfsParams(
    backend='openeye',
    mode='macrocycle',
    max_confs=300,
    flipper=True
)
```

---

## Plugin Development

MolForge supports custom plugin actors for extending the pipeline with domain-specific functionality. Plugins are automatically discovered and integrate seamlessly with the core framework.

### Quick Start

Custom actors are placed in `molforge/actor_plugins/` and automatically loaded at runtime. Each plugin consists of a parameter class and an actor class:

```python
from dataclasses import dataclass
from molforge.actors.base import BaseActor
from molforge.actors.params.base import BaseParams

@dataclass
class MyPluginParams(BaseParams):
    my_parameter: str = 'default_value'

class MyPlugin(BaseActor):
    __step_name__ = 'my_plugin'  # Used in pipeline configuration
    __param_class__ = MyPluginParams

    def process(self, data):
        # Your processing logic here
        return data
```

### Using Plugins in Pipelines

Reference plugins by their `__step_name__` in the pipeline steps, and configure them using `plugin_params`:

```python
from molforge import MolForge, ForgeParams

params = ForgeParams(
    steps=['source', 'chembl', 'curate', 'my_plugin'],  # Add plugin to steps
    plugin_params={
        'my_plugin': MyPluginParams(
            my_parameter='custom_value'
        )
    }
)

forge = MolForge(params)
df = forge.forge("CHEMBL234")
```

### Common Patterns

#### Input and Output Specification

Declare required input columns and output columns for pipeline validation:

```python
class MyPlugin(BaseActor):
    @property
    def required_columns(self):
        return ['curated_smiles']  # Columns needed from previous actors

    @property
    def output_columns(self):
        return ['my_property']  # Columns added by this actor
```

#### Handling Optional Dependencies

Use graceful error handling for optional dependencies:

```python
@dataclass
class MyPluginParams(BaseParams):
    def _validate_params(self):
        try:
            import optional_library
        except ImportError:
            raise ImportError(
                "optional_library required. Install with: pip install optional_library"
            )
```

#### Pipeline Logging

Use `self.log()` for consistent logging integration:

```python
def process(self, data):
    self.log(f"Processing {len(data)} molecules")
    # ... processing logic ...
    self.log(f"Completed processing", level='DEBUG')
    return data
```

### Testing Plugins

Test plugins standalone before pipeline integration:

```python
import pandas as pd

# Create test data
test_data = pd.DataFrame({
    'curated_smiles': ['CCO', 'c1ccccc1', 'CC(=O)O']
})

# Initialize and test
params = MyPluginParams(my_parameter='test_value')
plugin = MyPlugin(params)
result = plugin.process(test_data)

print(result)
```

### Complete Example

See `molforge/actor_plugins/example.py` for a fully documented plugin template covering:
- Parameter validation with `_validate_params()`
- Actor initialization with `__post_init__()`
- Robust error handling for individual molecules
- Custom metadata in `_create_output()`
- Integration with pipeline configuration

The example plugin demonstrates calculating molecular descriptors with RDKit and serves as a reference for all plugin patterns.

---

## Complete Pipeline Example

```python
from molforge import MolForge, ForgeParams

# Configure complete pipeline
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
    tokens_params={
        'dynamically_update_vocab': True
    },

    # Property filtering
    distributions_params={
        'global_statistical_threshold': 2.0,
        'plot_distributions': True
    },

    # Conformer generation
    confs_params={
        'backend': 'rdkit',
        'max_confs': 200,
        'rms_threshold': 0.5
    }
)

# Execute pipeline
forge = MolForge(params)
df = forge.forge("CHEMBL234")
print(f"Processed {len(df)} molecules")
```

---

## Example Output

Running the complete pipeline example above generates organized output with logs, data, and visualizations:

### Output Structure

```text
MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/
├── CHEMBL234_ID60ef8cb9.csv          # Final curated dataset
├── CHEMBL234_ID60ef8cb9.log          # Execution log
├── CHEMBL234_ID60ef8cb9_config.json  # Pipeline configuration
├── curation_results.json             # Metadata and statistics
├── vocab.json                        # SMILES vocabulary
├── vocab_curated.json                # Post-curation vocabulary
└── distributions/                     # Property distribution plots
    ├── molecular_weight.png
    ├── logp.png
    ├── num_atoms.png
    ├── fsp3.png
    └── ... (18 property plots)
```

### Distribution Analysis

The pipeline automatically generates distribution plots for molecular properties, helping identify outliers and validate filtering thresholds:

![Property Distributions](docs/images/example_distributions.png)

<details>
<summary>📋 View complete execution log</summary>

```log
2025-12-15 17:56:53 | [  PIPELINE   ] | INFO | Logging to MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/CHEMBL234_ID60ef8cb9.log
2025-12-15 17:56:53 | [  PIPELINE   ] | INFO | Starting pipe.
2025-12-15 17:56:53 | [  PIPELINE   ] | INFO | Input (CHEMBL234): ChEMBL ID from ChEMBL ID: CHEMBL234 (0 rows)
2025-12-15 17:56:53 | [  PIPELINE   ] | INFO | [1/6] Starting source.
2025-12-15 17:56:53 | [     SQL     ] | INFO | 
==================================================
|   Fetching activity data for CHEMBL234	|
|   Database: ./data/chembl/36.db
==================================================
2025-12-15 17:56:54 | [     SQL     ] | INFO | Total entries available: 20595
2025-12-15 17:56:54 | [     SQL     ] | INFO | Fetching all 20595 entries...
2025-12-15 17:56:54 | [     SQL     ] | INFO | 20595 entries retrieved.
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [1/6] Finished source | success | 20595 rows | 0.45s
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [2/6] Starting chembl.
2025-12-15 17:56:54 | [   CHEMBL    ] | WARNING | RAW | Configured curation is not optimal.
2025-12-15 17:56:54 | [   CHEMBL    ] | WARNING | Recommended configuration for better data retention:
  standard_type   :	        'IC50' -> 'Ki'
  bao_format      :	 'BAO_0000357' -> 'BAO_0000219'
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | RAW
|    | standard_type   | assay_type   | target_organism   | bao_format   | entries    |   % |
|---:|:----------------|:-------------|:------------------|:-------------|:-----------|----:|
|  0 | Ki              | B            | Homo sapiens      | BAO_0000219  | 4684/20595 |  23 |
|  1 | AC50            | B            | Homo sapiens      | BAO_0000357  | 2128/20595 |  10 |
|  2 | Ki              | B            | Homo sapiens      | BAO_0000357  | 1750/20595 |   8 |
|  3 | Ki              | B            | Homo sapiens      | BAO_0000019  | 1414/20595 |   7 |
|  4 | Ki              | B            | Homo sapiens      | BAO_0000249  | 1243/20595 |   6 |
2025-12-15 17:56:54 | [   CHEMBL    ] | WARNING | VALID | Configured curation is not optimal.
2025-12-15 17:56:54 | [   CHEMBL    ] | WARNING | Recommended configuration for better data retention:
  standard_type   :	        'IC50' -> 'Ki'
  bao_format      :	 'BAO_0000357' -> 'BAO_0000219'
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | VALID
|    | standard_type   | assay_type   | target_organism   | bao_format   | entries   |   % |
|---:|:----------------|:-------------|:------------------|:-------------|:----------|----:|
|  0 | Ki              | B            | Homo sapiens      | BAO_0000219  | 3972/8588 |  46 |
|  1 | Ki              | B            | Homo sapiens      | BAO_0000357  | 1337/8588 |  16 |
|  2 | Ki              | B            | Homo sapiens      | BAO_0000019  | 1082/8588 |  13 |
|  3 | EC50            | F            | Homo sapiens      | BAO_0000219  | 427/8588  |   5 |
|  4 | AC50            | B            | Homo sapiens      | BAO_0000357  | 401/8588  |   5 |
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Conditional curation: 65/20595.
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Standardization: IC50 → pIC50 (log10 conversion).
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | 0 suspicious SMILES detected. mistakes_only=True, error_margin=0.0001
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Dropped 0 entries from 0 suspicious SMILES.
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Removed 1 (canonical) SMILES entries with std > 0.5 or range > 0.5.
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Aggregated 0 duplicate SMILES entries. Final: 61/65 rows.
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [2/6] Finished chembl | success | 61 rows | 0.27s
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [3/6] Starting curate.
2025-12-15 17:56:54 | [   CURATE    ] | INFO | Processing 61 molecules with single process.
2025-12-15 17:56:54 | [   CURATE    ] | INFO | Molecular curation: 61/61 successful. drop_na=True
2025-12-15 17:56:54 | [   CURATE    ] | INFO | Post-curation duplicates will be handled by ChEMBLCurator.handle_duplicate_smiles().
2025-12-15 17:56:54 | [   CHEMBL    ] | WARNING | Multiple ChEMBL IDs exist for (canonical) SMILES: 
| curated_smiles                                     |   pIC50 |     min |     max |      std |   count | molecule_chembl_ids                | document_years   |   earliest_year |   earliest_idx |
|:---------------------------------------------------|--------:|--------:|--------:|---------:|--------:|:-----------------------------------|:-----------------|----------------:|---------------:|
| C1=CC=C(CC[C@H]2CN(CC3=NC4=CC=CC=C4N3)CCO2)C=C1    | 4.30724 | 4.11691 | 4.49757 | 0.269172 |       2 | ['CHEMBL3335535', 'CHEMBL3335536'] | [2014.0, 2014.0] |            2014 |              0 |
| FC(F)(F)OC1=CC=C(CN2CCO[C@H](CCC3=CC=CC=C3)C2)C=C1 | 5.08778 | 5.08778 | 5.08778 | 0        |       2 | ['CHEMBL3335538', 'CHEMBL3335554'] | [2014.0, 2014.0] |            2014 |              0 |
2025-12-15 17:56:54 | [   CHEMBL    ] | WARNING | Taking first entry(s) for the names of curated data.
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Removed 0 (canonical) SMILES entries with std > 0.5 or range > 0.5.
2025-12-15 17:56:54 | [   CHEMBL    ] | INFO | Aggregated 2 duplicate SMILES entries. Final: 59/61 rows.
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [3/6] Finished curate | success | 59 rows | 0.06s
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [4/6] Starting tokens.
2025-12-15 17:56:54 | [   TOKENS    ] | INFO | Saved vocabulary to /home/luke/VScode/phd_luke_y1/TORSMILES/TORSMILES/MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/vocab.json
2025-12-15 17:56:54 | [   TOKENS    ] | INFO | Vocabulary constructed: 23 tokens: ['#', '$', '(', ')', '*', '/', '1', '2', '3', '4', '5', '6', '=', 'C', 'Cl', 'F', 'N', 'O', 'S', '[C@@H]', '[C@@]', '[C@H]', '^']
2025-12-15 17:56:54 | [   TOKENS    ] | INFO | Tokenized 59 sequences.
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [4/6] Finished tokens | success | 59 rows | 0.00s
2025-12-15 17:56:54 | [  PIPELINE   ] | INFO | [5/6] Starting distributions.
2025-12-15 17:56:54 | [DISTRIBUTIONS] | INFO | Using existing 'tokens' column.
2025-12-15 17:56:54 | [DISTRIBUTIONS] | INFO | Computing properties for 59 molecules (single process).
2025-12-15 17:56:54 | [DISTRIBUTIONS] | INFO | Computing token-based properties from 'tokens' column.
2025-12-15 17:56:55 | [DISTRIBUTIONS] | INFO | 
INITIAL DISTRIBUTION
| property              |   count |       mean |        std |        min |        max |     median |       skew |       kurt | symmetric [*]   | normal [**]   |   p_value | thresholds                              |
|:----------------------|--------:|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|:----------------|:--------------|----------:|:----------------------------------------|
| num_atoms             |      59 |  29.8475   |  6.07656   |  19        |  46        |  30        |  0.3301    | -0.494638  | ✓               | ✓             |  0.423    | ≥μ-2.0σ (17.69)      , ≤μ+2.0σ (42.00)  |
| num_rings             |      59 |   4.16949  |  0.722843  |   2        |   6        |   4        |  0.0139673 |  1.01656   | ✓               | ✓             |  0.231    | ≥μ-2.0σ (2.72)       , ≤μ+2.0σ (5.62)   |
| size_largest_ring     |      59 |   6.30508  |  0.464396  |   6        |   7        |   6        |  0.846642  | -1.2832    | ✓               | ✗             |  6.97e-07 | ≥μ-2.0σ (5.38)       , ≤μ+2.0σ (7.23)   |
| num_tokens            |      59 |  57.1017   | 12.4522    |  36        |  89        |  56        |  0.366139  | -0.525894  | ✓               | ✓             |  0.346    | ≥μ-2.0σ (32.20)      , ≤μ+2.0σ (82.01)  |
| tokens_atom_ratio     |      59 |   1.90969  |  0.0933192 |   1.69697  |   2.15789  |   1.89286  |  0.345711  |  0.0356023 | ✓               | ✓             |  0.455    | ≥μ-2.0σ (1.72)       , ≤μ+2.0σ (2.10)   |
| c_atom_ratio          |      59 |   0.779826 |  0.063982  |   0.612903 |   0.9      |   0.793103 | -0.818254  |  0.273873  | ✓               | ✗             |  0.0263   | ≥μ-2.0σ (0.65)       , ≤μ+2.0σ (0.91)   |
| longest_aliph_c_chain |      59 |   1.45763  |  1.27742   |   0        |   5        |   1        |  0.845997  |  0.186994  | ✓               | ✗             |  0.0241   | ≥μ-2.0σ (-1.10)      , ≤μ+2.0σ (4.01)   |
| molecular_weight      |      59 | 418.646    | 88.7198    | 256.257    | 644.754    | 428.536    |  0.124238  | -0.719739  | ✓               | ✓             |  0.351    | ≥μ-2.0σ (241.21)     , ≤μ+2.0σ (596.09) |
| logp                  |      59 |   4.19365  |  1.30618   |   1.2235   |   7.3578   |   3.9511   |  0.407158  | -0.419138  | ✓               | ✓             |  0.349    | ≥μ-2.0σ (1.58)       , ≤μ+2.0σ (6.81)   |
| tpsa                  |      59 |  52.7115   | 24.0964    |   6.48     | 116.87     |  51.37     |  0.465414  |  0.104584  | ✓               | ✓             |  0.255    | ≥μ-2.0σ (4.52)       , ≤μ+2.0σ (100.90) |
| num_rotatable_bonds   |      59 |   5.18644  |  3.04265   |   0        |  12        |   5        |  0.154515  | -0.902814  | ✓               | ✓             |  0.0926   | ≥μ-2.0σ (-0.90)      , ≤μ+2.0σ (11.27)  |
| num_h_donors          |      59 |   0.898305 |  0.75874   |   0        |   3        |   1        |  0.646883  |  0.309195  | ✓               | ✓             |  0.0756   | ≥μ-2.0σ (-0.62)      , ≤μ+2.0σ (2.42)   |
| num_h_acceptors       |      59 |   4.52542  |  1.75535   |   2        |  10        |   4        |  0.883331  |  0.6375    | ✓               | ✗             |  0.0101   | ≥μ-2.0σ (1.01)       , ≤μ+2.0σ (8.04)   |
| fsp3                  |      59 |   0.343616 |  0.103207  |   0        |   0.7      |   0.346154 |  0.39271   |  4.09819   | ✗               | ✗             |  0.000967 | ≥μ-2.0σ (0.14)       , ≤μ+2.0σ (0.55)   |
| num_aromatic_rings    |      59 |   2.62712  |  0.66691   |   1        |   4        |   3        |  0.229582  | -0.41293   | ✓               | ✓             |  0.658    | ≥μ-2.0σ (1.29)       , ≤μ+2.0σ (3.96)   |
| num_stereocenters     |      59 |   0.576271 |  1.13264   |   0        |   6        |   0        |  2.68162   |  8.3679    | ✗               | ✗             |  1.47e-12 | ≥μ-2.0σ (-1.69)      , ≤μ+2.0σ (2.84)   |
| num_heteroatoms       |      59 |   6.67797  |  2.5626    |   2        |  13        |   7        |  0.51764   | -0.339115  | ✓               | ✓             |  0.222    | ≥μ-2.0σ (1.55)       , ≤μ+2.0σ (11.80)  |
| heteroatom_ratio      |      59 |   0.220174 |  0.063982  |   0.1      |   0.387097 |   0.206897 |  0.818254  |  0.273873  | ✓               | ✗             |  0.0263   | ≥μ-2.0σ (0.09)       , ≤μ+2.0σ (0.35)   |
2025-12-15 17:56:55 | [DISTRIBUTIONS] | INFO | [*] Practical normality: Skewness and Kurtosis tests. ✓ = approximately normal (|skew|<1, |kurt|<2); ✗ = notably non-normal.
2025-12-15 17:56:55 | [DISTRIBUTIONS] | INFO | [**] Normality test: D'Agostino-Pearson test. p < 0.05 suggests non-normal distribution.
2025-12-15 17:56:55 | [DISTRIBUTIONS] | INFO | Generating distribution plots with 59 molecules...
2025-12-15 17:56:55 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_atoms.png
2025-12-15 17:56:55 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_rings.png
2025-12-15 17:56:56 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/size_largest_ring.png
2025-12-15 17:56:56 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_tokens.png
2025-12-15 17:56:56 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/tokens_atom_ratio.png
2025-12-15 17:56:56 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/c_atom_ratio.png
2025-12-15 17:56:57 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/longest_aliph_c_chain.png
2025-12-15 17:56:57 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/molecular_weight.png
2025-12-15 17:56:57 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/logp.png
2025-12-15 17:56:58 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/tpsa.png
2025-12-15 17:56:58 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_rotatable_bonds.png
2025-12-15 17:56:58 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_h_donors.png
2025-12-15 17:56:58 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_h_acceptors.png
2025-12-15 17:56:59 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/fsp3.png
2025-12-15 17:56:59 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_aromatic_rings.png
2025-12-15 17:56:59 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_stereocenters.png
2025-12-15 17:56:59 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/num_heteroatoms.png
2025-12-15 17:57:00 | [DISTRIBUTIONS] | DEBUG |   Saved plot: MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/distributions/heteroatom_ratio.png
2025-12-15 17:57:00 | [DISTRIBUTIONS] | INFO | 
FILTER RESULTS
| property              | passed   |   removed | pass_rate   |
|:----------------------|:---------|----------:|:------------|
| num_aromatic_rings    | 53/59    |         6 | 89.8%       |
| heteroatom_ratio      | 55/59    |         4 | 93.2%       |
| c_atom_ratio          | 55/59    |         4 | 93.2%       |
| fsp3                  | 55/59    |         4 | 93.2%       |
| logp                  | 56/59    |         3 | 94.9%       |
| num_heteroatoms       | 56/59    |         3 | 94.9%       |
| num_stereocenters     | 56/59    |         3 | 94.9%       |
| num_rings             | 56/59    |         3 | 94.9%       |
| tpsa                  | 56/59    |         3 | 94.9%       |
| num_h_donors          | 57/59    |         2 | 96.6%       |
| num_h_acceptors       | 57/59    |         2 | 96.6%       |
| longest_aliph_c_chain | 57/59    |         2 | 96.6%       |
| tokens_atom_ratio     | 57/59    |         2 | 96.6%       |
| num_tokens            | 57/59    |         2 | 96.6%       |
| num_atoms             | 57/59    |         2 | 96.6%       |
| molecular_weight      | 58/59    |         1 | 98.3%       |
| num_rotatable_bonds   | 58/59    |         1 | 98.3%       |
| size_largest_ring     | 59/59    |         0 | 100.0%      |
2025-12-15 17:57:00 | [DISTRIBUTIONS] | INFO | Distribution filters: 42/59 molecules retained (17 removed, 71.2% overall pass rate).
2025-12-15 17:57:00 | [DISTRIBUTIONS] | INFO | Distribution curation complete: 42 molecules retained.
2025-12-15 17:57:00 | [DISTRIBUTIONS] | INFO | Curation results saved to MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/curation_results.json
2025-12-15 17:57:00 | [DISTRIBUTIONS] | INFO | Curated vocabulary saved to MOLFORGE_OUT/CHEMBL234_ID60ef8cb9/vocab_curated.json
  Original: 23 tokens
  Curated: 22 tokens
  Removed: 1 tokens
2025-12-15 17:57:00 | [  PIPELINE   ] | INFO | [5/6] Finished distributions | success | 42 rows | 5.53s
2025-12-15 17:57:00 | [  PIPELINE   ] | INFO | [6/6] Starting confs.
2025-12-15 17:57:00 | [  CONFS-RDK  ] | INFO | Generating conformers for 42 SMILES
2025-12-15 17:57:00 | [     RDK     ] | INFO | Generating conformers for 42 molecules with 19 processes (14 chunks of 3).
2025-12-15 17:57:12 | [     RDK     ] | INFO | Progress: 1/14 chunks (7.1%) | Elapsed: 11.7s | ETA: 152.5s | Rate: 0 mol/s
2025-12-15 17:57:12 | [     RDK     ] | INFO | Progress: 2/14 chunks (14.3%) | Elapsed: 11.9s | ETA: 71.3s | Rate: 1 mol/s
2025-12-15 17:57:15 | [     RDK     ] | INFO | Progress: 3/14 chunks (21.4%) | Elapsed: 15.6s | ETA: 57.2s | Rate: 1 mol/s
2025-12-15 17:57:15 | [     RDK     ] | INFO | Progress: 4/14 chunks (28.6%) | Elapsed: 15.6s | ETA: 39.0s | Rate: 1 mol/s
2025-12-15 17:57:19 | [     RDK     ] | INFO | Progress: 5/14 chunks (35.7%) | Elapsed: 19.2s | ETA: 34.5s | Rate: 1 mol/s
2025-12-15 17:57:20 | [     RDK     ] | INFO | Progress: 6/14 chunks (42.9%) | Elapsed: 20.6s | ETA: 27.5s | Rate: 1 mol/s
2025-12-15 17:57:23 | [     RDK     ] | INFO | Progress: 7/14 chunks (50.0%) | Elapsed: 23.1s | ETA: 23.1s | Rate: 1 mol/s
2025-12-15 17:57:24 | [     RDK     ] | INFO | Progress: 8/14 chunks (57.1%) | Elapsed: 24.3s | ETA: 18.2s | Rate: 1 mol/s
2025-12-15 17:57:25 | [     RDK     ] | INFO | Progress: 9/14 chunks (64.3%) | Elapsed: 25.0s | ETA: 13.9s | Rate: 1 mol/s
2025-12-15 17:57:27 | [     RDK     ] | INFO | Progress: 10/14 chunks (71.4%) | Elapsed: 27.4s | ETA: 10.9s | Rate: 1 mol/s
2025-12-15 17:57:37 | [     RDK     ] | INFO | Progress: 11/14 chunks (78.6%) | Elapsed: 37.2s | ETA: 10.1s | Rate: 1 mol/s
2025-12-15 17:57:39 | [     RDK     ] | INFO | Progress: 12/14 chunks (85.7%) | Elapsed: 38.8s | ETA: 6.5s | Rate: 1 mol/s
2025-12-15 17:57:42 | [     RDK     ] | INFO | Progress: 13/14 chunks (92.9%) | Elapsed: 42.2s | ETA: 3.2s | Rate: 1 mol/s
2025-12-15 17:57:49 | [     RDK     ] | INFO | Progress: 14/14 chunks (100.0%) | Elapsed: 48.8s | ETA: 0.0s | Rate: 1 mol/s
2025-12-15 17:57:49 | [     RDK     ] | INFO | RDKit generation complete: 42/42 succeeded
2025-12-15 17:57:49 | [  CONFS-RDK  ] | INFO | Conformer generation complete: 42 succeeded, 0 failed, avg 121.0 conformers/molecule
2025-12-15 17:57:49 | [  PIPELINE   ] | INFO | [6/6] Finished confs | success | 42 rows | 48.83s
2025-12-15 17:57:49 | [  PIPELINE   ] | INFO | Pipeline completed | Total time: 55.15s
```
</details>

<details>
<summary>📋 View written configuration file (.json) </summary>
  
```yaml
{
    "verbose": true,
    "steps": [
        "source",
        "chembl",
        "curate",
        "tokens",
        "distributions",
        "confs"
    ],
    "return_none_on_fail": false,
    "output_root": "MOLFORGE_OUT",
    "write_output": true,
    "write_checkpoints": false,
    "console_log_level": "INFO",
    "file_log_level": "DEBUG",
    "override_actor_params": true,
    "source_params": {
        "verbose": true,
        "backend": "sql",
        "n": 1000,
        "search_all": true,
        "db_path": "./data/chembl/36.db",
        "version": 36,
        "auto_download": true,
        "download_dir": "./data/chembl"
    },
    "chembl_params": {
        "verbose": true,
        "standard_type": "IC50",
        "standard_units": "nM",
        "standard_relation": "=",
        "assay_type": "B",
        "target_organism": "Homo sapiens",
        "assay_format": "protein",
        "top_n": 5,
        "mutant_regex": "mutant|mutation|variant",
        "allosteric_regex": "allosteric",
        "C": 1000000000.0,
        "mistakes_only": true,
        "error_margin": 0.0001,
        "SMILES_column": "canonical_smiles",
        "std_threshold": 0.5,
        "range_threshold": 0.5,
        "bao_format": "BAO_0000357"
    },
    "curate_params": {
        "verbose": true,
        "mol_steps": [
            "desalt",
            "neutralize",
            "sanitize",
            "handleStereo"
        ],
        "smiles_steps": [
            "canonical",
            "kekulize"
        ],
        "neutralize_policy": "keep",
        "desalt_policy": "keep",
        "brute_force_desalt": false,
        "max_tautomers": 512,
        "step_timeout": 60,
        "SMILES_column": "canonical_smiles",
        "dropna": true,
        "check_duplicates": true,
        "duplicates_policy": "first",
        "stereo_policy": "assign",
        "assign_policy": "random",
        "random_seed": 42,
        "try_embedding": false,
        "only_unassigned": true,
        "only_unique": true,
        "max_isomers": 32,
        "stereo_params": {
            "stereo_policy": "assign",
            "assign_policy": "random",
            "max_isomers": 32,
            "try_embedding": false,
            "only_unassigned": true,
            "only_unique": true,
            "random_seed": 42
        },
        "desalt": true,
        "removeIsotope": false,
        "removeHs": false,
        "tautomers": false,
        "neutralize": true,
        "sanitize": true,
        "computeProps": false,
        "handleStereo": true,
        "canonical": true,
        "kekulize": true
    },
    "tokens_params": {
        "verbose": true,
        "vocab_file": null,
        "SMILES_column": "curated_smiles",
        "dynamically_update_vocab": true
    },
    "distributions_params": {
        "verbose": true,
        "SMILES_column": "curated_smiles",
        "tokens_column": "tokens",
        "thresholds": {
            "num_atoms": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_rings": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "size_largest_ring": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_tokens": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "tokens_atom_ratio": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "c_atom_ratio": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "longest_aliph_c_chain": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "molecular_weight": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "logp": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "tpsa": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_rotatable_bonds": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_h_donors": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_h_acceptors": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "fsp3": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_aromatic_rings": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_stereocenters": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "num_heteroatoms": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            },
            "heteroatom_ratio": {
                "verbose": true,
                "min_value": null,
                "max_value": null,
                "statistical_lower": -2.0,
                "statistical_upper": 2.0,
                "quantile_lower": null,
                "quantile_upper": null
            }
        },
        "global_statistical_threshold": 2.0,
        "global_quantile_threshold": null,
        "properties": "all",
        "curate_tokens": true,
        "filter_unknown_tokens": true,
        "token_frequency_threshold": null,
        "compute_properties": true,
        "dropna": true,
        "plot_distributions": true,
        "perform_pca": false
    },
    "confs_params": {
        "verbose": true,
        "backend": "rdkit",
        "max_confs": 200,
        "rms_threshold": 0.5,
        "SMILES_column": "curated_smiles",
        "names_column": "molecule_chembl_id",
        "dropna": true,
        "use_random_coords": true,
        "random_seed": 42,
        "num_threads": 0,
        "use_uff": true,
        "max_iterations": 200,
        "mode": "classic",
        "use_gpu": true,
        "mpi_np": -1,
        "strict": true,
        "flipper": false,
        "flipper_warts": false,
        "flipper_maxcenters": 4,
        "oeomega_path": null,
        "backend_kwargs": {}
    },
    "plugin_params": {},
    "log_time": "2025-12-15 17:56:53",
    "input_type": "ChEMBL ID",
    "input_id": "CHEMBL234",
    "input_source": "ChEMBL ID: CHEMBL234",
    "run_id": "CHEMBL234_ID60ef8cb9"
}
```
  
</details>

---

## Command Line Interface

```bash
# Run pipeline with configuration file
molforge run CHEMBL234 --config config.yaml

# Get actor information
molforge info actors

# Show architecture
molforge info architecture
```

---

## Contributing

This is a public beta release. Contributions welcome through [GitHub Issues](https://github.com/molML/molforge/issues).

## License

MIT License
