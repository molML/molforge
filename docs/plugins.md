# Plugin Development

[← Back to README](../README.md) · [Documentation index](../README.md#documentation)

MolForge supports custom plugin actors for extending the pipeline with domain-specific functionality. Plugins are automatically discovered and integrate seamlessly with the core framework.

## Discovery Rules

Plugins are auto-discovered from `molforge/actor_plugins/`. At load time the registry scans every `*.py` file in that directory (files starting with `_` are skipped) and registers any class that:

- defines `__step_name__` (the name used in pipeline `steps`), and
- defines `__param_class__` (its `BaseParams` subclass), and
- is callable.

A class may optionally declare `__dependencies__ = [...]`, a list of step names that must appear **earlier** in the pipeline; the registry validates this ordering.

## Writing a Plugin

Each plugin consists of a parameter dataclass and an actor class:

```python
from dataclasses import dataclass
from typing import List
import pandas as pd

from molforge.actors.base import BaseActor
from molforge.actors.params.base import BaseParams


@dataclass
class MyPluginParams(BaseParams):
    my_parameter: str = 'default_value'
    """Docstring describing this parameter (surfaced by `molforge info`)."""

    def _validate_params(self) -> None:
        # `_validate_params` is abstract on BaseParams and MUST be implemented.
        if not self.my_parameter:
            raise ValueError("my_parameter cannot be empty")


class MyPlugin(BaseActor):
    __step_name__ = 'my_plugin'          # Used in pipeline configuration
    __param_class__ = MyPluginParams
    # __dependencies__ = ['curate']      # optional: require an upstream step

    @property
    def required_columns(self) -> List[str]:
        return ['curated_smiles']        # columns needed from upstream actors

    @property
    def output_columns(self) -> List[str]:
        return ['my_property']           # columns this actor adds

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        self.log(f"Processing {len(data)} molecules")
        # ... processing logic ...
        return data
```

## Using Plugins in Pipelines

Reference a plugin by its `__step_name__` in `steps`, and configure it through `plugin_params`:

```python
from molforge import MolForge, ForgeParams

params = ForgeParams(
    steps=['source', 'chembl', 'curate', 'my_plugin'],
    plugin_params={'my_plugin': MyPluginParams(my_parameter='custom_value')}
)

forge = MolForge(params)
df = forge.forge("CHEMBL234")
```

See [Pipeline Configuration](configuration.md) for how `plugin_params` fits alongside the core `{step}_params` attributes.

## Handling Optional Dependencies

Validate optional third-party libraries in `_validate_params` so failures surface early:

```python
def _validate_params(self) -> None:
    try:
        import optional_library  # noqa: F401
    except ImportError:
        raise ImportError("optional_library required. Install with: pip install optional_library")
```

Use `self.log(...)` for consistent logging integration (`level='DEBUG'` for verbose detail).

## Template

See `molforge/actor_plugins/example.py` for a fully documented plugin template (the `properties` plugin below). It demonstrates parameter validation with `_validate_params()`, actor-specific setup in `__post_init__()`, robust per-molecule error handling, and custom metadata in `_create_output()`.

## Included Plugins

MolForge ships three demo plugins in `molforge/actor_plugins/`:

### `properties` — RDKit descriptors (`example.py`)

Computes RDKit molecular descriptors and adds one column per descriptor. Serves as the reference plugin template.

| Key Parameter | Default | Description |
|---------------|---------|-------------|
| `properties` | `None` → common set | RDKit descriptor names (from `rdkit.Chem.Descriptors`); defaults to `ExactMolWt`, `MolLogP`, `TPSA`, `NumHDonors`, `NumHAcceptors` |
| `smiles_column` | `'curated_smiles'` | Column containing the SMILES to process |
| `filter_invalid` | `False` | Drop molecules whose descriptor calculation failed |

**Output columns:** one per requested descriptor.

### `scaffold` — Bemis-Murcko scaffolds (`scaffold.py`)

Computes Bemis-Murcko scaffolds (and optionally generic carbon-skeleton scaffolds) for every molecule. Requires only RDKit. Declares `__dependencies__ = ['curate']`.

| Key Parameter | Default | Description |
|---------------|---------|-------------|
| `SMILES_column` | `'curated_smiles'` | Input SMILES column (output of `curate`) |
| `include_generic` | `True` | Also compute the generic (carbon-skeleton) scaffold |
| `include_chirality` | `False` | Preserve stereo annotations in the Murcko scaffold SMILES |
| `dropna` | `False` | Remove rows where `scaffold_success` is `False` |
| `acyclic_policy` | `'keep'` | `'keep'` retains acyclic molecules (empty scaffold); `'remove'` drops them |

**Output columns:** `scaffold_smiles`, `scaffold_generic_smiles`, `scaffold_success`.

### `split` — train/val/test split (`split.py`)

Partitions the dataset into train/val/test along two orthogonal axes. `unit` sets the atom of assignment: `scaffold` groups molecules by Bemis-Murcko scaffold and assigns whole clusters together, so molecules sharing a scaffold share a split; `molecule` treats each molecule independently. `method` sets the ordering used to fill test → val → train: `isolation` ranks units by weighted-mean ECFP4 Tanimoto distance to the rest of the dataset and fills the most isolated first, producing a chemically out-of-distribution test set; `random` orders units by a seeded permutation. Writes a JSON report card and a PNG report figure. Declares `__dependencies__ = ['scaffold']`.

| Key Parameter | Default | Description |
|---------------|---------|-------------|
| `unit` | `'scaffold'` | Atom of assignment: `scaffold` (whole clusters together) or `molecule` |
| `method` | `'isolation'` | Unit ordering: `isolation` (structural distance) or `random` (seeded permutation) |
| `scaffold_column` | `'scaffold_smiles'` | Scaffold column from the `scaffold` actor |
| `test_ratio` | `0.10` | Fraction of molecules assigned to the test set |
| `val_ratio` | `0.10` | Fraction assigned to validation (train receives the remainder) |
| `seed` | `42` | Seed for the `random` method's permutation |
| `max_units_for_isolation` | `10000` | Compute bound for `method='isolation'`; a larger unit count raises `ValueError` |
| `activity_column` | `None` | Continuous activity label for report diagnostics (auto-resolved from the curator's endpoint when `None`) |
| `ecfp4_radius` / `ecfp4_n_bits` | `2` / `2048` | Morgan fingerprint parameters |
| `max_nn_dist_n` | `10000` | Max dataset size for molecule-level ECFP4 nearest-neighbour distance stats |

**Output column:** `split` (`'train'` / `'val'` / `'test'`).

The [`scaffold` + `split` plugins are exercised in the example run](example.md), which shows the split report figure.
