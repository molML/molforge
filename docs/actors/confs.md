# GenerateConfs (`confs`)

[← Back to README](../../README.md) · [Documentation index](../../README.md#documentation)

Generates 3D conformer ensembles for molecules. Supports RDKit and OpenEye backends with backend-specific optimization parameters.

Typically consumes the output of the [`distributions` actor](distributions.md) (or any actor providing curated SMILES).

## Shared Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `backend` | `'rdkit'` | `'rdkit'` \| `'openeye'` | Conformer generation backend |
| `max_confs` | `200` | `int` | Maximum number of conformers to generate per molecule |
| `rms_threshold` | `0.5` | `float` | RMS threshold for conformer pruning (Angstroms) |
| `SMILES_column` | `'curated_smiles'` | `str` | DataFrame column containing SMILES strings |
| `names_column` | `'molecule_chembl_id'` | `str` | DataFrame column containing molecule identifiers |
| `dropna` | `True` | `bool` | Drop rows with missing SMILES before processing |
| `convert_to_rdkit` | `False` | `bool` | Convert conformers to RDKit `Mol` objects when extracting (OpenEye only; RDKit is native) |
| `timeout` | `3600` | `int` | Timeout for conformer generation execution, in seconds |

## RDKit Backend Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `use_random_coords` | `True` | `bool` | Use random coordinates for initial embedding (RDKit) |
| `random_seed` | `42` | `int` | Random seed for reproducibility (RDKit) |
| `use_uff` | `True` | `bool` | Use UFF force field for conformer optimization (RDKit) |
| `max_iterations` | `200` | `int` | Maximum optimization iterations per conformer (RDKit) |

## OpenEye Backend Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `mode` | `'classic'` | `'classic'` \| `'macrocycle'` \| `'rocs'` \| `'pose'` \| `'dense'` \| `'fastrocs'` | OMEGA generation mode (OpenEye) |
| `use_gpu` | `True` | `bool` | Enable GPU acceleration if available (OpenEye) |
| `mpi_np` | `-1` | `int` | Number of MPI processes, `-1` = auto-detect (OpenEye) |
| `strict` | `True` | `bool` | Only process molecules with complete (stereochemistry) information (OpenEye) |
| `flipper` | `False` | `bool` | Enable unassigned stereo-isomer flipping (OpenEye) |
| `flipper_warts` | `False` | `bool` | Add suffix to flipped stereoisomer names (OpenEye) |
| `flipper_maxcenters` | `4` | `int` | Maximum stereocenters to enumerate with flipper (OpenEye) |
| `oeomega_path` | `None` | `str` \| `None` | Path to OMEGA executable (OpenEye); `None` = auto-detect |

## Example

```python
from molforge import GenerateConfsParams

# RDKit backend
confs_params = GenerateConfsParams(backend='rdkit', max_confs=200, rms_threshold=0.5, use_uff=True)

# OpenEye backend
confs_params = GenerateConfsParams(backend='openeye', max_confs=500, mode='classic', use_gpu=True)

# OpenEye for macrocycles
confs_params = GenerateConfsParams(backend='openeye', mode='macrocycle', max_confs=300, flipper=True)
```

Using the OpenEye backend requires the OpenEye Toolkit — see [Installation](../../README.md#installation).

## Related

- [Base Parameters](../configuration.md#base-parameters) inherited by all actors
- [Pipeline configuration](../configuration.md) — how `confs_params` is set
- [Using this actor standalone](../standalone-actors.md) (pass a `PipelineContext` to control output paths)
- Previous step: [CurateDistribution (`distributions`)](distributions.md)
