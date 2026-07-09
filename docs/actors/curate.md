# CurateMol (`curate`)

[← Back to README](../../README.md) · [Documentation index](../../README.md#documentation)

Standardizes molecules through configurable curation steps including desalting, neutralization, stereochemistry handling, and SMILES canonicalization.

Typically consumes the output of the [`chembl` actor](chembl.md) and feeds the [`tokens` actor](tokens.md).

## Core Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `mol_steps` | See below | `List[str]` | Ordered list of mol-level curation steps to apply |
| `smiles_steps` | `['canonical', 'kekulize']` | `List[str]` | SMILES-level steps, handled separately and applied at the end of the pipe |
| `SMILES_column` | `'canonical_smiles'` | `str` | DataFrame column containing SMILES strings |
| `dropna` | `True` | `bool` | Drop molecules that failed curation, keeping only rows where `curation_success` is `True` |
| `check_duplicates` | `True` | `bool` | Check for and remove duplicate curated SMILES (post-curation) |
| `duplicates_policy` | `'first'` | `'first'` \| `'last'` \| `False` | Which duplicate to keep when dropping duplicate curated SMILES (passed to pandas `drop_duplicates`). Only applied in the fallback path when no `ChEMBLCurator` actor is available |

Default `mol_steps`: `['desalt', 'removeIsotope', 'removeHs', 'tautomers', 'neutralize', 'sanitize', 'handleStereo', 'computeProps']`

## Step-Specific Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `neutralize_policy` | `'keep'` | `'keep'` \| `'remove'` | `'keep'` returns already-neutral molecules unchanged without running neutralization reactions; `'remove'` always runs them (also neutralizing e.g. zwitterions). Charged molecules are neutralized under both; those that cannot reach zero charge fail the step |
| `desalt_policy` | `'keep'` | `'keep'` \| `'remove'` | `'keep'` strips salts from multi-fragment molecules and keeps the desalted molecule; `'remove'` drops any multi-fragment molecule without desalting. Single-fragment molecules pass unchanged |
| `brute_force_desalt` | `False` | `bool` | Use brute-force desalting to strip all fragments but the largest |
| `max_tautomers` | `512` | `int` | Maximum number of tautomers to enumerate per molecule |
| `step_timeout` | `60` | `int` | Per-step timeout in seconds per molecule (mainly for very large/complex molecules and CIP assignment) |

## Stereochemistry Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `stereo_policy` | `'assign'` | `'keep'` \| `'remove'` \| `'assign'` \| `'enumerate'` | How to handle stereochemistry: keep, remove, assign, or enumerate stereoisomers |
| `assign_policy` | `'random'` | `'first'` \| `'random'` \| `'lowest'` | Which stereoisomer to assign when `stereo_policy` is `'assign'` |
| `random_seed` | `42` | `int` | Random seed for reproducible stereoisomer assignment/enumeration |
| `try_embedding` | `False` | `bool` | Passed to RDKit `StereoEnumerationOptions.tryEmbedding`; when `True`, enumerated stereoisomers that cannot be 3D-embedded are discarded |
| `only_unassigned` | `True` | `bool` | Only enumerate/assign stereochemistry at unassigned stereocenters |
| `only_unique` | `True` | `bool` | Passed to RDKit's `StereoEnumerationOptions(unique=...)`; when `True`, duplicate enumerated stereoisomers are filtered out (applies when `stereo_policy='enumerate'`) |
| `max_isomers` | `32` | `int` | Maximum number of stereoisomers to enumerate per molecule |

## Available Curation Steps

**Molecule steps** (`mol_steps`): `desalt`, `removeIsotope`, `removeHs`, `tautomers`, `neutralize`, `sanitize`, `handleStereo`, `computeProps`

**SMILES steps** (`smiles_steps`): `canonical`, `kekulize`

## Example

```python
from molforge import CurateMolParams

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

## Related

- [Base Parameters](../configuration.md#base-parameters) inherited by all actors
- [Pipeline configuration](../configuration.md) — how `curate_params` is set
- [Using this actor standalone](../standalone-actors.md)
- Previous step: [ChEMBLCurator (`chembl`)](chembl.md)
- Next step: [TokenizeData (`tokens`)](tokens.md)
