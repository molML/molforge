# ChEMBLCurator (`chembl`)

[← Back to README](../../README.md) · [Documentation index](../../README.md#documentation)

Curates ChEMBL bioactivity data with automatic type-driven behavior. Curation logic adapts based on the `standard_type` category (potency, kinetic, ADMET, activity).

Typically consumes the output of the [`source` actor](source.md) and feeds the [`curate` actor](curate.md).

## Core Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `standard_type` | `'IC50'` | `str` | ChEMBL `standard_type` to curate (drives type-specific curation behavior) |
| `standard_units` | `None` | `str` \| `None` | Target `standard_units`; auto-set from the type registry when `None` |
| `standard_relation` | `'='` | `str` \| `None` | Activity relation to keep (e.g. `'='`, `'>'`, `'<'`) |
| `assay_type` | `'B'` | `'B'` \| `'F'` \| `'A'` \| `'T'` | ChEMBL assay type: `'B'` (binding), `'F'` (functional), `'A'` (ADME), or `'T'` (toxicity) |
| `target_organism` | `'Homo sapiens'` | `str` | Target organism to filter assays by |
| `assay_format` | `'protein'` | `'general'` \| `'protein'` \| `'cell'` \| `'organism'` \| `'tissue'` \| `'microsome'` | Assay format; mapped to a BAO format code |

## Quality Control Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `top_n` | `5` | `int` | Number of top experimental-condition groupings shown in the logged RAW/VALID analysis tables (reporting only; does not affect curation) |
| `mistakes_only` | `True` | `bool` | If `True`, flag a SMILES as suspicious only when its repeated standardized values differ by ~3 log units (or ~1 for Ki/Kd); if `False`, also flag SMILES whose repeated values are (near-)identical. Flagged SMILES are dropped |
| `error_margin` | `0.0001` | `float` | Numerical tolerance used when comparing values for mistake detection |
| `std_threshold` | `0.5` | `float` | Maximum per-SMILES standard deviation of standardized values for a duplicate group to be kept during aggregation |
| `range_threshold` | `0.5` | `float` | Maximum per-SMILES value range (max − min of standardized values) for a duplicate group to be kept during aggregation |

## Advanced Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `mutant_regex` | `'mutant\|mutation\|variant'` | `str` | Regex used to detect and exclude mutant/variant assays |
| `allosteric_regex` | `'allosteric'` | `str` | Regex used to detect and exclude allosteric assays |
| `C` | `1e9` | `float` | Molarity conversion constant, auto-set from `standard_units`. Only change for non-standard unit conversion |
| `SMILES_column` | `'canonical_smiles'` | `str` | DataFrame column containing SMILES strings |
| `bao_format` | `None` (auto-set) | `str` \| `None` | Derived BioAssay Ontology format code; overwritten from `assay_format` (do not set manually) |

## Supported Standard Types

- **Potency**: `IC50`, `EC50`, `AC50`, `Ki`, `Kd`, `XC50`
- **Kinetic**: `kon`, `k_off`, `koff`, `Kon`, `Koff`
- **ADMET**: `Log BB`, `Log D`, `Log P`, `Clearance`, `T1/2`, `Solubility`, `Permeability`, `Caco-2`, `MDCK`
- **Activity**: `Inhibition`, `Percent Effect`, `Activity`, `Potency`

Use `ChEMBLCuratorParams.list_supported_types()` to list all supported types.

## Example

```python
from molforge import ChEMBLCuratorParams

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

## Related

- [Base Parameters](../configuration.md#base-parameters) inherited by all actors
- [Pipeline configuration](../configuration.md) — how `chembl_params` is set
- [Using this actor standalone](../standalone-actors.md)
- Previous step: [ChEMBLSource (`source`)](source.md)
- Next step: [CurateMol (`curate`)](curate.md)
