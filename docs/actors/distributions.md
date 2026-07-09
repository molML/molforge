# CurateDistribution (`distributions`)

[← Back to README](../../README.md) · [Documentation index](../../README.md#documentation)

Filters molecules based on molecular property distributions using statistical or quantile thresholds. A global threshold can be set for all specified properties for ease of use. Note that distribution thresholds apply to the set of molecules at that stage of the pipeline; these are typically non-normal and low sample size (~1k mols). For large-scale curation it is best to calibrate thresholds on a larger chemical space.

Typically consumes the output of the [`tokens` actor](tokens.md) (it reuses that actor's vocabulary for token curation).

## Output columns: `distribution_success` / `distribution_failures`

Every molecule is flagged with two output columns: `distribution_success` (`bool`) is `True` only when the molecule passes **all** active property thresholds **and** the token filters, and `distribution_failures` (`str`) lists the violated filters for failing molecules as a semicolon-separated string (empty for molecules that pass). Reasons use the form `"<property>:low"` / `"<property>:high"` for a value below/above a property bound, `"tokens:unknown"` for out-of-vocabulary tokens, and `"tokens:rare"` for rare tokens — e.g. `"logp:high; num_atoms:low; tokens:rare"`. By default (`dropna=True`) failing molecules are dropped and survivors carry `distribution_success=True` with an empty `distribution_failures`; set `dropna=False` to retain every molecule and use these columns to inspect or post-filter the pass/fail decisions yourself.

## Core Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `SMILES_column` | `'curated_smiles'` | `str` | DataFrame column containing SMILES strings |
| `tokens_column` | `'tokens'` | `str` | DataFrame column containing tokenized SMILES |
| `properties` | `None` | `List[str]` \| `'all'` \| `None` | Which properties to compute and filter on: `'all'` for every supported property, an explicit list of names, or `None` to use only the properties that have an active threshold configured |
| `compute_properties` | `True` | `bool` | Compute molecular properties before applying distribution filters |
| `dropna` | `True` | `bool` | If `True`, drop molecules that fail any active distribution/token filter; if `False`, retain all molecules and record pass/fail in `distribution_success` / `distribution_failures` |

## Threshold Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `thresholds` | `{}` | `Dict[str, PropertyThreshold]` | Per-property threshold configuration, keyed by property name |
| `global_statistical_threshold` | `None` | `float` \| `None` | Statistical threshold applied to all properties lacking explicit thresholds (e.g. `2.0` for ±2σ) |
| `global_quantile_threshold` | `None` | `float` \| `None` | Quantile threshold applied to all properties lacking explicit thresholds (e.g. `0.025` for 2.5%–97.5%) |

## Token Curation Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `curate_tokens` | `True` | `bool` | Enable token-based curation of molecules |
| `filter_unknown_tokens` | `True` | `bool` | Remove molecules containing tokens absent from the vocabulary (requires `curate_tokens=True`) |
| `token_frequency_threshold` | `None` | `float` \| `None` | Rarity cutoff as a percentage (strictly between 0 and 100): tokens present in fewer than this percent of molecules are "rare", and molecules containing any rare token are removed (requires `curate_tokens=True`) |

## Vocabulary handling

Two distinct vocabularies are involved:

- **Base vocabulary** — owned by the [`tokens` actor](tokens.md). Either built fresh from the input data or loaded from `vocab_file`, and grown with new tokens only when `dynamically_update_vocab=True`. It is the reference token set for the run.
- **Curated vocabulary** (`vocab_curated.json`) — produced by `distributions`. It is rebuilt from the molecules that pass curation (`distribution_success == True`) and describes the curated dataset.

How the parameters relate the two:

| Base vocabulary | `filter_unknown_tokens` | Result |
|---|---|---|
| built fresh from data | any | every token is known; curated ⊆ base |
| provided, `dynamically_update_vocab=True` | any | base grows to cover new tokens; curated ⊆ base |
| provided, `dynamically_update_vocab=False` | `True` | out-of-vocabulary molecules are removed; curated ⊆ base |
| provided, `dynamically_update_vocab=False` | `False` | out-of-vocabulary tokens are neither filtered nor added; the curated vocabulary may exceed the base vocabulary. This combination is flagged with a warning. |

## Curation results (`curation_results.json`)

`distributions` writes `curation_results.json` alongside the run output. Beyond the per-property thresholds and statistics, it records:

- **`summary`** — the curation funnel: the input count, the retained count, and the split of removals into a property stage and a token stage.
- **`token_curation`** — for each token filter (`unknown_tokens`, `rare_tokens`), the number of molecules removed and the identities of the tokens that drove those removals.

The console log mirrors the funnel (`Distribution curation: N → M retained …`) and lists the unknown and rare token identities.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `plot_distributions` | `True` | `bool` | Generate distribution plots for the computed properties |
| `perform_pca` | `False` | `bool` | Run PCA on the computed properties for visualization |

## Available Properties

`num_atoms`, `num_rings`, `size_largest_ring`, `num_tokens`, `tokens_atom_ratio`, `c_atom_ratio`, `longest_aliph_carbon`, `molecular_weight`, `logp`, `tpsa`, `num_rotatable_bonds`, `num_h_donors`, `num_h_acceptors`, `fsp3`, `num_aromatic_rings`, `num_stereocenters`, `num_heteroatoms`, `heteroatom_ratio`

## PropertyThreshold Configuration

`PropertyThreshold` supports three types of thresholds (priority: absolute > statistical > quantile). It is exported from the top level: `from molforge import PropertyThreshold`.

| Parameter | Type | Description |
|-----------|------|-------------|
| `min_value` | `float` \| `None` | Absolute lower bound (fixed minimum value); highest priority threshold type |
| `max_value` | `float` \| `None` | Absolute upper bound (fixed maximum value); highest priority threshold type |
| `statistical_lower` | `float` \| `None` | Lower bound as mean + n×std (e.g. `-2.0` for mean − 2σ) |
| `statistical_upper` | `float` \| `None` | Upper bound as mean + n×std (e.g. `2.0` for mean + 2σ) |
| `quantile_lower` | `float` \| `None` | Lower bound as a percentile (e.g. `0.025` for the 2.5th percentile) |
| `quantile_upper` | `float` \| `None` | Upper bound as a percentile (e.g. `0.975` for the 97.5th percentile) |

## Example

```python
from molforge import CurateDistributionParams, PropertyThreshold

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

## Related

- [Base Parameters](../configuration.md#base-parameters) inherited by all actors
- [Pipeline configuration](../configuration.md) — how `distributions_params` is set
- [Using this actor standalone](../standalone-actors.md) (token curation needs the `tokens` vocabulary via context)
- Previous step: [TokenizeData (`tokens`)](tokens.md)
- Next step: [GenerateConfs (`confs`)](confs.md)
- See the [example run](../example.md) for distribution plots produced by this actor
