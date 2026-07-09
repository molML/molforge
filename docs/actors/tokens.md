# TokenizeData (`tokens`)

[← Back to README](../../README.md) · [Documentation index](../../README.md#documentation)

Tokenizes SMILES strings into token sequences for machine-learning applications. Supports dynamic vocabulary updating, which is useful for building a single vocabulary across several forge inputs (e.g. an array of protein targets).

Typically consumes the output of the [`curate` actor](curate.md) and feeds the [`distributions` actor](distributions.md), which also performs token-based curation. See [Vocabulary handling](distributions.md#vocabulary-handling) for how the base vocabulary built here relates to the curated vocabulary.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `vocab_file` | `None` | `str` \| `None` | Path to a vocabulary JSON file: loaded if it exists, otherwise a new vocabulary is built from the data and saved to this path. When `None`, uses (or creates) `vocab.json` in the run directory |
| `SMILES_column` | `'curated_smiles'` | `str` | DataFrame column containing SMILES strings to tokenize |
| `dynamically_update_vocab` | `True` | `bool` | When a non-empty vocabulary is loaded from file, extend it with any new tokens found in the data and re-save it. Ignored when a fresh vocabulary is built (which always includes all tokens from the data) |

## Example

```python
from molforge import TokenizeDataParams

# Dynamic vocabulary
tokens_params = TokenizeDataParams(dynamically_update_vocab=True)

# Fixed vocabulary
tokens_params = TokenizeDataParams(vocab_file='vocab.json', dynamically_update_vocab=False)
```

## Related

- [Base Parameters](../configuration.md#base-parameters) inherited by all actors
- [Pipeline configuration](../configuration.md) — how `tokens_params` is set
- [Using this actor standalone](../standalone-actors.md)
- Previous step: [CurateMol (`curate`)](curate.md)
- Next step: [CurateDistribution (`distributions`)](distributions.md)
