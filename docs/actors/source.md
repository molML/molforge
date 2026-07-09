# ChEMBLSource (`source`)

[← Back to README](../../README.md) · [Documentation index](../../README.md#documentation)

Retrieves bioactivity data from ChEMBL. Supports SQL (local database) and API (web service) backends. The SQL backend is significantly faster than the API once downloaded (<1s vs. >30s per protein, depending on the target).

Its output typically feeds the [`chembl` curator](chembl.md).

## Shared Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `backend` | `'sql'` | `'sql'` \| `'api'` | ChEMBL data source backend |
| `n` | `1000` | `int` | Maximum number of activity entries to fetch, applied only when `search_all` is `False` (ignored when `search_all=True`; in the API backend it also bounds the per-request page size) |
| `search_all` | `True` | `bool` | Whether to fetch all entries instead of limiting to `n` |

## SQL Backend Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `db_path` | `None` | `str` \| `None` | Path to the local ChEMBL SQLite database (SQL backend) |
| `version` | `'latest'` | `str` \| `int` | ChEMBL version to use (SQL backend) |
| `auto_download` | `True` | `bool` | If `True`, downloads the database if not found (SQL backend) |
| `download_dir` | `"./data/chembl"` | `str` | Directory to store the downloaded database (SQL backend) |

The API backend currently uses only the shared parameters.

## Example

```python
from molforge import ChEMBLSourceParams

# SQL backend
source_params = ChEMBLSourceParams(backend='sql', version=36, auto_download=True)

# API backend
source_params = ChEMBLSourceParams(backend='api', n=5000, search_all=False)
```

## Related

- [Base Parameters](../configuration.md#base-parameters) inherited by all actors
- [Pipeline configuration](../configuration.md) — how `source_params` is set
- [Using this actor standalone](../standalone-actors.md) (needs a `PipelineContext` with `input_id`)
- Next step: [ChEMBLCurator (`chembl`)](chembl.md)
