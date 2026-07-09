# Pipeline Configuration

[← Back to README](../README.md) · [Documentation index](../README.md#documentation)

The `ForgeParams` class (a user-friendly alias for `PipeParams`) configures the overall pipeline behavior and individual actor parameters.

## Base Parameters

All actor parameter classes inherit from `BaseParams`, which provides common functionality.

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `verbose` | `True` | `bool` | Emit informational (verbose) logging for this component |

## Pipeline Parameters

| Parameter | Default | Type | Description |
|-----------|---------|------|-------------|
| `steps` | `['source', 'chembl', 'curate', 'tokens', 'distributions']` | `List[str]` | Ordered list of actor step names to execute (defaults to the standard 5-step pipeline when `None`) |
| `return_none_on_fail` | `False` | `bool` | On step failure, return `None` in place of the step output; with `False`, the partial step output is returned. Step errors are caught internally and surfaced through the return value |
| `output_root` | `"MOLFORGE_OUT"` | `str` | Root directory for pipeline outputs |
| `write_output` | `True` | `bool` | Write the final pipeline output DataFrame to disk (per-step checkpoints are controlled separately by `write_checkpoints`) |
| `write_checkpoints` | `False` | `bool` | Write intermediate checkpoints between steps |
| `console_log_level` | `'INFO'` | `str` | Logging level for console output |
| `file_log_level` | `'DEBUG'` | `str` | Logging level for the log file |
| `override_actor_params` | `True` | `bool` | When `True`, copy each pipeline-level field value onto every actor's params (and any nested param dataclass) for any field whose name matches — e.g. the shared `verbose` flag. Matching top-level fields are always overwritten, regardless of the actor's own value |

## Setting Actor Parameters

Actor parameters are set via the `{step}_params` attributes: `source_params`, `chembl_params`, `curate_params`, `tokens_params`, `distributions_params`, `confs_params`. Plugin actors are configured via `plugin_params` (a dict keyed by step name). Each may be given as a params object or a plain dict.

```python
from molforge import ForgeParams, ChEMBLSourceParams

ForgeParams(
    steps=['source', 'chembl'],
    source_params=ChEMBLSourceParams(backend='sql'),
    chembl_params={'standard_type': 'IC50', 'standard_units': 'nM'}  # dict also accepted
)
```

Each core actor documents its own parameter tables and examples:

- [ChEMBLSource (`source`)](actors/source.md)
- [ChEMBLCurator (`chembl`)](actors/chembl.md)
- [CurateMol (`curate`)](actors/curate.md)
- [TokenizeData (`tokens`)](actors/tokens.md)
- [CurateDistribution (`distributions`)](actors/distributions.md)
- [GenerateConfs (`confs`)](actors/confs.md)

Plugin actors and their `plugin_params` configuration are covered in [Plugin Development](plugins.md). Config files (YAML/JSON) are covered in the [Python API](python-api.md) and [CLI](cli.md) pages.
