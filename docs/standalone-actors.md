# Using Actors Standalone

[← Back to README](../README.md) · [Documentation index](../README.md#documentation)

Every actor can be used on its own, outside a `MolForge` pipeline. Each actor is constructed from its `*Params` object and a logger, then run on a DataFrame. This is useful for testing, notebooks, or reusing a single step (e.g. molecule curation) in another workflow.

## Two Execution Patterns

Actors inherit from `BaseActor`, constructed as `Actor(params, logger=None)` (`logger=None` prints to the terminal via the actor's `verbose` flag). There are two ways to run one:

1. **Direct** — `actor.process(df)` runs the core logic and returns the processed DataFrame. It does **not** run input validation, produce an `ActorOutput`, or provide a pipeline context.
2. **Full wrapper** — `actor(ActorInput(data=df, context=context))` is the officially supported execution path (the same one the pipeline uses). It validates input, logs, supplies a `PipelineContext` (output directory, cross-actor access), and returns an `ActorOutput` with `.data`, `.success`, and `.metadata`.

### Minimal direct example

```python
import pandas as pd
from molforge import CurateMolParams
from molforge.actors import CurateMol

df = pd.DataFrame({'canonical_smiles': ['CC(=O)Oc1ccccc1C(=O)O', 'c1ccccc1']})

params = CurateMolParams(mol_steps=['desalt', 'neutralize', 'sanitize'], dropna=True)
actor = CurateMol(params, logger=None)          # logger=None → print to terminal

curated = actor.process(df)                       # returns a DataFrame
```

### Full wrapper example (with a context)

```python
import pandas as pd
from molforge import GenerateConfsParams
from molforge.actors import GenerateConfs
from molforge.actors.protocol import ActorInput
from molforge.configuration.context import PipelineContext

df = pd.DataFrame({
    'curated_smiles': ['CCO', 'c1ccccc1'],
    'molecule_chembl_id': ['m1', 'm2'],
})

params = GenerateConfsParams(backend='rdkit', max_confs=10, dropna=False)
actor = GenerateConfs(params, logger=None)

context = PipelineContext(run_id='standalone', output_dir='./standalone_out')
output = actor(ActorInput(data=df, context=context))   # returns an ActorOutput

result_df = output.data
print(output.success, output.metadata)
```

## Which Actors Need a Context

Some actors only transform the DataFrame and run fine with `process(df)` directly. Others read cross-actor state through the `PipelineContext` and should be run with the full wrapper (or inside a pipeline).

| Actor | Direct `process(df)` | Notes |
|-------|----------------------|-------|
| [`chembl`](actors/chembl.md) | Yes | Pure DataFrame curation; no context used |
| [`curate`](actors/curate.md) | Yes | Standalone it falls back to internal duplicate handling (`duplicates_policy`); inside a pipeline it reads the `chembl` actor via context for aggregation |
| [`tokens`](actors/tokens.md) | Yes | Without a context there is no run directory to write `vocab.json`; set `vocab_file` to persist the vocabulary |
| `properties` plugin | Yes | Pure descriptor computation |
| [`scaffold` plugin](plugins.md#scaffold--bemis-murcko-scaffolds-scaffoldpy) | Yes | Pure scaffold computation |
| [`source`](actors/source.md) | No | Fetches data rather than transforming a DataFrame; requires `input_id` from the context (`required_context_keys`). Use the wrapper with `PipelineContext(..., input_id=...)` |
| [`distributions`](actors/distributions.md) | Partial | Property-distribution filtering runs on the DataFrame, but token curation (`curate_tokens` / `filter_unknown_tokens`) reads the [`tokens` vocabulary via context](actors/distributions.md#vocabulary-handling); run it through the pipeline (or a context with the `tokens` actor registered) if you rely on it |
| [`confs`](actors/confs.md) | Partial | Works without a context but writes conformer/report files to default paths; pass a `PipelineContext` with `output_dir` to control the output location |
| [`split` plugin](plugins.md#split--trainvaltest-split-splitpy) | Partial | Reads the `chembl` result via context to auto-resolve the activity column and writes its report to the context output dir; set `activity_column` explicitly to run standalone |

For anything relying on cross-actor context, the simplest correct option is to run the full pipeline via the [Python API](python-api.md), which wires the context for you.

## Related

- [Python API](python-api.md) — running the full pipeline
- [Plugin Development](plugins.md) — plugins are actors too, and follow the same `process(df)` / `ActorInput` contract
- [Pipeline Configuration](configuration.md) — the `*Params` objects each actor takes
