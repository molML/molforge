"""
Integration smoke tests for the two-axis train/val/test splitter plugin.

The splitter depends on scaffold assignment, so each test first runs the
``scaffold`` actor to populate ``scaffold_smiles``, then runs ``Splitter``.

The splitter has two orthogonal axes: ``unit`` in {scaffold, molecule} and
``method`` in {isolation, random}.  The central guarantee under test for
``unit='scaffold'`` is *scaffold integrity*: no scaffold may appear in more
than one split, regardless of method.
"""

import json

import pytest
import pandas as pd

from molforge.actor_plugins.scaffold import ComputeScaffolds, ComputeScaffoldsParams
from molforge.actor_plugins.split import Splitter, SplitParams
from molforge.actors.protocol import ActorInput
from molforge.configuration.context import PipelineContext
from molforge.configuration.logger import PipelineLogger


@pytest.fixture
def split_data():
    """Curated-SMILES DataFrame with shared scaffolds, acyclic, and invalid rows."""
    return pd.DataFrame({
        'curated_smiles': [
            'c1ccccc1',            # benzene        \
            'Cc1ccccc1',           # toluene         >  share scaffold c1ccccc1
            'CCc1ccccc1',          # ethylbenzene   /
            'c1ccc2ccccc2c1',      # naphthalene     \  share naphthalene scaffold
            'Cc1ccc2ccccc2c1',     # methylnaphthalene/
            'C1CCCCC1',            # cyclohexane     \  share cyclohexane scaffold
            'CC1CCCCC1',           # methylcyclohexane/
            'c1ccncc1',            # pyridine
            'c1ccc(-c2ccccc2)cc1', # biphenyl
            'CCO',                 # ethanol  - acyclic -> train
            'CCCCC',               # pentane  - acyclic -> train
            'not_a_smiles',        # invalid  - failed scaffold -> train
        ],
        'molecule_id': [
            'benzene', 'toluene', 'ethylbenzene',
            'naphthalene', 'methylnaphthalene',
            'cyclohexane', 'methylcyclohexane',
            'pyridine', 'biphenyl',
            'ethanol', 'pentane', 'invalid',
        ],
    })


@pytest.fixture
def output_dir(tmp_path):
    out_dir = tmp_path / "test_output"
    out_dir.mkdir(exist_ok=True)
    return out_dir


@pytest.fixture
def logger(output_dir):
    return PipelineLogger(
        name='test_split_plugin',
        log_file=str(output_dir / 'test_split_plugin.log'),
        console_level='INFO',
    )


@pytest.fixture
def pipeline_context(output_dir):
    return PipelineContext(
        run_id='test_split_plugin',
        output_dir=str(output_dir),
    )


def _run_scaffold(df, logger, context):
    """Run the scaffold actor and return the DataFrame with scaffold_smiles."""
    params = ComputeScaffoldsParams(SMILES_column='curated_smiles')
    actor = ComputeScaffolds(params, logger=logger)
    out = actor(ActorInput(data=df, context=context))
    assert out.success
    return out.data


def _run_split(scaffolded, logger, context, **param_kwargs):
    """Run the splitter with the given parameter overrides; return the output."""
    # Silence the small-dataset advisory unless a test overrides it.
    param_kwargs.setdefault('low_n_threshold', 0)
    params = SplitParams(
        SMILES_column='curated_smiles',
        scaffold_column='scaffold_smiles',
        test_ratio=0.10,
        val_ratio=0.10,
        **param_kwargs,
    )
    splitter = Splitter(params, logger=logger)
    return splitter(ActorInput(data=scaffolded, context=context))


@pytest.mark.integration
class TestSplitter:
    """Smoke tests for the two-axis splitter plugin."""

    def test_split_labels_integrity_and_report(self, split_data, logger, pipeline_context):
        """Default scaffold/isolation: valid labels, scaffold integrity, report."""
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        output = _run_split(scaffolded, logger, pipeline_context)

        assert output.success
        df = output.data
        assert 'split' in df.columns

        # 1. Only valid split labels are produced.
        assert set(df['split']).issubset({'train', 'val', 'test'})

        # 2. SCAFFOLD INTEGRITY: no scaffold appears in more than one split.
        ring = df[df['scaffold_smiles'].fillna("") != ""]
        split_per_scaffold = ring.groupby('scaffold_smiles')['split'].nunique()
        assert (split_per_scaffold == 1).all(), (
            "scaffold leakage detected: a scaffold spans multiple splits"
        )

        # 3. Acyclic / failed-scaffold rows always land in train.
        acyclic_or_failed = df[df['scaffold_smiles'].fillna("") == ""]
        assert (acyclic_or_failed['split'] == 'train').all()
        # sanity: our fixture's acyclic/invalid rows are represented here
        assert {'ethanol', 'pentane', 'invalid'}.issubset(set(acyclic_or_failed['molecule_id']))

        # 4. The report card JSON is written and parseable.
        report_path = output.endpoint
        assert report_path is not None
        with open(report_path) as fh:
            report = json.load(fh)
        assert report['n_total'] == len(df)
        assert report['assignment_method'] == 'scaffold/isolation'
        # No overlap of test/val scaffolds with train in the report either.
        assert report["scaffold_overlap_test_train"] == 0
        assert report["scaffold_overlap_val_train"] == 0

    def test_report_figure_when_matplotlib_available(self, split_data, logger, pipeline_context):
        """The PNG report figure is written iff matplotlib is importable."""
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        output = _run_split(scaffolded, logger, pipeline_context)
        assert output.success

        try:
            import matplotlib  # noqa: F401
            has_mpl = True
        except ImportError:
            has_mpl = False

        png_path = str(output.endpoint).replace('.json', '.png')
        import os
        if has_mpl:
            assert os.path.exists(png_path), "expected report figure PNG to be written"
        # When matplotlib is missing the plugin degrades gracefully (no PNG);
        # nothing to assert beyond the successful run above.

    @pytest.mark.parametrize("unit", ["scaffold", "molecule"])
    @pytest.mark.parametrize("method", ["isolation", "random"])
    def test_all_combinations_valid(self, split_data, logger, pipeline_context, unit, method):
        """All four unit x method combinations produce valid labels and ratios."""
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        output = _run_split(scaffolded, logger, pipeline_context, unit=unit, method=method)
        assert output.success

        df = output.data
        n = len(df)
        assert set(df['split']).issubset({'train', 'val', 'test'})

        with open(output.endpoint) as fh:
            report = json.load(fh)
        assert report['assignment_method'] == f"{unit}/{method}"

        # Ratios are within rounding of the requested targets.  With a tiny
        # dataset each set gets at least one molecule (greedy fill floors at 1),
        # so allow a generous tolerance driven by the single-molecule quantum.
        tol = 1.0 / n + 1e-9
        assert abs(report['test_ratio'] - 0.10) <= 0.10 + tol
        assert abs(report['val_ratio'] - 0.10) <= 0.10 + tol
        assert report['n_train'] + report['n_val'] + report['n_test'] == n
        assert report['n_train'] > 0 and report['n_val'] > 0 and report['n_test'] > 0

    @pytest.mark.parametrize("method", ["isolation", "random"])
    def test_scaffold_unit_no_leakage(self, split_data, logger, pipeline_context, method):
        """unit='scaffold' has zero scaffold leakage under both methods."""
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        output = _run_split(scaffolded, logger, pipeline_context, unit='scaffold', method=method)
        assert output.success

        df = output.data
        ring = df[df['scaffold_smiles'].fillna("") != ""]
        split_per_scaffold = ring.groupby('scaffold_smiles')['split'].nunique()
        assert (split_per_scaffold == 1).all()

        with open(output.endpoint) as fh:
            report = json.load(fh)
        assert report['scaffold_overlap_test_train'] == 0
        assert report['scaffold_overlap_val_train'] == 0
        assert report['scaffold_overlap_test_val'] == 0

    def test_isolation_is_deterministic(self, split_data, logger, pipeline_context):
        """isolation ignores the seed and yields identical splits across runs."""
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        out_a = _run_split(scaffolded, logger, pipeline_context,
                           method='isolation', seed=1)
        out_b = _run_split(scaffolded, logger, pipeline_context,
                           method='isolation', seed=999)
        assert list(out_a.data['split']) == list(out_b.data['split'])

    def test_random_seed_reproducibility(self, split_data, logger, pipeline_context):
        """random is reproducible per seed and generally differs across seeds."""
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)

        same_1 = _run_split(scaffolded, logger, pipeline_context,
                            unit='molecule', method='random', seed=42)
        same_2 = _run_split(scaffolded, logger, pipeline_context,
                            unit='molecule', method='random', seed=42)
        assert list(same_1.data['split']) == list(same_2.data['split'])

        # A different seed generally produces a different permutation.  Scan a
        # handful of seeds to make the "generally different" claim robust.
        baseline = list(same_1.data['split'])
        differs = any(
            list(_run_split(scaffolded, logger, pipeline_context,
                            unit='molecule', method='random', seed=s).data['split'])
            != baseline
            for s in (1, 7, 13, 100, 2024)
        )
        assert differs, "no tested alternate seed changed the random split"

    def test_over_limit_raises(self, split_data, logger, pipeline_context):
        """method='isolation' over max_units_for_isolation raises ValueError.

        ``process`` raises so the error propagates to the framework, which halts
        the run; the actor never silently falls back to another ordering.
        """
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        params = SplitParams(
            SMILES_column='curated_smiles',
            scaffold_column='scaffold_smiles',
            method='isolation',
            max_units_for_isolation=1,
            low_n_threshold=0,
        )
        splitter = Splitter(params, logger=logger)
        splitter._context = pipeline_context
        with pytest.raises(ValueError, match="max_units_for_isolation"):
            splitter.process(scaffolded)

        # The framework wrapper surfaces the same failure as an unsuccessful
        # ActorOutput, terminating the pipeline rather than falling back.
        output = splitter(ActorInput(data=scaffolded, context=pipeline_context))
        assert not output.success
        assert output.metadata.get('exception_type') == 'ValueError'

    def test_low_n_threshold_warns(self, split_data, logger, pipeline_context, caplog):
        """A dataset below low_n_threshold emits an advisory WARNING but runs."""
        import logging
        scaffolded = _run_scaffold(split_data, logger, pipeline_context)
        with caplog.at_level(logging.WARNING):
            output = _run_split(
                scaffolded, logger, pipeline_context,
                low_n_threshold=1_000_000,   # force the advisory
            )
        assert output.success
        assert any("low_n_threshold" in rec.getMessage() for rec in caplog.records)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
