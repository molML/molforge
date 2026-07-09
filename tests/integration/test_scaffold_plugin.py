"""
Integration smoke tests for the Bemis-Murcko scaffold plugin actor.

Mirrors the fixture/style of ``test_confs.py``: a small hardcoded SMILES
DataFrame, an actor constructed from a params object + logger, and execution
via the actor's ``__call__`` with an ``ActorInput`` carrying a
``PipelineContext``.
"""

import pytest
import pandas as pd

from molforge.actor_plugins.scaffold import ComputeScaffolds, ComputeScaffoldsParams
from molforge.actors.protocol import ActorInput
from molforge.configuration.context import PipelineContext
from molforge.configuration.logger import PipelineLogger


@pytest.fixture
def scaffold_data():
    """Small curated-SMILES DataFrame mixing ring and acyclic molecules."""
    return pd.DataFrame({
        'curated_smiles': [
            'CCO',                 # ethanol   - acyclic  -> scaffold ""
            'CCCC',                # butane    - acyclic  -> scaffold ""
            'c1ccccc1',            # benzene   - ring     -> non-empty scaffold
            'Cc1ccccc1',           # toluene   - ring     -> same scaffold as benzene
            'c1ccc2ccccc2c1',      # naphthalene - ring   -> distinct scaffold
        ],
        'molecule_id': ['ethanol', 'butane', 'benzene', 'toluene', 'naphthalene'],
    })


@pytest.fixture
def output_dir(tmp_path):
    out_dir = tmp_path / "test_output"
    out_dir.mkdir(exist_ok=True)
    return out_dir


@pytest.fixture
def logger(output_dir):
    return PipelineLogger(
        name='test_scaffold_plugin',
        log_file=str(output_dir / 'test_scaffold_plugin.log'),
        console_level='INFO',
    )


@pytest.fixture
def pipeline_context(output_dir):
    return PipelineContext(
        run_id='test_scaffold_plugin',
        output_dir=str(output_dir),
    )


@pytest.mark.integration
class TestComputeScaffolds:
    """Smoke tests for the scaffold computation plugin."""

    def test_output_columns_and_acyclic_vs_ring(self, scaffold_data, logger, pipeline_context):
        """All three output columns exist; acyclic -> '' (success), ring -> non-empty."""
        params = ComputeScaffoldsParams(
            SMILES_column='curated_smiles',
            include_generic=True,
        )
        actor = ComputeScaffolds(params, logger=logger)
        output = actor(ActorInput(data=scaffold_data, context=pipeline_context))

        assert output.success
        df = output.data

        for col in ('scaffold_smiles', 'scaffold_generic_smiles', 'scaffold_success'):
            assert col in df.columns, f"missing output column {col}"

        # Acyclic molecule: empty scaffold string, but success is True.
        eth = df[df['molecule_id'] == 'ethanol'].iloc[0]
        assert eth['scaffold_smiles'] == ""
        assert bool(eth['scaffold_success']) is True

        # Ring molecule: non-empty scaffold, success True.
        benz = df[df['molecule_id'] == 'benzene'].iloc[0]
        assert isinstance(benz['scaffold_smiles'], str) and benz['scaffold_smiles'] != ""
        assert bool(benz['scaffold_success']) is True

        # Generic scaffold populated when include_generic=True.
        assert benz['scaffold_generic_smiles'] is not None

        # Molecules sharing a ring system share the same Murcko scaffold.
        tol = df[df['molecule_id'] == 'toluene'].iloc[0]
        assert tol['scaffold_smiles'] == benz['scaffold_smiles']

    def test_include_generic_false_keeps_column_as_none(self, scaffold_data, logger, pipeline_context):
        """With include_generic=False the generic column exists but holds None."""
        params = ComputeScaffoldsParams(
            SMILES_column='curated_smiles',
            include_generic=False,
        )
        actor = ComputeScaffolds(params, logger=logger)
        output = actor(ActorInput(data=scaffold_data, context=pipeline_context))

        df = output.data
        assert 'scaffold_generic_smiles' in df.columns
        assert df['scaffold_generic_smiles'].isna().all()

    def test_acyclic_policy_remove_drops_acyclic_rows(self, scaffold_data, logger, pipeline_context):
        """acyclic_policy='remove' drops rows whose scaffold is empty."""
        params = ComputeScaffoldsParams(
            SMILES_column='curated_smiles',
            acyclic_policy='remove',
        )
        actor = ComputeScaffolds(params, logger=logger)
        output = actor(ActorInput(data=scaffold_data, context=pipeline_context))

        df = output.data
        # No empty scaffolds remain; acyclic molecules were dropped.
        assert (df['scaffold_smiles'] != "").all()
        assert set(df['molecule_id']) == {'benzene', 'toluene', 'naphthalene'}


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
