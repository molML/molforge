"""
Unit tests for CurateDistribution token filtering (filter_unknown_tokens).

These target ``_filter_by_unknown_tokens`` directly with a constructed
``Vocabulary`` rather than wiring a full TokenizeData actor + pipeline context,
which keeps the test fast while still exercising the real filtering logic.
"""

import pytest
import pandas as pd

from molforge.actors.distributions import CurateDistribution
from molforge.actors.params.distributions import CurateDistributionParams, PropertyThreshold
from molforge.configuration.steps import Steps
from molforge.utils.actortools.vocabulary import Vocabulary


@pytest.fixture
def token_df():
    """DataFrame with a tokens column; one molecule contains an 'N' token."""
    return pd.DataFrame({
        'curated_smiles': ['CCO', 'c1ccccc1', 'CCN'],
        'tokens': [
            ['C', 'C', 'O'],        # all-standard tokens
            ['c', '1', 'c', 'c', 'c', 'c', 'c', '1'],
            ['C', 'C', 'N'],        # contains 'N'
        ],
        'molecule_id': ['ethanol', 'benzene', 'ethylamine'],
    })


@pytest.fixture
def actor():
    """A CurateDistribution actor (default params enable token curation)."""
    params = CurateDistributionParams(
        SMILES_column='curated_smiles',
        tokens_column='tokens',
        curate_tokens=True,
        filter_unknown_tokens=True,
    )
    # logger=None -> actor prints (verbose); acceptable for a unit test.
    return CurateDistribution(params, logger=None)


class TestFilterUnknownTokens:
    """Tests for CurateDistribution._filter_by_unknown_tokens."""

    def test_fixed_vocab_missing_token_removes_molecules(self, actor, token_df):
        """A fixed vocab missing 'N' drops molecules whose tokens include 'N'."""
        vocab = Vocabulary()
        # Deliberately omit 'N' so 'CCN' contains an out-of-vocabulary token.
        vocab.update(['C', 'O', 'c', '1'])

        result = actor._filter_by_unknown_tokens(token_df, vocab)

        assert len(result) == 2
        assert set(result['molecule_id']) == {'ethanol', 'benzene'}
        assert 'ethylamine' not in set(result['molecule_id'])

    def test_all_tokens_known_removes_nothing(self, actor, token_df):
        """When every token is in the vocabulary, no molecule is removed."""
        vocab = Vocabulary()
        vocab.update(['C', 'O', 'c', '1', 'N'])  # includes every token present

        result = actor._filter_by_unknown_tokens(token_df, vocab)

        assert len(result) == len(token_df)
        assert set(result['molecule_id']) == {'ethanol', 'benzene', 'ethylamine'}


class _StubTokenizeActor:
    """Minimal stand-in for TokenizeData that only exposes a vocabulary."""

    def __init__(self, vocab):
        self._vocab = vocab

    def get_vocabulary(self):
        return self._vocab


class _StubContext:
    """Minimal pipeline context exposing an upstream TokenizeData actor + output dir."""

    def __init__(self, output_dir, td_actor):
        self.output_dir = output_dir
        self._td = td_actor

    def get_actor(self, attr_name):
        return self._td if attr_name == Steps.TOKENS else None

    def get_actor_result(self, attr_name):
        return None

    def get(self, key, default=None):
        return default

    def set(self, key, value):
        pass


@pytest.fixture
def flag_df():
    """Six molecules exercising every failure mode (property low/high, OOV, rare)."""
    return pd.DataFrame({
        'molecule_id':    ['pass', 'atoms_high', 'atoms_low', 'logp_high', 'oov', 'rare'],
        'curated_smiles': ['CCO', 'CCO', 'CCO', 'CCO', 'CCO', 'CCO'],
        'tokens': [
            ['C', 'C', 'O'],   # pass
            ['C', 'O'],        # fails num_atoms:high
            ['C', 'O'],        # fails num_atoms:low
            ['C', 'O'],        # fails logp:high
            ['C', 'O', 'Q'],   # 'Q' absent from vocab (unknown) and appears in 1/6 (rare)
            ['C', 'O', 'Z'],   # 'Z' appears in 1/6 (rare)
        ],
        'num_atoms': [10.0, 40.0, 3.0, 10.0, 10.0, 10.0],
        'logp':      [2.0, 2.0, 2.0, 5.0, 2.0, 2.0],
    })


def _make_flag_actor(dropna, tmp_path):
    params = CurateDistributionParams(
        SMILES_column='curated_smiles',
        tokens_column='tokens',
        properties=['num_atoms', 'logp'],
        thresholds={
            'num_atoms': PropertyThreshold(min_value=5, max_value=30),
            'logp': PropertyThreshold(max_value=4.0),
        },
        compute_properties=False,
        curate_tokens=True,
        filter_unknown_tokens=True,
        token_frequency_threshold=60.0,  # 'Q' and 'Z' each appear in 1/6 molecules -> rare
        plot_distributions=False,
        perform_pca=False,
        dropna=dropna,
    )
    actor = CurateDistribution(params, logger=None)
    vocab = Vocabulary()
    vocab.update(['C', 'O', 'Z'])  # 'Q' deliberately absent
    actor._context = _StubContext(str(tmp_path), _StubTokenizeActor(vocab))
    return actor


class TestDropnaFlagging:
    """Tests for the dropna / distribution_success / distribution_failures wiring."""

    def test_dropna_true_drops_failures_and_flags_survivors(self, flag_df, tmp_path):
        actor = _make_flag_actor(dropna=True, tmp_path=tmp_path)
        result = actor.process(flag_df)

        # New columns present.
        assert 'distribution_success' in result.columns
        assert 'distribution_failures' in result.columns

        # Only the fully-passing molecule survives.
        assert set(result['molecule_id']) == {'pass'}
        # Survivors are all flagged success with empty failure strings.
        assert result['distribution_success'].all()
        assert (result['distribution_failures'] == '').all()

    def test_dropna_false_retains_all_with_correct_flags(self, flag_df, tmp_path):
        actor = _make_flag_actor(dropna=False, tmp_path=tmp_path)
        result = actor.process(flag_df)

        # No rows dropped.
        assert len(result) == len(flag_df)

        success = dict(zip(result['molecule_id'], result['distribution_success']))
        failures = dict(zip(result['molecule_id'], result['distribution_failures']))

        # Exactly the passing molecule succeeds.
        assert success['pass'] is True or success['pass'] == True  # noqa: E712
        assert failures['pass'] == ''
        for mid in ['atoms_high', 'atoms_low', 'logp_high', 'oov', 'rare']:
            assert not success[mid], f"{mid} should fail"

        # Reasons are accurate per failure mode.
        # Every violated filter is listed per molecule (evaluated on all rows).
        assert failures['atoms_high'] == 'num_atoms:high'
        assert failures['atoms_low'] == 'num_atoms:low'
        assert failures['logp_high'] == 'logp:high'
        assert failures['oov'] == 'tokens:unknown; tokens:rare'
        assert failures['rare'] == 'tokens:rare'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
