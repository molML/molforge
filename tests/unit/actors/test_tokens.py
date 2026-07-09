# Baseline test for TokenizeData

import pytest
import pandas as pd
from molforge.actors.tokens import TokenizeData
from molforge.actors.params.tokens import TokenizeDataParams
from molforge.configuration.steps import Steps

def test_tokenize_init():
    """Test TokenizeData initializes."""
    params = TokenizeDataParams()
    actor = TokenizeData(params, logger=None)
    assert Steps.TOKENS.upper() in actor._nametag

def test_tokenize_basic():
    """Test TokenizeData processes DataFrame."""
    from molforge.actors.protocol import ActorInput
    from molforge.configuration.context import PipelineContext
    import tempfile

    df = pd.DataFrame({
        'smiles': ['CCO', 'c1ccccc1', 'CC(C)C']
    })

    with tempfile.TemporaryDirectory() as tmpdir:
        params = TokenizeDataParams(verbose=False, vocab_file=None, SMILES_column='smiles')
        actor = TokenizeData(params, logger=None)

        # Create context for actor
        context = PipelineContext(run_id='test', output_dir=tmpdir)
        actor_input = ActorInput(data=df, context=context)

        result = actor(actor_input)

        assert result.success
        assert isinstance(result.data, pd.DataFrame)
        assert 'tokens' in result.data.columns
        assert 'seqlen' in result.data.columns
        assert len(result.data) == len(df)
