# Simple baseline test for CurateMol

import pytest
import pandas as pd
from molforge.actors.curate import CurateMol
from molforge.actors.params.curate import CurateMolParams
from molforge.configuration.steps import Steps

def test_curate_mol_init():
    """Test CurateMol initializes."""
    params = CurateMolParams()
    actor = CurateMol(params, logger=None)
    assert Steps.CURATE.upper() in actor._nametag

def test_curate_mol_basic():
    """Test CurateMol processes DataFrame."""
    from molforge.actors.protocol import ActorInput
    from molforge.configuration.context import PipelineContext

    # Create test data
    df = pd.DataFrame({
        'smiles': ['CCO', 'c1ccccc1', 'CC(C)C']
    })

    params = CurateMolParams(verbose=False, SMILES_column='smiles')
    actor = CurateMol(params, logger=None)

    # Create context for actor
    context = PipelineContext(run_id='test', output_dir='/tmp/test_curate')
    actor_input = ActorInput(data=df, context=context)

    result = actor(actor_input)

    assert result.success
    assert isinstance(result.data, pd.DataFrame)
    assert 'curated_smiles' in result.data.columns
    assert len(result.data) > 0
