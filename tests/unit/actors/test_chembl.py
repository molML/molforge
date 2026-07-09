# Baseline test for ChEMBLCurator

import pytest
import pandas as pd
from molforge.actors.chembl import ChEMBLCurator
from molforge.actors.params.chembl import ChEMBLCuratorParams
from molforge.configuration.steps import Steps

def test_chembl_curator_init():
    """Test ChEMBLCurator initializes."""
    params = ChEMBLCuratorParams()
    actor = ChEMBLCurator(params, logger=None)
    assert Steps.CHEMBL.upper() in actor._nametag

def test_chembl_curator_basic():
    """Test ChEMBLCurator processes DataFrame."""
    from molforge.actors.protocol import ActorInput
    from molforge.configuration.context import PipelineContext

    # ChEMBL-like data with all required columns
    df = pd.DataFrame({
        'canonical_smiles': ['CCO', 'c1ccccc1'],
        'standard_value': [100, 200],
        'standard_units': ['nM', 'nM'],
        'standard_type': ['IC50', 'IC50'],
        'standard_relation': ['=', '='],
        'standard_flag': [True, True],
        'potential_duplicate': [False, False],
        'pchembl_value': [7.0, 6.7],
        'data_validity_comment': [None, None],
        'document_year': [2020, 2021],
        'document_chembl_id': ['CHEMBL1', 'CHEMBL2'],
        'assay_type': ['B', 'B'],
        'assay_description': ['test assay', 'test assay'],
        'target_organism': ['Homo sapiens', 'Homo sapiens'],
        'bao_format': ['BAO_0000357', 'BAO_0000357'],
        'molecule_chembl_id': ['CHEMBL100', 'CHEMBL200']
    })

    params = ChEMBLCuratorParams(verbose=False)
    actor = ChEMBLCurator(params, logger=None)

    # Create context for actor
    context = PipelineContext(run_id='test', output_dir='/tmp/test_chembl')
    actor_input = ActorInput(data=df, context=context)

    result = actor(actor_input)

    assert result.success
    assert isinstance(result.data, pd.DataFrame)
    assert len(result.data) > 0
