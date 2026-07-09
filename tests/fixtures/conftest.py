# pytest configuration and fixtures

import pytest
import pandas as pd

@pytest.fixture
def sample_smiles_df():
    """Small sample DataFrame with SMILES."""
    return pd.DataFrame({
        'canonical_smiles': ['CCO', 'c1ccccc1', 'CC(C)C'],
        'molecule_chembl_id': ['CHEMBL545', 'CHEMBL277', 'CHEMBL123']
    })

@pytest.fixture
def sample_chembl_id():
    """Sample ChEMBL target ID."""
    return "CHEMBL123"
