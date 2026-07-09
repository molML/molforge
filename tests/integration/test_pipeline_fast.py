# Fast integration test for pipeline (source -> chembl -> curate -> tokens)
#
# This runs the REAL ChEMBL SQL pipeline exactly ONCE with a LIMITED query
# (search_all=False, n=100) so it stays fast, and is marked @pytest.mark.database
# so the default `-m "not database"` suite skips it (never downloads the DB).

from pathlib import Path

import pytest
import pandas as pd
from molforge import (
    ForgeParams, ConstructPipe,
    ChEMBLSourceParams, ChEMBLCuratorParams,
    CurateMolParams, TokenizeDataParams
)

# Location of the local ChEMBL 36 SQL database (symlinked into the worktree).
REPO_ROOT = Path(__file__).resolve().parents[2]
CHEMBL_DB = REPO_ROOT / "data" / "chembl" / "36.db"


@pytest.fixture
def output_dir(tmp_path):
    """Unified output directory for all test artifacts."""
    out_dir = tmp_path / "test_output"
    out_dir.mkdir(exist_ok=True)
    return out_dir


@pytest.fixture
def fast_params(output_dir):
    """Fast test configuration (minimal steps, limited query) with unified output dir."""
    if not CHEMBL_DB.exists():
        pytest.skip(f"ChEMBL 36 database not found at {CHEMBL_DB}; skipping real-DB pipeline test.")

    return ForgeParams(
        verbose=True,
        override_actor_params=True,
        return_none_on_fail=True,
        write_checkpoints=False,
        output_root=str(output_dir),

        steps=['source', 'chembl', 'curate', 'tokens'],

        curate_params=CurateMolParams(
            SMILES_column='canonical_smiles',
            mol_steps=['desalt', 'removeIsotope', 'neutralize', 'sanitize', 'handleStereo'],
            smiles_steps=['canonical'],
            desalt_policy='keep',
            neutralize_policy='keep',
            stereo_policy='assign',
            dropna=True,
        ),
        tokens_params=TokenizeDataParams(),
        source_params=ChEMBLSourceParams(
            backend='sql',
            version=36,
            db_path=str(CHEMBL_DB),
            auto_download=False,   # never download in tests
            search_all=False,      # LIMITED query ...
            n=100,                 # ... to 100 activities
        ),
    )


@pytest.mark.database
def test_pipeline_fast_execution(fast_params):
    """Test full pipeline executes without errors (single run, limited query)."""
    pipeline = ConstructPipe(fast_params)

    # Use small target for fast test
    result = pipeline('CHEMBL234')

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert 'tokens' in result.columns          # From tokens step
    assert 'curated_smiles' in result.columns  # From curate step
    assert len(result) > 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
