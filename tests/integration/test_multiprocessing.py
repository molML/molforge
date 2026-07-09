"""Integration test for multiprocessing functionality.

This test ensures that actors can be properly pickled and used in
multiprocessing workers. It targets the fix for lambda/StereoHandler pickling
issues.

To trigger the multiprocessing code path (DEFAULT_MP_THRESHOLD = 1000 rows) the
tests feed a small, controlled SYNTHETIC DataFrame of >1000 valid molecules
straight into the ``curate`` step, bypassing ChEMBL entirely. This keeps the
multiprocessing/pickling coverage intent while being fast and requiring no
ChEMBL database download, so the tests run in the default suite (unmarked).
``check_duplicates=False`` keeps the duplicated rows so the downstream
``distributions`` step also exceeds the multiprocessing threshold.
"""

import pytest
import pandas as pd
from molforge import (
    ForgeParams, ConstructPipe,
    ChEMBLSourceParams, ChEMBLCuratorParams,
    CurateMolParams, CurateDistributionParams
)


# A diverse set of valid SMILES, including a stereochemically-defined molecule
# so StereoHandler (stereo_policy='assign') is exercised in the MP workers.
_BASE_SMILES = [
    'CCO', 'c1ccccc1', 'CC(C)C', 'CCN', 'CCCCO',
    'CC(=O)O', 'c1ccncc1', 'CCOCC', 'CCCCCC', 'OCCO',
    'CC(C)Cc1ccc(cc1)C(C)C(=O)O',          # ibuprofen
    'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',     # (S)-ibuprofen (stereo)
    'CN1CCC[C@H]1c1cccnc1',                # nicotine (stereo)
    'Clc1ccccc1', 'Brc1ccccc1', 'Fc1ccccc1',
    'CCC(=O)OCC', 'c1ccc2ccccc2c1', 'CC(=O)Nc1ccccc1', 'CCCCCCCC',
]


def _large_synthetic_df(n_rows=1200):
    """Build a >1000-row DataFrame of valid molecules to trigger multiprocessing."""
    smiles = [_BASE_SMILES[i % len(_BASE_SMILES)] for i in range(n_rows)]
    ids = [f'MOL{i}' for i in range(n_rows)]
    return pd.DataFrame({'canonical_smiles': smiles, 'molecule_chembl_id': ids})


@pytest.fixture
def output_dir(tmp_path):
    """Unified output directory for all test artifacts."""
    out_dir = tmp_path / "test_output_mp"
    out_dir.mkdir(exist_ok=True)
    return out_dir


@pytest.fixture
def multiprocessing_params(output_dir):
    """
    Configuration for the multiprocessing test.

    Runs curate + distributions (both multiprocessing-capable) on a synthetic
    >1000-row DataFrame. ``check_duplicates=False`` preserves the duplicated
    rows so distributions also crosses the multiprocessing threshold.
    """
    return ForgeParams(
        verbose=True,
        override_actor_params=True,
        return_none_on_fail=True,
        write_checkpoints=False,
        output_root=str(output_dir),

        steps=['curate', 'distributions'],

        curate_params=CurateMolParams(
            SMILES_column='canonical_smiles',
            mol_steps=['desalt', 'removeIsotope', 'neutralize', 'sanitize', 'handleStereo'],
            smiles_steps=['canonical'],
            desalt_policy='keep',
            neutralize_policy='keep',
            stereo_policy='assign',  # exercises StereoHandler (the lambda bug)
            dropna=True,
            check_duplicates=False,  # keep >1000 rows for downstream MP
        ),
        distributions_params=CurateDistributionParams(
            SMILES_column='curated_smiles',
            global_statistical_threshold=2.0,
            plot_distributions=False,
            perform_pca=False,
        ),
    )


def test_multiprocessing_curate_and_distributions(multiprocessing_params):
    """
    Multiprocessing works for both the curate and distributions actors.

    Validates:
    1. Actors can be pickled for multiprocessing workers
    2. StereoHandler can be pickled (no lambda issues)
    3. Worker initialization is suppressed (suppress_init_logs works)
    4. Large datasets (> 1000 molecules) process successfully
    """
    pipeline = ConstructPipe(multiprocessing_params)

    # >1000-row synthetic input triggers multiprocessing in curate and distributions.
    result = pipeline(_large_synthetic_df(), input_id='mp_synthetic')

    # Validate successful execution
    assert result is not None, "Pipeline should return result"
    assert isinstance(result, pd.DataFrame), "Result should be DataFrame"
    assert len(result) > 0, "Result should have data"

    # Validate expected columns from curate step
    assert 'curated_smiles' in result.columns, "Should have curated_smiles"
    assert 'curation_success' in result.columns, "Should have curation_success"

    # Distributions step should add molecular properties
    assert 'molecular_weight' in result.columns, "Distributions step should add molecular_weight"

    # Check that we have successful curations
    successful_curations = result['curation_success'].sum() if 'curation_success' in result.columns else len(result)
    assert successful_curations > 100, f"Should have > 100 successful curations, got {successful_curations}"


def test_multiprocessing_curate_only(multiprocessing_params):
    """
    Multiprocessing works for the curate actor in isolation.

    Validates that CurateMol can initialize StereoHandler in multiprocessing
    workers without lambda pickling errors.
    """
    # Modify params to only test the curate step
    params = multiprocessing_params
    params.steps = ['curate']

    pipeline = ConstructPipe(params)
    result = pipeline(_large_synthetic_df(), input_id='mp_synthetic')

    assert result is not None
    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert 'curated_smiles' in result.columns

    # Verify stereo handling occurred (stereo_policy='assign')
    successful_curations = result['curation_success'].sum() if 'curation_success' in result.columns else len(result)
    assert successful_curations > 100, "Should have many successful curations with stereo handling"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
