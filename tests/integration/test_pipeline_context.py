"""Integration tests for pipeline context and actor communication.

Design notes
------------
The full real-ChEMBL pipeline (source -> chembl -> curate -> tokens -> confs)
is expensive: it needs the ChEMBL SQL database and RDKit conformer generation.
To keep the suite fast and focused we run that full pipeline EXACTLY ONCE, via a
module-scoped fixture (``real_pipeline_run``) using a LIMITED query
(``search_all=False, n=100``) and a small ``max_confs``. Every real-DB assertion
test reads from that single shared run instead of re-executing the pipeline.

Those real-DB tests are marked ``@pytest.mark.database`` so the default
``-m "not database"`` suite skips them and never triggers a multi-GB download.
The fixture additionally skips gracefully if the ChEMBL DB is absent, and sets
``auto_download=False`` so it can never download.

Pure context / actor-communication mechanics that do NOT need ChEMBL are
exercised as fast, unmarked tests using a tiny synthetic DataFrame that starts
the pipeline at the ``curate`` step (bypassing source/chembl entirely).
"""

from types import SimpleNamespace
from pathlib import Path

import pytest
import pandas as pd

from molforge import (
    ForgeParams, ConstructPipe,
    ChEMBLSourceParams, ChEMBLCuratorParams,
    CurateMolParams, TokenizeDataParams,
    GenerateConfsParams
)
from molforge.actors.protocol import ActorOutput
from molforge.configuration.context import PipelineContext
from molforge.configuration.steps import Steps


# Location of the local ChEMBL 36 SQL database (symlinked into the worktree).
REPO_ROOT = Path(__file__).resolve().parents[2]
CHEMBL_DB = REPO_ROOT / "data" / "chembl" / "36.db"

# Small max_confs keeps RDKit embedding fast while still exercising the
# multi-conformer code path (n_conformers > 0 for valid molecules).
FAST_MAX_CONFS = 3


# --------------------------------------------------------------------------- #
# Shared curate/tokens/confs params (reused by real-DB and synthetic tests)
# --------------------------------------------------------------------------- #

def _curate_params(smiles_column='canonical_smiles'):
    return CurateMolParams(
        SMILES_column=smiles_column,
        mol_steps=['desalt', 'removeIsotope', 'neutralize', 'sanitize', 'handleStereo'],
        smiles_steps=['canonical'],
        desalt_policy='keep',
        neutralize_policy='keep',
        stereo_policy='assign',
        dropna=True,
    )


def _tokens_params():
    return TokenizeDataParams(
        vocab_file=None,  # Will create dynamically
        SMILES_column='curated_smiles',
        dynamically_update_vocab=True,
    )


def _confs_params():
    return GenerateConfsParams(
        backend='rdkit',
        max_confs=FAST_MAX_CONFS,
        rms_threshold=0.5,
        random_seed=42,
        use_uff=True,
        max_iterations=200,
    )


# --------------------------------------------------------------------------- #
# Real-ChEMBL pipeline: run ONCE, module scoped, limited query
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def real_pipeline_run(tmp_path_factory):
    """Run the full real-ChEMBL pipeline a SINGLE time with a limited query.

    Returns a namespace bundling the executed pipeline, its result DataFrame,
    context, run_id and output_dir so individual tests can assert against one
    shared run rather than each re-running the (expensive) pipeline.

    Skips if the ChEMBL SQL DB is absent so a fresh checkout without the DB does
    not attempt a download.
    """
    if not CHEMBL_DB.exists():
        pytest.skip(f"ChEMBL 36 database not found at {CHEMBL_DB}; skipping real-DB pipeline test.")

    output_dir = tmp_path_factory.mktemp("pipeline_context_real")

    params = ForgeParams(
        verbose=True,
        override_actor_params=True,
        return_none_on_fail=False,   # we want to see errors
        write_checkpoints=True,      # exercise checkpoint writing
        output_root=str(output_dir),

        steps=['source', 'chembl', 'curate', 'tokens', 'confs'],

        source_params=ChEMBLSourceParams(
            backend='sql',
            version=36,
            db_path=str(CHEMBL_DB),
            auto_download=False,     # never download in tests
            search_all=False,        # LIMITED query ...
            n=100,                   # ... to 100 activities
        ),
        chembl_params=ChEMBLCuratorParams(
            standard_type='IC50',
            standard_units='nM',
            standard_relation='=',
            assay_type='B',
            target_organism='Homo sapiens',
            assay_format='protein',
            std_threshold=0.5,
            range_threshold=0.5,
        ),
        curate_params=_curate_params(),
        tokens_params=_tokens_params(),
        confs_params=_confs_params(),
    )

    pipeline = ConstructPipe(params)
    result = pipeline('CHEMBL234')

    return SimpleNamespace(
        pipeline=pipeline,
        result=result,
        context=pipeline.context,
        run_id=pipeline.run_id,
        output_dir=Path(pipeline.output_dir),
    )


@pytest.mark.database
def test_pipeline_context_initialization(real_pipeline_run):
    """Pipeline properly initializes context and registers all actors."""
    context = real_pipeline_run.context

    assert isinstance(context, PipelineContext)
    assert context.get('input_id') == 'CHEMBL234'
    assert context.get('input_type') == 'ChEMBL ID'

    # Verify actors were registered
    assert context.get_actor(Steps.SOURCE) is not None
    assert context.get_actor(Steps.CHEMBL) is not None
    assert context.get_actor(Steps.CURATE) is not None
    assert context.get_actor(Steps.TOKENS) is not None
    assert context.get_actor(Steps.CONFS) is not None


@pytest.mark.database
def test_actor_communication(real_pipeline_run):
    """Actors communicate via context; each records success + key metadata."""
    context = real_pipeline_run.context

    # ChEMBL data source success tracking
    cs_result = context.get_actor_result(Steps.SOURCE)
    assert isinstance(cs_result, ActorOutput)
    assert cs_result.success is True
    assert 'n_activities' in cs_result.metadata
    assert cs_result.metadata['n_activities'] > 0

    # ChEMBL curator success tracking
    cc_result = context.get_actor_result(Steps.CHEMBL)
    assert isinstance(cc_result, ActorOutput)
    assert cc_result.success is True
    assert 'n_curated' in cc_result.metadata
    assert 'output_label' in cc_result.metadata
    assert 'standard_type' in cc_result.metadata
    assert cc_result.metadata['n_curated'] > 0

    # CurateMol success tracking
    cm_result = context.get_actor_result(Steps.CURATE)
    assert isinstance(cm_result, ActorOutput)
    assert cm_result.success is True
    assert 'n_curated' in cm_result.metadata or 'n_processed' in cm_result.metadata

    # TokenizeData success tracking and endpoint
    td_result = context.get_actor_result(Steps.TOKENS)
    assert isinstance(td_result, ActorOutput)
    assert td_result.success is True
    assert 'vocab_size' in td_result.metadata
    assert 'n_tokenized' in td_result.metadata
    assert td_result.endpoint is not None
    assert Path(td_result.endpoint).exists()

    # Conformer generation success tracking
    gc_result = context.get_actor_result(Steps.CONFS)
    assert isinstance(gc_result, ActorOutput)
    assert gc_result.success is True
    assert gc_result.metadata['backend'] == 'rdkit'
    assert 'molecules_processed' in gc_result.metadata
    assert 'success_count' in gc_result.metadata
    assert 'total_conformers' in gc_result.metadata
    assert 'avg_conformers_per_molecule' in gc_result.metadata
    assert gc_result.metadata['success_count'] > 0
    assert gc_result.metadata['success_count'] <= gc_result.metadata['molecules_processed']

    # Verify conformers were generated
    gc_actor = context.get_actor(Steps.CONFS)
    mols = list(gc_actor.extract_molecules())
    assert len(mols) > 0
    assert all(mol.GetNumConformers() > 0 for mol in mols)

    # Verify sequential processing
    all_results = [cs_result, cc_result, cm_result, td_result, gc_result]
    assert all(r.success for r in all_results), "All actors should succeed"


@pytest.mark.database
def test_pipeline_checkpoints(real_pipeline_run):
    """Pipeline writes one checkpoint per step with expected content."""
    run_id = real_pipeline_run.run_id
    output_dir = real_pipeline_run.output_dir

    # ChEMBL source checkpoint verification
    cs_file = output_dir / f"{run_id}_source.csv"
    assert cs_file.exists()
    cs_df = pd.read_csv(cs_file)
    required_cs_cols = ['molecule_chembl_id', 'standard_value', 'canonical_smiles']
    assert all(col in cs_df.columns for col in required_cs_cols), \
        f"Missing columns in source output. Has: {list(cs_df.columns)[:10]}"
    assert len(cs_df) > 0

    # ChEMBL curator checkpoint
    chembl_file = output_dir / f"{run_id}_chembl.csv"
    assert chembl_file.exists()
    chembl_df = pd.read_csv(chembl_file)
    assert 'activity_id' in chembl_df.columns
    assert 'pIC50' in chembl_df.columns  # Output label column

    # Check other checkpoint files exist with basic content validation
    checkpoint_files = {
        'curate': ['curated_smiles', 'curation_success'],
        'tokens': ['tokens', 'seqlen'],
        'confs': ['conformer_success', 'n_conformers'],
    }

    for step, required_cols in checkpoint_files.items():
        file_path = output_dir / f"{run_id}_{step}.csv"
        assert file_path.exists(), f"Checkpoint file for {step} is missing"

        df = pd.read_csv(file_path)
        assert all(col in df.columns for col in required_cols), \
            f"Missing required columns in {step} checkpoint"
        assert len(df) > 0, f"Checkpoint file {step} is empty"

    # Verify checkpoint files were all created (source, chembl, curate, tokens, confs)
    checkpoint_count = len(list(output_dir.glob(f"{run_id}_*.csv")))
    assert checkpoint_count == 5, f"Expected 5 checkpoint files, found {checkpoint_count}"


@pytest.mark.database
def test_pipeline_output_validation(real_pipeline_run):
    """Pipeline output contains expected data from every step."""
    result = real_pipeline_run.result

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0

    # Check columns from each step
    assert 'standard_value' in result.columns   # From ChEMBL source
    assert 'activity_id' in result.columns       # From ChEMBL curator
    assert 'curated_smiles' in result.columns    # From CurateMol
    assert 'tokens' in result.columns            # From TokenizeData
    assert 'conformer_success' in result.columns  # From GenerateConfs
    assert 'n_conformers' in result.columns       # From GenerateConfs

    # Check data validity
    assert not result['curated_smiles'].isnull().any()
    assert not result['tokens'].isnull().any()
    assert result['conformer_success'].sum() > 0
    assert (result['n_conformers'] > 0).sum() > 0


# --------------------------------------------------------------------------- #
# Fast, DB-free tests: context/actor-communication mechanics + error handling
# --------------------------------------------------------------------------- #

@pytest.fixture
def test_output_dir(tmp_path):
    """Temporary directory for test outputs."""
    output_dir = tmp_path / "pipeline_test"
    output_dir.mkdir()
    return str(output_dir)


def _synthetic_df():
    """Tiny hardcoded molecule DataFrame to drive the pipeline without ChEMBL."""
    return pd.DataFrame({
        'canonical_smiles': ['CCO', 'c1ccccc1', 'CC(C)C', 'CCN', 'CCCCO'],
        'molecule_chembl_id': ['CHEMBL1', 'CHEMBL2', 'CHEMBL3', 'CHEMBL4', 'CHEMBL5'],
    })


def _df_pipeline_params(output_root):
    """ForgeParams for a fast DataFrame-fed pipeline (curate -> tokens)."""
    return ForgeParams(
        verbose=False,
        override_actor_params=True,
        return_none_on_fail=False,
        write_checkpoints=False,
        output_root=output_root,
        steps=['curate', 'tokens'],
        curate_params=_curate_params(),
        tokens_params=_tokens_params(),
    )


def test_dataframe_pipeline_context(test_output_dir):
    """Context mechanics work with a synthetic DataFrame input (no ChEMBL DB)."""
    pipeline = ConstructPipe(_df_pipeline_params(test_output_dir))
    result = pipeline(_synthetic_df(), input_id='synthetic')

    context = pipeline.context
    assert isinstance(context, PipelineContext)
    assert context.get('input_id') == 'synthetic'
    assert context.get('input_type') == 'DataFrame'

    # Actors registered and reachable via context
    assert context.get_actor(Steps.CURATE) is not None
    assert context.get_actor(Steps.TOKENS) is not None

    # Per-actor ActorOutput success tracking
    cm_result = context.get_actor_result(Steps.CURATE)
    assert isinstance(cm_result, ActorOutput)
    assert cm_result.success is True

    td_result = context.get_actor_result(Steps.TOKENS)
    assert isinstance(td_result, ActorOutput)
    assert td_result.success is True
    assert 'vocab_size' in td_result.metadata

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0
    assert 'curated_smiles' in result.columns
    assert 'tokens' in result.columns


def test_pipeline_state_consistency(test_output_dir):
    """Each call gets a fresh context reflecting that call's input (no ChEMBL DB)."""
    pipeline = ConstructPipe(_df_pipeline_params(test_output_dir))

    pipeline(_synthetic_df(), input_id='batch_a')
    context1 = pipeline.context

    pipeline(_synthetic_df(), input_id='batch_b')
    context2 = pipeline.context

    # Verify we got a fresh context per run
    assert context1 is not context2
    assert context1.get('input_id') == 'batch_a'
    assert context2.get('input_id') == 'batch_b'


def test_pipeline_error_handling(test_output_dir):
    """Params validation raises FileNotFoundError for a missing DB with no download."""
    with pytest.raises(FileNotFoundError) as exc_info:
        ForgeParams(
            steps=['source', 'chembl', 'curate', 'tokens'],
            override_actor_params=True,
            output_root=test_output_dir,
            return_none_on_fail=False,
            source_params=ChEMBLSourceParams(
                backend='sql',
                auto_download=False,                        # no download
                db_path="/nonexistent/path/chembl.db",      # invalid path
                version=36,
            ),
            chembl_params=ChEMBLCuratorParams(),
            curate_params=_curate_params(),
            tokens_params=_tokens_params(),
        )

    error_msg = str(exc_info.value)
    assert "ChEMBL database not found at" in error_msg
    assert "/nonexistent/path/chembl.db" in error_msg
    assert "auto_download=True" in error_msg


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
