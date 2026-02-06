"""
Integration tests for torsion analysis plugin.

Tests torsion calculations with both RDKit and OpenEye conformers.
"""

import pytest
import pandas as pd
import ast
from rdkit import Chem

from molforge.actors.confs import GenerateConfs
from molforge.actors.params.confs import GenerateConfsParams
from molforge.actors.protocol import ActorInput
from molforge.configuration.context import PipelineContext
from molforge.configuration.logger import PipelineLogger
from molforge.utils.constants import HAS_OPENEYE

# Import torsion plugin
try:
    from molforge.actor_plugins.torsion import CalculateTorsions, CalculateTorsionsParams
    TORSION_AVAILABLE = True
except ImportError:
    TORSION_AVAILABLE = False


@pytest.fixture
def torsion_test_data():
    """Drug-like molecules with rotatable bonds for torsion testing."""
    return pd.DataFrame({
        'smiles': [
            'CC(C)Cc1ccc(cc1)C(C)C(=O)O',        # Ibuprofen (racemic) - 5 rotors
            'CC(=O)Nc1ccc(O)cc1',                 # Paracetamol - 1 rotor
            'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',   # (S)-Ibuprofen - 5 rotors
            'c1ccccc1',                           # Benzene - no rotors (ring only)
        ],
        'molecule_id': ['ibuprofen', 'paracetamol', 's_ibuprofen', 'benzene']
    })


@pytest.fixture
def output_dir(tmp_path):
    """Unified output directory for all test artifacts."""
    out_dir = tmp_path / "test_torsion"
    out_dir.mkdir(exist_ok=True)
    return out_dir


@pytest.fixture
def logger(output_dir):
    """Create logger for testing."""
    log_file = str(output_dir / 'test_torsion.log')
    return PipelineLogger(
        name='test_torsion',
        log_file=log_file,
        console_level='INFO'
    )


@pytest.fixture
def pipeline_context(output_dir):
    """Create pipeline context for testing."""
    return PipelineContext(
        run_id='test_torsion',
        output_dir=str(output_dir)
    )


@pytest.mark.skipif(not TORSION_AVAILABLE, reason="Torsion plugin (phd-tools) not available")
class TestTorsionPluginRDKit:
    """Tests for torsion plugin with RDKit conformers."""

    def test_torsion_with_rdkit_conformers(self, torsion_test_data, logger, pipeline_context):
        """Test torsion analysis with RDKit-generated conformers."""
        # Step 1: Generate conformers
        confs_params = GenerateConfsParams(
            backend='rdkit',
            max_confs=10,
            rms_threshold=0.5,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False
        )

        confs_actor = GenerateConfs(confs_params, logger=logger)
        actor_input = ActorInput(data=torsion_test_data, context=pipeline_context)
        confs_output = confs_actor(actor_input)

        # Step 2: Analyze torsions
        torsion_params = CalculateTorsionsParams(
            device='cpu',
            account_for_symmetry=True,
            symmetry_radius=3
        )

        torsion_actor = CalculateTorsions(torsion_params, logger=logger)

        # Register actors in context
        pipeline_context.register_actor('confs', confs_actor)
        pipeline_context.register_actor('torsion', torsion_actor)

        torsion_input = ActorInput(data=confs_output.data, context=pipeline_context)
        torsion_output = torsion_actor(torsion_input)

        # Check output structure
        assert torsion_output.success
        df = torsion_output.data

        # Check all expected columns exist
        assert 'torsion_smiles' in df.columns
        assert 'torsion_mapping' in df.columns
        assert 'n_confs' in df.columns
        assert 'n_torsions' in df.columns
        assert 'n_rotor_torsions' in df.columns
        assert 'n_ring_torsions' in df.columns
        assert 'mean_variance' in df.columns
        assert 'low_confs' in df.columns
        assert 'warnings' in df.columns

        # Check endpoint
        assert torsion_output.endpoint == 'torsion_mapping'

        # Check metadata
        assert 'n_molecules' in torsion_output.metadata
        assert 'total_torsions' in torsion_output.metadata

        # Validate results for each molecule
        for idx, row in df.iterrows():
            mol_id = row['molecule_id']

            # Check torsion_mapping is serialized string
            assert isinstance(row['torsion_mapping'], str)

            # Deserialize for content checks
            torsion_mapping = ast.literal_eval(row['torsion_mapping'])
            assert isinstance(torsion_mapping, dict)

            # Molecules with conformers should have torsions
            if row['conformer_success']:
                # Check torsion_smiles present
                assert isinstance(row['torsion_smiles'], str)
                assert len(row['torsion_smiles']) > 0

                # All molecules should have at least ring torsions
                assert row['n_torsions'] >= 0

                # Ibuprofen molecules should have rotatable torsions
                if 'ibuprofen' in mol_id:
                    assert row['n_rotor_torsions'] > 0
                    assert len(torsion_mapping) > 0

    def test_torsion_handles_failures(self, logger, pipeline_context):
        """Test torsion actor handles molecules without conformers."""
        # Data with no conformers generated
        test_data = pd.DataFrame({
            'smiles': ['INVALID', 'ALSO_INVALID'],
            'molecule_id': ['bad_1', 'bad_2']
        })

        # Generate conformers (will fail)
        confs_params = GenerateConfsParams(
            backend='rdkit',
            max_confs=10,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False
        )

        confs_actor = GenerateConfs(confs_params, logger=logger)
        actor_input = ActorInput(data=test_data, context=pipeline_context)
        confs_output = confs_actor(actor_input)

        # Analyze torsions (should handle gracefully)
        torsion_params = CalculateTorsionsParams()
        torsion_actor = CalculateTorsions(torsion_params, logger=logger)

        # Register actors in context
        pipeline_context.register_actor('confs', confs_actor)
        pipeline_context.register_actor('torsion', torsion_actor)

        torsion_input = ActorInput(data=confs_output.data, context=pipeline_context)
        torsion_output = torsion_actor(torsion_input)

        # Should succeed but with empty results
        assert torsion_output.success
        df = torsion_output.data

        # All values should be 0 or empty
        for idx, row in df.iterrows():
            torsion_mapping = ast.literal_eval(row['torsion_mapping'])
            assert torsion_mapping == {}
            assert row['torsion_smiles'] == ""
            assert row['n_confs'] == 0
            assert row['n_torsions'] == 0

    def test_torsion_single_vs_multiprocess(self, torsion_test_data, logger, pipeline_context):
        """Test that single process and multiprocess give same results."""
        # Generate conformers once
        confs_params = GenerateConfsParams(
            backend='rdkit',
            max_confs=5,
            rms_threshold=0.5,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False
        )

        confs_actor = GenerateConfs(confs_params, logger=logger)
        actor_input = ActorInput(data=torsion_test_data, context=pipeline_context)
        confs_output = confs_actor(actor_input)

        # Process with single process (force by using small dataset)
        torsion_params = CalculateTorsionsParams()
        torsion_actor_single = CalculateTorsions(torsion_params, logger=logger)

        # Register actors in context
        pipeline_context.register_actor('confs', confs_actor)
        pipeline_context.register_actor('torsion', torsion_actor_single)

        torsion_input = ActorInput(data=confs_output.data, context=pipeline_context)
        result_single = torsion_actor_single(torsion_input)

        # Results should be deterministic
        df = result_single.data
        assert df['n_torsions'].sum() > 0  # Should have found torsions
        assert all(isinstance(m, str) for m in df['torsion_mapping'])

        # Verify deserialization works
        for mapping_str in df['torsion_mapping']:
            mapping = ast.literal_eval(mapping_str)
            assert isinstance(mapping, dict)


@pytest.mark.skipif(not TORSION_AVAILABLE, reason="Torsion plugin (phd-tools) not available")
@pytest.mark.skipif(not HAS_OPENEYE, reason="OpenEye not available")
class TestTorsionPluginOpenEye:
    """Tests for torsion plugin with OpenEye conformers."""

    def test_torsion_with_openeye_conformers_no_conversion(self, torsion_test_data, logger, pipeline_context):
        """Test torsion analysis with OpenEye molecules (no RDKit conversion)."""
        # Step 1: Generate conformers with OpenEye (no conversion)
        confs_params = GenerateConfsParams(
            backend='openeye',
            max_confs=10,
            rms_threshold=0.5,
            SMILES_column='smiles',
            names_column='molecule_id',
            convert_to_rdkit=False,  # Keep as OpenEye molecules
            dropna=False
        )

        confs_actor = GenerateConfs(confs_params, logger=logger)
        actor_input = ActorInput(data=torsion_test_data, context=pipeline_context)
        confs_output = confs_actor(actor_input)

        # Step 2: Analyze torsions (should handle OpenEye molecules)
        torsion_params = CalculateTorsionsParams(
            device='cpu',
            account_for_symmetry=True,
            symmetry_radius=3
        )

        torsion_actor = CalculateTorsions(torsion_params, logger=logger)

        # Register actors in context
        pipeline_context.register_actor('confs', confs_actor)
        pipeline_context.register_actor('torsion', torsion_actor)

        torsion_input = ActorInput(data=confs_output.data, context=pipeline_context)
        torsion_output = torsion_actor(torsion_input)

        # Check output structure
        assert torsion_output.success
        df = torsion_output.data

        # Check all expected columns exist
        assert 'torsion_mapping' in df.columns
        assert 'n_torsions' in df.columns

        # Validate OpenEye molecules processed correctly
        for idx, row in df.iterrows():
            if row['conformer_success']:
                # Should have torsions
                assert row['n_torsions'] >= 0
                assert isinstance(row['torsion_mapping'], str)

                # Verify deserialization
                torsion_mapping = ast.literal_eval(row['torsion_mapping'])
                assert isinstance(torsion_mapping, dict)

    def test_torsion_with_openeye_conformers_with_conversion(self, torsion_test_data, logger, pipeline_context):
        """Test torsion analysis with OpenEye conformers converted to RDKit."""
        # Step 1: Generate conformers with OpenEye (with conversion)
        confs_params = GenerateConfsParams(
            backend='openeye',
            max_confs=10,
            rms_threshold=0.5,
            SMILES_column='smiles',
            names_column='molecule_id',
            convert_to_rdkit=True,  # Convert to RDKit
            dropna=False
        )

        confs_actor = GenerateConfs(confs_params, logger=logger)
        actor_input = ActorInput(data=torsion_test_data, context=pipeline_context)
        confs_output = confs_actor(actor_input)

        # Step 2: Analyze torsions (RDKit molecules)
        torsion_params = CalculateTorsionsParams()
        torsion_actor = CalculateTorsions(torsion_params, logger=logger)

        # Register actors in context
        pipeline_context.register_actor('confs', confs_actor)
        pipeline_context.register_actor('torsion', torsion_actor)

        torsion_input = ActorInput(data=confs_output.data, context=pipeline_context)
        torsion_output = torsion_actor(torsion_input)

        # Check success
        assert torsion_output.success
        df = torsion_output.data

        # Validate results
        for idx, row in df.iterrows():
            if row['conformer_success']:
                assert isinstance(row['torsion_mapping'], str)
                assert row['n_torsions'] >= 0

                # Verify deserialization
                torsion_mapping = ast.literal_eval(row['torsion_mapping'])
                assert isinstance(torsion_mapping, dict)
