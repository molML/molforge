"""
Integration tests for unified conformer generation actor.

Tests both RDKit and OpenEye backends with the unified GenerateConfs actor.
"""

import pytest
import pandas as pd
from rdkit import Chem


from molforge.actors.confs import GenerateConfs
from molforge.actors.params.confs import GenerateConfsParams
from molforge.actors.protocol import ActorInput
from molforge.configuration.context import PipelineContext
from molforge.configuration.logger import PipelineLogger
from molforge.utils.constants import HAS_OPENEYE


@pytest.fixture
def test_data():
    """Sample SMILES data for testing."""
    return pd.DataFrame({
        'smiles': [
            'CCO',  # Ethanol - simple
            'c1ccccc1',  # Benzene - aromatic
            'CC(C)C',  # Isobutane - branched
            'INVALID',  # Invalid SMILES
        ],
        'molecule_id': ['mol_1', 'mol_2', 'mol_3', 'mol_4']
    })


@pytest.fixture
def conformer_test_data():
    """Drug-like molecules with rotatable bonds for conformer generation testing."""
    return pd.DataFrame({
        'smiles': [
            'CC(C)Cc1ccc(cc1)C(C)C(=O)O',        # Ibuprofen (racemic)
            'CC(=O)Nc1ccc(O)cc1',                 # Paracetamol/Acetaminophen
            'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',   # (S)-Ibuprofen (with stereochemistry)
            'INVALID_SMILES',                     # Invalid SMILES
        ],
        'molecule_id': ['ibuprofen', 'paracetamol', 's_ibuprofen', 'invalid']
    })


@pytest.fixture
def output_dir(tmp_path):
    """Unified output directory for all test artifacts."""
    out_dir = tmp_path / "test_output"
    out_dir.mkdir(exist_ok=True)
    return out_dir

@pytest.fixture
def logger(output_dir):
    """Create logger for testing with unified output dir."""
    log_file = str(output_dir / 'test_confs_unified.log')
    return PipelineLogger(
        name='test_confs_unified',
        log_file=log_file,
        console_level='INFO'
    )




@pytest.fixture
def pipeline_context(output_dir):
    """Create pipeline context for testing with unified output dir."""
    return PipelineContext(
        run_id='test_confs_unified',
        output_dir=str(output_dir)
    )


class TestRDKitBackend:
    """Tests for RDKit backend via unified actor."""

    def test_rdkit_conformer_generation(self, test_data, logger, pipeline_context):
        """Test basic RDKit conformer generation."""
        params = GenerateConfsParams(
            backend='rdkit',
            max_confs=50,
            rms_threshold=0.5,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False
        )

        actor = GenerateConfs(params, logger=logger)
        actor_input = ActorInput(data=test_data, context=pipeline_context)

        output = actor(actor_input)

        # Check output structure
        assert output.success
        assert 'n_conformers' in output.data.columns
        assert 'conformer_success' in output.data.columns
        assert 'Status' in output.data.columns

        # Check metadata
        assert output.metadata['backend'] == 'rdkit'
        assert output.metadata['molecules_processed'] == 4
        assert output.metadata['success_count'] == 3  # 3 valid, 1 invalid

        # Check results
        df = output.data
        assert df.loc[df['molecule_id'] == 'mol_1', 'conformer_success'].iloc[0]
        assert df.loc[df['molecule_id'] == 'mol_2', 'conformer_success'].iloc[0]
        assert df.loc[df['molecule_id'] == 'mol_3', 'conformer_success'].iloc[0]
        assert not df.loc[df['molecule_id'] == 'mol_4', 'conformer_success'].iloc[0]

        # Check conformer counts
        assert df.loc[df['molecule_id'] == 'mol_1', 'n_conformers'].iloc[0] > 0
        assert df.loc[df['molecule_id'] == 'mol_4', 'n_conformers'].iloc[0] == 0

    def test_rdkit_extract_molecules(self, test_data, logger, pipeline_context):
        """Test molecule extraction from RDKit backend."""
        params = GenerateConfsParams(
            backend='rdkit',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=True
        )

        actor = GenerateConfs(params, logger=logger)
        actor_input = ActorInput(data=test_data, context=pipeline_context)

        output = actor(actor_input)

        # Extract molecules
        molecules = list(actor.extract_molecules())

        # Should have 3 successful molecules
        assert len(molecules) == 3

        # Check they are RDKit molecules with conformers
        for mol in molecules:
            assert isinstance(mol, Chem.Mol)
            assert mol.GetNumConformers() > 0
            assert mol.HasProp('_Name')

    def test_rdkit_dropna(self, test_data, logger, pipeline_context):
        """Test dropna parameter with RDKit backend."""
        params = GenerateConfsParams(
            backend='rdkit',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=True
        )

        actor = GenerateConfs(params, logger=logger)
        actor_input = ActorInput(data=test_data, context=pipeline_context)

        output = actor(actor_input)

        # Should only have successful molecules
        assert len(output.data) == 3
        assert output.data['conformer_success'].all()

    def test_rdkit_legacy_pattern(self, test_data, logger, pipeline_context):
        """Test that actor only accepts ActorInput (legacy removed)."""
        params = GenerateConfsParams(
            backend='rdkit',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False
        )

        actor = GenerateConfs(params, logger=logger)

        # New pattern: must use ActorInput
        from molforge.actors.protocol import ActorInput, ActorOutput
        actor_input = ActorInput(data=test_data, context=pipeline_context)
        result = actor(actor_input)

        # Should return ActorOutput (not raw DataFrame)
        assert isinstance(result, ActorOutput)
        assert 'conformer_success' in result.data.columns
        assert 'n_conformers' in result.data.columns
        assert 'Status' in result.data.columns

    def test_rdkit_pickle_file_integrity(self, conformer_test_data, logger, pipeline_context):
        """
        Verify RDKit pickle file integrity and consistency with DataFrame results.

        Tests:
        1. Pickle file exists at expected path
        2. File contains (binary, name) tuples in correct format
        3. Pickle molecule names match conformer_success=True rows in DataFrame
        4. Molecules can be reconstructed from pickle with accessible conformers
        5. Invalid SMILES are properly excluded from pickle file
        """
        import pickle
        from pathlib import Path

        params = GenerateConfsParams(
            backend='rdkit',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False
        )

        actor = GenerateConfs(params, logger=logger)
        actor_input = ActorInput(data=conformer_test_data, context=pipeline_context)

        output = actor(actor_input)

        # 1. Verify pickle file exists at expected path
        expected_pickle_path = Path(pipeline_context.output_dir) / "RDKit_conformers.pkl"
        assert expected_pickle_path.exists(), "Pickle file should exist"

        # 2. Read pickle file directly and verify format
        with open(expected_pickle_path, 'rb') as f:
            mol_data = pickle.load(f)

        # Verify it's a list of tuples
        assert isinstance(mol_data, list), "Pickle should contain a list"
        assert all(isinstance(item, tuple) and len(item) == 2 for item in mol_data), \
            "Each item should be a (binary, name) tuple"

        # 3. Extract successful names from DataFrame and compare with pickle
        successful_df_names = output.data[output.data['conformer_success']]['molecule_id'].tolist()
        pickle_names = [name for _, name in mol_data]

        assert set(pickle_names) == set(successful_df_names), \
            "Pickle file names should match conformer_success=True names in DataFrame"

        # 4. Verify invalid SMILES is excluded
        assert 'invalid' not in pickle_names, "Invalid SMILES should not be in pickle file"
        assert 'invalid' in output.data['molecule_id'].tolist(), "Invalid SMILES should be in DataFrame"
        invalid_row = output.data[output.data['molecule_id'] == 'invalid']
        assert not invalid_row['conformer_success'].iloc[0], "Invalid SMILES should have conformer_success=False"

        # 5. Reconstruct molecules from binary and verify conformers
        for mol_binary, name in mol_data:
            mol = Chem.Mol(mol_binary)
            assert mol is not None, f"Should be able to reconstruct molecule {name}"
            assert mol.GetNumConformers() > 0, f"Molecule {name} should have conformers"

            # Verify name is preserved
            if mol.HasProp('_Name'):
                assert mol.GetProp('_Name') == name


class TestOpenEyeBackend:
    """Tests for OpenEye backend via unified actor."""

    @pytest.mark.skipif(
        not HAS_OPENEYE,
        reason="OpenEye toolkit not available"
    )
    def test_openeye_conformer_generation(self, test_data, logger, pipeline_context):
        """Test basic OpenEye conformer generation."""
        params = GenerateConfsParams(
            backend='openeye',
            max_confs=50,
            rms_threshold=0.5,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=False,
            mode='classic'
        )

        actor = GenerateConfs(params, logger=logger)
        actor_input = ActorInput(data=test_data, context=pipeline_context)

        output = actor(actor_input)

        # Check output structure
        assert output.success
        assert 'n_conformers' in output.data.columns
        assert 'conformer_success' in output.data.columns

        # Check metadata
        assert output.metadata['backend'] == 'openeye'
        assert output.metadata['molecules_processed'] == 4

        # Check endpoint is file path
        assert isinstance(output.endpoint, str)
        assert output.endpoint.endswith('.oeb')

    @pytest.mark.skipif(
        not HAS_OPENEYE,
        reason="OpenEye toolkit not available"
    )
    def test_openeye_extract_molecules(self, test_data, logger, pipeline_context):
        """Test molecule extraction from OpenEye backend."""
        params = GenerateConfsParams(
            backend='openeye',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=True
        )

        actor = GenerateConfs(params, logger=logger)
        actor_input = ActorInput(data=test_data, context=pipeline_context)

        output = actor(actor_input)

        # Extract molecules
        molecules = list(actor.extract_molecules())

        # Should have successful molecules
        assert len(molecules) > 0

        # Check they are RDKit molecules with conformers
        for mol in molecules:
            assert isinstance(mol, Chem.Mol)
            assert mol.GetNumConformers() > 0

    @pytest.mark.skipif(
        not HAS_OPENEYE,
        reason="OpenEye toolkit not available"
    )
    def test_openeye_convert_to_rdkit_parameter(self, conformer_test_data, logger, pipeline_context):
        """
        Verify convert_to_rdkit parameter controls molecule type from extract_molecules().

        Tests:
        1. convert_to_rdkit=False yields OpenEye molecules with conformers
        2. convert_to_rdkit=True yields RDKit molecules with conformers
        3. Invalid SMILES are properly excluded
        """
        from molforge.utils.constants import oechem

        # Test Case 1: convert_to_rdkit=False (yields OpenEye molecules)
        params_oe = GenerateConfsParams(
            backend='openeye',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=True,
            convert_to_rdkit=False
        )

        actor_oe = GenerateConfs(params_oe, logger=logger)
        actor_input_oe = ActorInput(data=conformer_test_data, context=pipeline_context)

        output_oe = actor_oe(actor_input_oe)

        # Extract molecules and verify they are OpenEye molecules
        molecules_oe = list(actor_oe.extract_molecules())
        assert len(molecules_oe) > 0, "Should have successful molecules"

        for mol in molecules_oe:
            # Verify type is OpenEye molecule
            assert isinstance(mol, oechem.OEMol), "Expected OpenEye OEMol instance"
            assert type(mol).__module__ == 'oechem', \
                f"Expected oechem module, got {type(mol).__module__}"
            # Verify conformers are accessible (OpenEye uses GetMaxConfIdx)
            assert mol.GetMaxConfIdx() >= 0, "OpenEye molecule should have conformers"

        # Test Case 2: convert_to_rdkit=True (yields RDKit molecules)
        params_rdkit = GenerateConfsParams(
            backend='openeye',
            max_confs=20,
            SMILES_column='smiles',
            names_column='molecule_id',
            dropna=True,
            convert_to_rdkit=True
        )

        actor_rdkit = GenerateConfs(params_rdkit, logger=logger)
        actor_input_rdkit = ActorInput(data=conformer_test_data, context=pipeline_context)

        output_rdkit = actor_rdkit(actor_input_rdkit)

        # Extract molecules and verify they are RDKit molecules
        molecules_rdkit = list(actor_rdkit.extract_molecules())
        assert len(molecules_rdkit) > 0, "Should have successful molecules"

        for mol in molecules_rdkit:
            # Verify type is RDKit molecule
            assert isinstance(mol, Chem.Mol), "Expected RDKit molecule"
            assert type(mol).__module__ == 'rdkit.Chem.rdchem', \
                f"Expected RDKit molecule, got {type(mol).__module__}"
            # Verify conformers are accessible (RDKit uses GetNumConformers)
            assert mol.GetNumConformers() > 0, "RDKit molecule should have conformers"


class TestBackendRegistry:
    """Tests for backend registry system."""

    def test_list_backends(self):
        """Test listing available backends."""
        from molforge.backends.registry import BackendRegistry

        backends = BackendRegistry.list_backends()
        assert 'confs' in backends.keys()

        backends = backends['confs']
        assert 'rdkit' in backends
        assert 'openeye' in backends

    def test_get_backend(self):
        """Test getting backend by name."""
        from molforge.backends.confs import RDKitBackend, OpenEyeBackend
        from molforge.backends.registry import BackendRegistry
        
        rdkit_cls = BackendRegistry.get('confs', 'rdkit')
        assert rdkit_cls == RDKitBackend

        openeye_cls = BackendRegistry.get('confs', 'openeye')
        assert openeye_cls == OpenEyeBackend

    def test_invalid_backend(self, test_data, logger, pipeline_context):
        """Test invalid backend name raises error."""
        with pytest.raises(ValueError, match="Unknown backend 'invalid_backend' for type 'confs'"):
            params = GenerateConfsParams(
                backend='invalid_backend',
                SMILES_column='smiles',
                names_column='molecule_id'
            )
            actor = GenerateConfs(params, logger=logger)
