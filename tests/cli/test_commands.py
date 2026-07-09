"""Tests for CLI commands module."""

import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import argparse

import pytest

from molforge.cli import commands


class TestLoadConfig:
    """Test configuration loading."""
    
    def test_load_config_basic(self, tmp_path):
        """Test loading basic configuration."""
        config_path = tmp_path / "test_config.yaml"
        config_path.write_text("""
steps:
  - source
  - curate

source:
  backend: api

output:
  dir: ./output
  checkpoints: false
""")
        
        params = commands.load_config(str(config_path))
        assert params.steps == ['source', 'curate']
        assert params.output_root == './output'
    
    def test_load_config_missing_file(self):
        """Test loading non-existent config file."""
        with pytest.raises(FileNotFoundError):
            commands.load_config("/nonexistent/config.yaml")
    
    def test_load_config_invalid_yaml(self, tmp_path):
        """Test loading invalid YAML."""
        config_path = tmp_path / "invalid.yaml"
        config_path.write_text("invalid: yaml: content:")
        
        with pytest.raises(Exception):
            commands.load_config(str(config_path))


class TestGetTemplate:
    """Test template retrieval."""
    
    def test_get_template_basic(self):
        """Test getting basic template."""
        template = commands.get_template('basic')
        assert 'steps:' in template
        assert 'source' in template
        assert 'chembl' in template
    
    def test_get_template_full(self):
        """Test getting full template."""
        template = commands.get_template('full')
        assert 'steps:' in template
        assert 'distributions' in template
        assert 'tokens' in template
    
    def test_get_template_conformers(self):
        """Test getting conformers template."""
        template = commands.get_template('conformers')
        assert 'confs:' in template
        assert 'backend: rdkit' in template
    
    def test_get_template_invalid(self):
        """Test getting invalid template."""
        with pytest.raises(ValueError) as exc_info:
            commands.get_template('nonexistent')
        assert 'Unknown template' in str(exc_info.value)


class TestInitConfig:
    """Test config initialization command."""
    
    def test_init_config_basic(self, tmp_path):
        """Test initializing basic config."""
        output_path = tmp_path / "config.yaml"
        args = argparse.Namespace(
            template='basic',
            output=str(output_path)
        )
        
        result = commands.init_config(args)
        assert result == 0
        assert output_path.exists()
        
        content = output_path.read_text()
        assert 'steps:' in content
        assert 'source' in content
    
    def test_init_config_creates_directories(self, tmp_path):
        """Test that init_config creates parent directories."""
        output_path = tmp_path / "subdir" / "config.yaml"
        args = argparse.Namespace(
            template='full',
            output=str(output_path)
        )
        
        result = commands.init_config(args)
        assert result == 0
        assert output_path.exists()
        assert output_path.parent.exists()
    
    @patch('molforge.cli.commands.HAS_YAML', False)
    def test_init_config_no_yaml(self, capsys):
        """Test init_config without PyYAML."""
        args = argparse.Namespace(template='basic', output='test.yaml')
        result = commands.init_config(args)
        assert result == 1
        captured = capsys.readouterr()
        assert 'PyYAML' in captured.err


class TestValidateConfig:
    """Test config validation command."""
    
    def test_validate_config_valid(self, tmp_path):
        """Test validating valid config."""
        config_path = tmp_path / "valid.yaml"
        config_path.write_text("""
steps:
  - source
  - curate

source:
  backend: sql

output:
  dir: ./output
""")
        
        args = argparse.Namespace(
            config=str(config_path),
            quiet=True
        )
        
        with patch('molforge.cli.commands.load_config') as mock_load:
            mock_load.return_value = Mock(steps=['source', 'curate'])
            result = commands.validate_config(args)
            assert result == 0
    
    def test_validate_config_missing(self):
        """Test validating missing config."""
        args = argparse.Namespace(
            config='/nonexistent.yaml',
            quiet=True
        )
        
        result = commands.validate_config(args)
        assert result == 1
    
    @patch('molforge.cli.commands.HAS_YAML', False)
    def test_validate_config_no_yaml(self, capsys):
        """Test validate without PyYAML."""
        args = argparse.Namespace(config='test.yaml', quiet=True)
        result = commands.validate_config(args)
        assert result == 1
        captured = capsys.readouterr()
        assert 'PyYAML' in captured.err


class TestRunSingle:
    """Test single pipeline execution."""
    
    @patch('molforge.cli.commands.MolForge')
    def test_run_single_success(self, mock_forge_class):
        """Test successful single run."""
        mock_forge = Mock()
        mock_forge.forge.return_value = (Mock(), True, None)
        mock_forge._pipe.output_dir = '/output'
        
        result = commands.run_single(mock_forge, 'CHEMBL234', quiet=True)
        assert result == 0
        mock_forge.forge.assert_called_once()
    
    @patch('molforge.cli.commands.MolForge')
    def test_run_single_failure(self, mock_forge_class):
        """Test failed single run."""
        mock_forge = Mock()
        mock_forge.forge.return_value = (None, False, 'curate')
        
        result = commands.run_single(mock_forge, 'CHEMBL234', quiet=True)
        assert result == 1
    
    @patch('molforge.cli.commands.MolForge')
    def test_run_single_exception(self, mock_forge_class):
        """Test single run with exception."""
        mock_forge = Mock()
        mock_forge.forge.side_effect = RuntimeError("Test error")
        
        result = commands.run_single(mock_forge, 'CHEMBL234', quiet=True)
        assert result == 1


class TestRunBatch:
    """Test batch pipeline execution."""
    
    @patch('molforge.cli.commands.MolForge')
    def test_run_batch_all_success(self, mock_forge_class):
        """Test batch run with all successes."""
        mock_forge = Mock()
        mock_forge.forge.return_value = (Mock(), True, None)
        
        inputs = ['CHEMBL234', 'CHEMBL279']
        result = commands.run_batch(mock_forge, inputs, quiet=True)
        assert result == 0
        assert mock_forge.forge.call_count == 2
    
    @patch('molforge.cli.commands.MolForge')
    def test_run_batch_partial_failure(self, mock_forge_class):
        """Test batch run with partial failures."""
        mock_forge = Mock()
        mock_forge.forge.side_effect = [
            (Mock(), True, None),
            (None, False, 'curate')
        ]
        
        inputs = ['CHEMBL234', 'CHEMBL279']
        result = commands.run_batch(mock_forge, inputs, quiet=True)
        assert result == 1


class TestShowInfo:
    """Test info display command."""
    
    def test_show_info_version(self, capsys):
        """Test showing version info."""
        args = argparse.Namespace(
            version=True,
            actors=False,
            backends=False,
            examples=False
        )
        
        with patch('molforge.cli.commands.print_version'):
            result = commands.show_info(args)
            assert result == 0
    
    def test_show_info_actors(self, capsys):
        """Test showing actors info."""
        args = argparse.Namespace(
            version=False,
            actors=True,
            backends=False,
            examples=False
        )
        
        with patch('molforge.cli.commands.print_actors_info'):
            result = commands.show_info(args)
            assert result == 0
    
    def test_show_info_backends(self, capsys):
        """Test showing backends info."""
        args = argparse.Namespace(
            version=False,
            actors=False,
            backends=True,
            examples=False,
            actor=None
        )
        
        with patch('molforge.cli.commands.print_backends_info'):
            result = commands.show_info(args)
            assert result == 0


class TestPrintFunctions:
    """Test info printing functions."""
    
    @patch('molforge.cli.commands.ActorRegistry')
    def test_print_actors_info(self, mock_registry, capsys):
        """Test printing actors information."""
        mock_registry.get_all_step_names.return_value = ['source', 'curate']
        mock_registry.get_actor_info.return_value = {
            'class': Mock(__name__='TestActor'),
            'is_plugin': False
        }
        
        commands.print_actors_info()
        captured = capsys.readouterr()
        assert len(captured.out) > 0
    
    @patch('molforge.cli.commands.BackendRegistry')
    def test_print_backends_info(self, mock_registry, capsys):
        """Test printing backends information."""
        mock_registry.list_backends.return_value = {
            'source': ['sql', 'api'],
            'confs': ['rdkit', 'openeye']
        }
        
        commands.print_backends_info()
        captured = capsys.readouterr()
        assert len(captured.out) > 0
    
    def test_print_examples(self, capsys):
        """Test printing usage examples."""
        commands.print_examples()
        captured = capsys.readouterr()
        assert 'molforge' in captured.out
        assert len(captured.out) > 0


class TestRunPipeline:
    """Test main run_pipeline command.

    ``build_params`` is patched in these tests because constructing a real
    ForgeParams with a 'source' step defaults to the SQL backend, whose
    validation performs a network lookup / potential multi-GB ChEMBL DB
    download. Patching it keeps the fast suite free of any network I/O while
    still exercising ``run_pipeline`` control flow and wiring.
    """

    @patch('molforge.cli.commands.MolForge')
    @patch('molforge.cli.commands.build_params')
    def test_run_pipeline_with_config(self, mock_build, mock_forge_class, tmp_path):
        """Test running pipeline with config file routed through build_params."""
        config_path = tmp_path / "config.yaml"
        config_path.write_text("steps: [source]")

        mock_params = Mock()
        mock_params.output_root = './output'
        mock_build.return_value = mock_params

        mock_forge = Mock()
        mock_forge.forge.return_value = (Mock(), True, None)
        mock_forge._pipe.output_dir = './output'
        mock_forge_class.return_value = mock_forge

        args = argparse.Namespace(
            config=str(config_path),
            inputs=['CHEMBL234'],
            output=None,
            steps=None,
            set=None,
            batch=False,
            quiet=True
        )

        result = commands.run_pipeline(args)
        assert result == 0
        # config, steps and overrides are all funneled through build_params.
        _, kwargs = mock_build.call_args
        assert kwargs['config_path'] == str(config_path)

    @patch('molforge.cli.commands.MolForge')
    @patch('molforge.cli.commands.build_params')
    def test_run_pipeline_no_config(self, mock_build, mock_forge_class):
        """Test running pipeline without config file (steps + set overrides)."""
        mock_params = Mock()
        mock_params.output_root = './output'
        mock_build.return_value = mock_params

        mock_forge = Mock()
        mock_forge.forge.return_value = (Mock(), True, None)
        mock_forge._pipe.output_dir = './output'
        mock_forge_class.return_value = mock_forge

        args = argparse.Namespace(
            config=None,
            inputs=['CHEMBL234'],
            output=None,
            steps='source,curate',
            set=['curate.dropna=false'],
            batch=False,
            quiet=True
        )

        result = commands.run_pipeline(args)
        assert result == 0
        _, kwargs = mock_build.call_args
        assert kwargs['steps'] == 'source,curate'
        assert kwargs['overrides'] == ['curate.dropna=false']


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
