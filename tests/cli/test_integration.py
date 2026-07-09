"""Integration tests for MolForge CLI."""

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

from molforge.cli.main import main, create_parser


class TestCLIMain:
    """Test main CLI entry point."""
    
    def test_main_no_args(self, capsys):
        """Test CLI with no arguments shows help."""
        result = main([])
        captured = capsys.readouterr()
        assert result == 0
        assert 'molforge' in captured.out.lower()
    
    def test_main_help(self, capsys):
        """Test --help flag."""
        with pytest.raises(SystemExit) as exc_info:
            main(['--help'])
        assert exc_info.value.code == 0
    
    def test_main_version(self, capsys):
        """Test --version flag."""
        result = main(['--version'])
        captured = capsys.readouterr()
        assert result == 0
        assert 'Version' in captured.out or 'version' in captured.out.lower()
    
    def test_main_diagram(self, capsys):
        """Test --diagram flag."""
        result = main(['--diagram'])
        captured = capsys.readouterr()
        assert result == 0
        assert len(captured.out) > 0


class TestCLIRunCommand:
    """Test run command."""
    
    @patch('molforge.cli.commands.run_pipeline')
    def test_run_command_basic(self, mock_run):
        """Test basic run command."""
        mock_run.return_value = 0
        result = main(['run', 'CHEMBL234'])
        assert result == 0
        mock_run.assert_called_once()
    
    @patch('molforge.cli.commands.run_pipeline')
    def test_run_command_with_config(self, mock_run):
        """Test run command with config file."""
        mock_run.return_value = 0
        result = main(['run', 'CHEMBL234', '--config', 'test.yaml'])
        assert result == 0
        
        args = mock_run.call_args[0][0]
        assert args.config == 'test.yaml'
    
    @patch('molforge.cli.commands.run_pipeline')
    def test_run_command_batch(self, mock_run):
        """Test run command with batch flag."""
        mock_run.return_value = 0
        result = main(['run', 'CHEMBL234', 'CHEMBL279', '--batch'])
        assert result == 0
        
        args = mock_run.call_args[0][0]
        assert args.batch is True
        assert len(args.inputs) == 2
    
    @patch('molforge.cli.commands.run_pipeline')
    def test_run_command_with_output(self, mock_run):
        """Test run command with output directory."""
        mock_run.return_value = 0
        result = main(['run', 'CHEMBL234', '--output', './results'])
        assert result == 0
        
        args = mock_run.call_args[0][0]
        assert args.output == './results'


class TestCLIInitCommand:
    """Test init command."""
    
    @patch('molforge.cli.commands.init_config')
    def test_init_command_default(self, mock_init):
        """Test init command with defaults.

        The template now defaults to None so that init_config can decide between
        launching the interactive wizard and copying the 'basic' template.
        """
        mock_init.return_value = 0
        result = main(['init'])
        assert result == 0

        args = mock_init.call_args[0][0]
        assert args.template is None
        assert args.interactive is False
    
    @patch('molforge.cli.commands.init_config')
    def test_init_command_template(self, mock_init):
        """Test init command with specific template."""
        mock_init.return_value = 0
        result = main(['init', '--template', 'full'])
        assert result == 0
        
        args = mock_init.call_args[0][0]
        assert args.template == 'full'
    
    @patch('molforge.cli.commands.init_config')
    def test_init_command_output(self, mock_init):
        """Test init command with output path."""
        mock_init.return_value = 0
        result = main(['init', '--output', 'custom.yaml'])
        assert result == 0
        
        args = mock_init.call_args[0][0]
        assert args.output == 'custom.yaml'


class TestCLIInfoCommand:
    """Test info command."""
    
    @patch('molforge.cli.commands.show_info')
    def test_info_command_default(self, mock_info):
        """Test info command with defaults."""
        mock_info.return_value = 0
        result = main(['info'])
        assert result == 0
    
    @patch('molforge.cli.commands.show_info')
    def test_info_command_actors(self, mock_info):
        """Test info command for actors."""
        mock_info.return_value = 0
        result = main(['info', 'actors'])
        assert result == 0
        
        args = mock_info.call_args[0][0]
        assert args.actors is True
    
    @patch('molforge.cli.commands.show_info')
    def test_info_command_backends(self, mock_info):
        """Test info command for backends."""
        mock_info.return_value = 0
        result = main(['info', 'backends'])
        assert result == 0
        
        args = mock_info.call_args[0][0]
        assert args.backends is True
    
    @patch('molforge.cli.commands.show_info')
    def test_info_command_examples(self, mock_info):
        """Test info command for examples."""
        mock_info.return_value = 0
        result = main(['info', 'examples'])
        assert result == 0
        
        args = mock_info.call_args[0][0]
        assert args.examples is True


class TestCLIValidateCommand:
    """Test validate command."""
    
    @patch('molforge.cli.commands.validate_config')
    def test_validate_command_basic(self, mock_validate):
        """Test basic validate command."""
        mock_validate.return_value = 0
        result = main(['validate', 'config.yaml'])
        assert result == 0
        
        args = mock_validate.call_args[0][0]
        assert args.config == 'config.yaml'
    
    @patch('molforge.cli.commands.validate_config')
    def test_validate_command_quiet(self, mock_validate):
        """Test validate command with quiet flag."""
        mock_validate.return_value = 0
        result = main(['validate', 'config.yaml', '--quiet'])
        assert result == 0
        
        args = mock_validate.call_args[0][0]
        assert args.quiet is True


class TestCLIErrorHandling:
    """Test error handling."""
    
    @patch('molforge.cli.commands.run_pipeline')
    def test_keyboard_interrupt(self, mock_run, capsys):
        """Test handling of KeyboardInterrupt."""
        mock_run.side_effect = KeyboardInterrupt()
        result = main(['run', 'CHEMBL234'])
        assert result == 1
        captured = capsys.readouterr()
        assert 'cancelled' in captured.out.lower()
    
    @patch('molforge.cli.commands.run_pipeline')
    def test_unexpected_exception(self, mock_run, capsys):
        """Test handling of unexpected exceptions."""
        mock_run.side_effect = RuntimeError("Test error")
        result = main(['run', 'CHEMBL234'])
        assert result == 1
        captured = capsys.readouterr()
        assert 'error' in captured.err.lower()


class TestArgumentParser:
    """Test argument parser."""
    
    def test_parser_creation(self):
        """Test that parser can be created."""
        parser = create_parser()
        assert parser is not None
    
    def test_parser_run_command(self):
        """Test parsing run command."""
        parser = create_parser()
        args = parser.parse_args(['run', 'CHEMBL234'])
        assert args.command == 'run'
        assert args.inputs == ['CHEMBL234']
    
    def test_parser_init_command(self):
        """Test parsing init command."""
        parser = create_parser()
        args = parser.parse_args(['init', '--template', 'full'])
        assert args.command == 'init'
        assert args.template == 'full'
    
    def test_parser_info_command(self):
        """Test parsing info command."""
        parser = create_parser()
        args = parser.parse_args(['info', 'actors'])
        assert args.command == 'info'
        assert args.target == 'actors'
    
    def test_parser_validate_command(self):
        """Test parsing validate command."""
        parser = create_parser()
        args = parser.parse_args(['validate', 'config.yaml'])
        assert args.command == 'validate'
        assert args.config == 'config.yaml'
    
    def test_parser_global_flags(self):
        """Test parsing global flags."""
        parser = create_parser()
        
        args = parser.parse_args(['--version'])
        assert args.version is True
        
        args = parser.parse_args(['--diagram'])
        assert args.diagram is True


class TestCLIEndToEnd:
    """End-to-end CLI tests (if possible without actual molforge execution)."""
    
    @patch('molforge.cli.commands.MolForge')
    def test_full_workflow_init_and_run(self, mock_forge_class, tmp_path):
        """Test full workflow: init config then run."""
        # Initialize config
        config_path = tmp_path / "test.yaml"
        result = main(['init', '--output', str(config_path)])
        assert result == 0
        assert config_path.exists()
        
        # Mock forge execution
        mock_forge = mock_forge_class.return_value
        mock_forge.forge.return_value = (None, True, None)
        mock_forge._pipe.output_dir = './output'
        
        # Run pipeline (would fail without mocking)
        # This is just to test argument flow
        with patch('molforge.cli.commands.load_config'):
            result = main(['run', 'test_input', '--config', str(config_path)])
            # Result depends on mock setup


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
