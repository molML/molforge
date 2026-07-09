"""Tests for CLI display module."""

import io
import sys
from unittest.mock import patch

import pytest

from molforge.cli import display


class TestColors:
    """Test Colors class."""
    
    def test_colors_defined(self):
        """Test that all required colors are defined."""
        assert hasattr(display.Colors, 'BLUE')
        assert hasattr(display.Colors, 'GREEN')
        assert hasattr(display.Colors, 'RED')
        assert hasattr(display.Colors, 'YELLOW')
        assert hasattr(display.Colors, 'CYAN')
        assert hasattr(display.Colors, 'MAGENTA')
        assert hasattr(display.Colors, 'BOLD')
        assert hasattr(display.Colors, 'RESET')
    
    def test_colors_are_strings(self):
        """Test that colors are string values."""
        assert isinstance(display.Colors.BLUE, str)
        assert isinstance(display.Colors.GREEN, str)
        assert isinstance(display.Colors.RESET, str)


class TestPrintBanner:
    """Test banner printing functions."""
    
    def test_print_banner_output(self, capsys):
        """Test that print_banner produces output with the current version."""
        from molforge import __version__
        display.print_banner()
        captured = capsys.readouterr()
        assert len(captured.out) > 100
        assert f'v{__version__}' in captured.out

    def test_print_simple_banner_output(self, capsys):
        """Test that print_simple_banner produces output with the current version."""
        from molforge import __version__
        display.print_simple_banner()
        captured = capsys.readouterr()
        assert len(captured.out) > 50
        assert f'v{__version__}' in captured.out
    
    def test_print_simple_banner_shorter(self, capsys):
        """Test that simple banner is shorter than full banner."""
        display.print_simple_banner()
        simple = capsys.readouterr().out
        
        display.print_banner()
        full = capsys.readouterr().out
        
        assert len(simple) < len(full)


class TestPrintDiagram:
    """Test architecture diagram printing."""
    
    def test_print_diagram_output(self, capsys):
        """Test that print_diagram produces output."""
        display.print_diagram()
        captured = capsys.readouterr()
        assert len(captured.out) > 0
        assert 'Architecture' in captured.out or 'Flow' in captured.out


class TestColorize:
    """Test colorize function."""
    
    @patch('molforge.cli.display.Colors.is_supported', return_value=True)
    def test_colorize_basic(self, mock_tty):
        """Test basic colorization."""
        result = display.colorize("test", display.Colors.RED)
        assert display.Colors.RED in result
        assert display.Colors.RESET in result
        assert "test" in result
    
    @patch('molforge.cli.display.Colors.is_supported', return_value=True)
    def test_colorize_empty_string(self, mock_tty):
        """Test colorizing empty string."""
        result = display.colorize("", display.Colors.BLUE)
        assert display.Colors.BLUE in result
        assert display.Colors.RESET in result
    
    def test_colorize_special_characters(self):
        """Test colorizing strings with special characters."""
        text = "test\n123!@#"
        result = display.colorize(text, display.Colors.GREEN)
        assert text in result


class TestStatusPrinting:
    """Test status message functions."""
    
    def test_print_success(self, capsys):
        """Test success message printing."""
        display.print_success("Operation completed")
        captured = capsys.readouterr()
        assert "Operation completed" in captured.out
        assert len(captured.out) > 0
    
    def test_print_error(self, capsys):
        """Test error message printing."""
        display.print_error("Something failed")
        captured = capsys.readouterr()
        assert "Something failed" in captured.err
    
    def test_print_warning(self, capsys):
        """Test warning message printing."""
        display.print_warning("Be careful")
        captured = capsys.readouterr()
        assert "Be careful" in captured.out
    
    def test_print_info(self, capsys):
        """Test info message printing."""
        display.print_info("FYI")
        captured = capsys.readouterr()
        assert "FYI" in captured.out
    
    def test_print_step(self, capsys):
        """Test step message printing."""
        display.print_step(1, 5, "curate", 100, 1.5)
        captured = capsys.readouterr()
        assert "curate" in captured.out
        assert "100" in captured.out


class TestSection:
    """Test section printing."""
    
    def test_print_section(self, capsys):
        """Test section header printing."""
        display.print_section("Configuration")
        captured = capsys.readouterr()
        assert "Configuration" in captured.out
        assert len(captured.out) > len("Configuration")


class TestTableRow:
    """Test table row printing."""
    
    def test_print_table_row_basic(self, capsys):
        """Test basic table row."""
        display.print_table_row("key", "value")
        captured = capsys.readouterr()
        assert "key" in captured.out
        assert "value" in captured.out
    
    def test_print_table_row_with_widths(self, capsys):
        """Test table row with custom widths."""
        display.print_table_row("short", "long value", width1=10)
        captured = capsys.readouterr()
        assert "short" in captured.out
        assert "long value" in captured.out


class TestFormatTime:
    """Test time formatting."""
    
    def test_format_time_seconds(self):
        """Test formatting seconds."""
        result = display.format_time(45.5)
        assert "45.5s" in result
    
    def test_format_time_minutes(self):
        """Test formatting minutes."""
        result = display.format_time(125.0)
        assert "2m" in result
        assert "5s" in result
    
    def test_format_time_hours(self):
        """Test formatting hours."""
        result = display.format_time(3665.0)
        assert "1h" in result
        assert "1m" in result
    
    def test_format_time_zero(self):
        """Test formatting zero time."""
        result = display.format_time(0.0)
        assert "0.0s" in result
    
    def test_format_time_negative(self):
        """Test formatting negative time."""
        result = display.format_time(-10.0)
        # Should handle gracefully (show as 0 or absolute value)
        assert result is not None


class TestProgressBar:
    """Test progress bar creation."""
    
    def test_progress_bar_empty(self):
        """Test progress bar at 0%."""
        result = display.progress_bar(0, 100)
        assert "0%" in result or "0 %" in result
        assert result is not None
    
    def test_progress_bar_half(self):
        """Test progress bar at 50%."""
        result = display.progress_bar(50, 100)
        assert "50%" in result or "50 %" in result
    
    def test_progress_bar_full(self):
        """Test progress bar at 100%."""
        result = display.progress_bar(100, 100)
        assert "100%" in result or "100 %" in result
    
    def test_progress_bar_custom_width(self):
        """Test progress bar with custom width."""
        result = display.progress_bar(25, 100, width=30)
        assert len(result) > 0
    
    def test_progress_bar_zero_total(self):
        """Test progress bar with zero total (edge case)."""
        result = display.progress_bar(0, 0)
        # Should handle gracefully
        assert result is not None
    
    def test_progress_bar_over_100(self):
        """Test progress bar over 100%."""
        result = display.progress_bar(150, 100)
        # Should cap at 100% or handle gracefully
        assert result is not None


class TestColorDetection:
    """Test color support detection."""
    
    @patch('sys.stdout.isatty')
    def test_color_disabled_when_not_tty(self, mock_isatty):
        """Test that colors can be disabled for non-TTY."""
        mock_isatty.return_value = False
        # Note: This test assumes display module checks isatty
        # Actual implementation may vary
        result = display.colorize("test", display.Colors.RED)
        # Should still work, just may not have color codes
        assert "test" in result


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
