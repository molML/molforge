import logging
import sys
from typing import Optional

class PipelineLogger:
    """
    Centralized logger for pipeline actors that supports both file logging and console output.
    Provides colored output for console and timestamped, method-tagged logs.
    """
    
    # ANSI color codes for console output
    COLORS = {
        'DEBUG': '\033[36m',      # Cyan
        'INFO': '\033[32m',       # Green  
        'WARNING': '\033[33m',    # Yellow
        'ERROR': '\033[31m',      # Red
        'CRITICAL': '\033[35m',   # Magenta
        'RESET': '\033[0m'        # Reset
    }
    
    def __init__(self, 
                 name: str = "pipeline",
                 log_file: Optional[str] = None,
                 console_level: str = "INFO",
                 file_level: str = "INFO",
                 use_colors: bool = True):
        """
        Initialize pipeline logger.
        
        Args:
            name: Logger name
            log_file: Path to log file (if None, only console logging)
            console_level: Minimum level for console output
            file_level: Minimum level for file output
            use_colors: Whether to use colored console output
        """
        self.name = name
        self.log_file = log_file
        self.console_level = console_level
        self.file_level = file_level
        self.use_colors = use_colors
        
        self.fmt = '%(asctime)s    | %(method)s | %(levelname)s | %(message)s'
        self.datefmt = '%Y-%m-%d %H:%M:%S'

        self._setup_logger()
    
    def _setup_logger(self):
        """Setup the actual logging infrastructure."""
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.DEBUG)
        self.use_colors_enabled = self.use_colors and sys.stdout.isatty()
        
        # Clear any existing handlers
        self.logger.handlers.clear()
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, self.console_level.upper()))
        
        if self.use_colors_enabled:
            console_formatter = ColoredFormatter(
                fmt=self.fmt,
                datefmt=self.datefmt
            )
        else:
            console_formatter = BaseFormatter(
                fmt=self.fmt,
                datefmt=self.datefmt
            )
        
        console_handler.setFormatter(console_formatter)
        self.logger.addHandler(console_handler)
        
        # File handler (if specified)
        if self.log_file:
            self._add_file_handler()
    
    def _add_file_handler(self):
        """Add file handler with current log_file path."""
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(getattr(logging, self.file_level.upper()))
        file_formatter = BaseFormatter(
            fmt=self.fmt,
            datefmt=self.datefmt
        )
        file_handler.setFormatter(file_formatter)
        self.logger.addHandler(file_handler)
    
    def _remove_file_handlers(self):
        """Remove all existing file handlers."""
        handlers_to_remove = [h for h in self.logger.handlers if isinstance(h, logging.FileHandler)]
        for handler in handlers_to_remove:
            handler.close()  # Important: close the handler to release file resources
            self.logger.removeHandler(handler)
    
    def update_logger(self, log_file: str = None):
        """Update logger to write to a new log file."""
        self._remove_file_handlers()
        
        if log_file:
            self.log_file = log_file
            self._add_file_handler()
            return True
        else:
            self.log_file = None
            return False
                
        
    def log(self, level: str, method: str, message: str, *args, **kwargs):
        """Internal logging method that adds method context."""
        # Create a LogRecord with custom 'method' field
        extra = {'method': method}
        getattr(self.logger, level.lower())(message, *args, extra=extra, **kwargs)
    
class BaseFormatter(logging.Formatter):
    """Base formatter that abbreviates log levels and handles method field."""

    def format(self, record):
        # Abbreviate log level to 4 characters for consistent spacing
        original_levelname = record.levelname
        record.levelname = f"{original_levelname:<8s}"  # Left-aligned, 8 chars

        # Format the message
        formatted = super().format(record)

        # Restore original level name
        record.levelname = original_levelname

        return formatted


class ColoredFormatter(BaseFormatter):
    """Formatter that adds colors to console output."""

    def format(self, record):
        # Get the base formatted message
        formatted = super().format(record)

        # Add color to the entire line
        color = PipelineLogger.COLORS.get(record.levelname, '')
        reset = PipelineLogger.COLORS['RESET']
        
        if color:
            formatted = f"{color}{formatted}{reset}"

        return formatted