"""
ChEMBL data source backends with automatic registration.

This module provides backends for ChEMBL data retrieval and registers
them with the BackendRegistry via auto-discovery.
"""

from pathlib import Path
from ..discovery import auto_register_backends

# Auto-discover and register all backends in this package
auto_register_backends('source', Path(__file__).parent, verbose=False)

# Import for external use
from .sql import ChEMBLSQLBackend
from .api import ChEMBLAPIBackend

__all__ = ['ChEMBLSQLBackend', 'ChEMBLAPIBackend']
