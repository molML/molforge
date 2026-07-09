"""
Conformer generation backends with automatic registration.

This module provides backends for conformer generation and registers
them with the BackendRegistry via auto-discovery.
"""

from pathlib import Path
from ..discovery import auto_register_backends

# Auto-discover and register all backends in this package
auto_register_backends('confs', Path(__file__).parent, verbose=False)

# Import for external use
from .rdkit import RDKitBackend
from .openeye import OpenEyeBackend

__all__ = ['RDKitBackend', 'OpenEyeBackend']
