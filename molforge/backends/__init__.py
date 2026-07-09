"""
Backend system for pluggable actor implementations.

This package provides a registry system for managing different backend
implementations of actor functionality (e.g., different data sources,
conformer generators, etc.).
"""

from .registry import BackendRegistry, Backend

__all__ = ['BackendRegistry', 'Backend']
