"""
Automatic backend discovery system for MolForge.

This module provides functionality to automatically discover and register
backends based on naming conventions, eliminating the need for manual
registration calls in __init__.py files.

Naming Convention:
- File: <backend_name>.py (e.g., sql.py, rdkit.py)
- Class: Must end with 'Backend' (e.g., ChEMBLSQLBackend, RDKitBackend)
- One backend class per file

Example:
    >>> from pathlib import Path
    >>> auto_register_backends('confs', Path(__file__).parent / 'confs')
    ✓ Registered: confs/rdkit → RDKitBackend
    ✓ Registered: confs/openeye → OpenEyeBackend
"""

import importlib
import inspect
from pathlib import Path
from typing import Optional

from .registry import BackendRegistry, Backend


def auto_register_backends(step_name: str, package_path: Path, verbose: bool = True):
    """
    Auto-discover and register backends in a package.

    Args:
        step_name: Step name for namespace (e.g., 'source', 'confs')
        package_path: Path to backend package directory
        verbose: Print registration messages

    Returns:
        Number of backends successfully registered

    Example:
        >>> # In backends/confs/__init__.py
        >>> from pathlib import Path
        >>> from ..discovery import auto_register_backends
        >>> auto_register_backends('confs', Path(__file__).parent)
    """
    if not package_path.exists() or not package_path.is_dir():
        if verbose:
            print(f"⚠ Backend directory not found: {package_path}")
        return 0

    package_name = package_path.name
    registered_count = 0

    # Scan for Python files
    for file in sorted(package_path.glob('*.py')):
        # Skip special files
        if file.stem in ['__init__', 'base', 'discovery']:
            continue

        backend_name = file.stem  # e.g., 'sql', 'rdkit'

        try:
            # Import module
            module = importlib.import_module(
                f'molforge.backends.{package_name}.{backend_name}'
            )

            # Find backend class (must end with 'Backend')
            backend_class = _find_backend_class(module)

            if backend_class:
                # Register with BackendRegistry
                BackendRegistry.register(step_name, backend_name, backend_class)
                registered_count += 1

                if verbose:
                    desc = getattr(backend_class, 'description', '')
                    desc_str = f" - {desc}" if desc else ""
                    print(f"✓ Registered: {step_name}/{backend_name} → {backend_class.__name__}{desc_str}")
            else:
                if verbose:
                    print(f"⚠ No backend class found in {file.name}")

        except Exception as e:
            if verbose:
                print(f"✗ Failed to register {backend_name}: {e}")

    return registered_count


def _find_backend_class(module):
    """
    Find backend class in module.

    Looks for a class that:
    - Ends with 'Backend'
    - Is defined in the module (not imported)
    - Inherits from Backend
    - Is not the Backend base class itself

    Args:
        module: Imported Python module

    Returns:
        Backend class or None if not found
    """
    for name, obj in inspect.getmembers(module, inspect.isclass):
        if (name.endswith('Backend') and
            obj.__module__ == module.__name__ and
            issubclass(obj, Backend) and
            obj is not Backend):
            return obj

    return None


def validate_backend_metadata(backend_class) -> tuple[bool, Optional[str]]:
    """
    Validate that a backend class has proper metadata.

    Checks for:
    - description attribute
    - required_params list
    - optional_params list

    Args:
        backend_class: Backend class to validate

    Returns:
        (is_valid, error_message) tuple
    """
    if not hasattr(backend_class, 'description'):
        return False, f"{backend_class.__name__} missing 'description' attribute"

    if not hasattr(backend_class, 'required_params'):
        return False, f"{backend_class.__name__} missing 'required_params' attribute"

    if not hasattr(backend_class, 'optional_params'):
        return False, f"{backend_class.__name__} missing 'optional_params' attribute"

    # Check types
    if not isinstance(backend_class.description, str):
        return False, f"{backend_class.__name__}.description must be a string"

    if not isinstance(backend_class.required_params, (list, tuple)):
        return False, f"{backend_class.__name__}.required_params must be a list"

    if not isinstance(backend_class.optional_params, (list, tuple)):
        return False, f"{backend_class.__name__}.optional_params must be a list"

    return True, None
