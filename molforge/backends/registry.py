"""
Generic backend registry system for MolForge actors.

This module provides a unified interface for registering and accessing
different backend implementations across all actor types (data sources,
conformer generators, property calculators, etc.).

All actors that support multiple backends should use the 'backend' parameter
in their params classes and register their backends here.
"""

from typing import Dict, Type, List, Optional, Any
from abc import ABC, abstractmethod


class Backend(ABC):
    """
    Base class for all backend implementations.

    Backends provide alternative implementations for actor functionality
    (e.g., RDKit vs OpenEye for conformers, SQL vs API for data sources).

    Class Attributes (should be defined in subclasses):
        __backend_name__ (str): Short identifier for logging tags (e.g., 'RDK', 'OEO', 'SQL', 'API')
        description (str): Human-readable description of this backend
        required_params (list): List of required parameter names
        optional_params (list): List of optional parameter names
    """

    # Class-level metadata (override in subclasses)
    __backend_name__: str = None
    """Short backend identifier for logging tags (e.g., 'RDK', 'OEO', 'SQL', 'API')"""

    description: str = ""
    """Human-readable description of this backend"""

    required_params: list = []
    """List of required parameter names"""

    optional_params: list = []
    """List of optional parameter names"""

    @abstractmethod
    def __init__(self, params: Any, logger: Any = None, context: Any = None):
        """
        Initialize backend with parameters.

        Args:
            params: Backend-specific parameters
            logger: Optional logger instance
            context: Optional pipeline context
        """
        pass

    @classmethod
    def validate_params(cls, params) -> None:
        """
        Validate params before instantiation.

        Raises:
            ValueError: If required parameters are missing

        Example:
            >>> ChEMBLSQLBackend.validate_params(params)
        """
        for param in cls.required_params:
            if not hasattr(params, param) or getattr(params, param) is None:
                raise ValueError(
                    f"{cls.__name__} requires '{param}' parameter"
                )


class BackendRegistry:
    """
    Central registry for all backend implementations.

    Similar to ActorRegistry but for backend implementations within actors.
    Supports multiple backend types (source, confs, etc.).
    """

    # Registry: backend_type -> {backend_name -> backend_class}
    _backends: Dict[str, Dict[str, Type[Backend]]] = {}

    @classmethod
    def register(cls, backend_type: str, name: str, backend_class: Type[Backend]) -> None:
        """
        Register a backend implementation.

        Args:
            backend_type: Category of backend (e.g., 'source', 'confs')
            name: Backend identifier (e.g., 'sql', 'api', 'rdkit', 'openeye')
            backend_class: Backend class (must inherit from Backend)

        Raises:
            TypeError: If backend_class doesn't inherit from Backend
            ValueError: If backend already registered

        Example:
            >>> BackendRegistry.register('source', 'sql', ChEMBLSQLBackend)
            >>> BackendRegistry.register('confs', 'rdkit', RDKitBackend)
        """
        if not issubclass(backend_class, Backend):
            raise TypeError(
                f"Backend class must inherit from Backend, "
                f"got {backend_class.__name__}"
            )

        backend_type = backend_type.lower()
        name = name.lower()

        if backend_type not in cls._backends:
            cls._backends[backend_type] = {}

        if name in cls._backends[backend_type]:
            raise ValueError(
                f"Backend '{name}' already registered for type '{backend_type}'. "
                f"Use a different name or unregister first."
            )

        cls._backends[backend_type][name] = backend_class

    @classmethod
    def get(cls, backend_type: str, name: str) -> Type[Backend]:
        """
        Get backend class by type and name.

        Args:
            backend_type: Category of backend (e.g., 'source', 'confs')
            name: Backend identifier (e.g., 'sql', 'rdkit')

        Returns:
            Backend class (not instantiated)

        Raises:
            ValueError: If backend type or name not found

        Example:
            >>> backend_cls = BackendRegistry.get('confs', 'rdkit')
            >>> backend = backend_cls(params, context)
        """
        backend_type = backend_type.lower()
        name = name.lower()

        if backend_type not in cls._backends:
            available_types = ', '.join(cls._backends.keys())
            raise ValueError(
                f"Unknown backend type '{backend_type}'. "
                f"Available types: {available_types}"
            )

        if name not in cls._backends[backend_type]:
            available = ', '.join(cls._backends[backend_type].keys())
            raise ValueError(
                f"Unknown backend '{name}' for type '{backend_type}'. "
                f"Available backends: {available}"
            )

        return cls._backends[backend_type][name]

    @classmethod
    def list_backends(cls, backend_type: Optional[str] = None) -> Dict[str, List[str]]:
        """
        List all available backends.

        Args:
            backend_type: If provided, only list backends for this type

        Returns:
            Dict mapping backend types to lists of backend names,
            or dict with single type if backend_type specified

        Example:
            >>> BackendRegistry.list_backends()
            {'source': ['sql', 'api'], 'confs': ['rdkit', 'openeye']}
            >>> BackendRegistry.list_backends('confs')
            {'confs': ['rdkit', 'openeye']}
        """
        if backend_type is not None:
            backend_type = backend_type.lower()
            if backend_type not in cls._backends:
                return {}
            return {backend_type: list(cls._backends[backend_type].keys())}

        return {
            btype: list(backends.keys())
            for btype, backends in cls._backends.items()
        }

    @classmethod
    def unregister(cls, backend_type: str, name: str) -> None:
        """
        Unregister a backend.

        Args:
            backend_type: Category of backend
            name: Backend identifier

        Raises:
            ValueError: If backend not found
        """
        backend_type = backend_type.lower()
        name = name.lower()

        if backend_type not in cls._backends:
            raise ValueError(f"Backend type '{backend_type}' not registered")

        if name not in cls._backends[backend_type]:
            raise ValueError(
                f"Backend '{name}' not registered for type '{backend_type}'"
            )

        del cls._backends[backend_type][name]

    @classmethod
    def has_backend(cls, backend_type: str, name: str) -> bool:
        """
        Check if a backend is registered.

        Args:
            backend_type: Category of backend
            name: Backend identifier

        Returns:
            True if backend is registered, False otherwise
        """
        backend_type = backend_type.lower()
        name = name.lower()

        return (
            backend_type in cls._backends and
            name in cls._backends[backend_type]
        )

    @classmethod
    def get_backend_type_info(cls, backend_type: str) -> Dict[str, Type[Backend]]:
        """
        Get all backends for a specific type.

        Args:
            backend_type: Category of backend

        Returns:
            Dict mapping backend names to backend classes

        Raises:
            ValueError: If backend type not found
        """
        backend_type = backend_type.lower()

        if backend_type not in cls._backends:
            available_types = ', '.join(cls._backends.keys())
            raise ValueError(
                f"Unknown backend type '{backend_type}'. "
                f"Available types: {available_types}"
            )

        return cls._backends[backend_type].copy()
