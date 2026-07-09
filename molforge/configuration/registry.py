"""
Actor registry with automatic plugin discovery.

Uses step names as the single identifier throughout the system.
"""

from typing import Any, List, Dict, Type, Optional, Tuple
from pathlib import Path
import importlib
import inspect

from ..actors import *
from ..actors.params import *

from ..backends import BackendRegistry
from .steps import Steps


class ActorRegistry:
    """
    Registry for pipeline actors with plugin support.

    Core actors are explicitly registered using Steps constants.
    Plugin actors are auto-discovered from actor_plugins/ directory.

    Step names are the single identifier used throughout (no attr_names).
    Dependencies reference step names that must appear earlier in the pipeline.
    """

    # Core actors (explicit registration using Steps constants)
    _core_actors: Dict[str, Dict[str, Any]] = {
        Steps.SOURCE: {
            'class': ChEMBLSource,
            'param_class': ChEMBLSourceParams,
            'is_plugin': False,
            'dependencies': None,
        },
        Steps.CHEMBL: {
            'class': ChEMBLCurator,
            'param_class': ChEMBLCuratorParams,
            'is_plugin': False,
            'dependencies': [Steps.SOURCE],  # Requires data source
        },
        Steps.CURATE: {
            'class': CurateMol,
            'param_class': CurateMolParams,
            'is_plugin': False,
            'dependencies': None,
        },
        Steps.TOKENS: {
            'class': TokenizeData,
            'param_class': TokenizeDataParams,
            'is_plugin': False,
            'dependencies': None,
        },
        Steps.DISTRIBUTIONS: {
            'class': CurateDistribution,
            'param_class': CurateDistributionParams,
            'is_plugin': False,
            'dependencies': None,
        },
        Steps.CONFS: {
            'class': GenerateConfs,
            'param_class': GenerateConfsParams,
            'is_plugin': False,
            'dependencies': None,
        }
    }

    # Plugin actors (dynamically discovered)
    _plugin_actors: Dict[str, Dict[str, Any]] = {}

    # Track if plugins have been loaded
    _plugins_loaded: bool = False

    @classmethod
    def discover_plugins(cls, plugin_dir: Optional[str] = None) -> None:
        """
        Discover and register plugin actors.

        Args:
            plugin_dir: Path to plugin directory. If None, uses default actor_plugins/
        """
        if plugin_dir is None:
            # Default: look for actor_plugins/ relative to this file
            current_dir = Path(__file__).parent
            plugin_dir = current_dir.parent / 'actor_plugins'
        else:
            plugin_dir = Path(plugin_dir)

        if not plugin_dir.exists():
            return

        # Scan for Python files
        for plugin_file in plugin_dir.glob('*.py'):
            if plugin_file.name.startswith('_'):
                continue

            try:
                cls._load_plugin_module(plugin_file)
            except Exception as e:
                print(f"Warning: Failed to load plugin {plugin_file.name}: {e}")

        cls._plugins_loaded = True

    @classmethod
    def _load_plugin_module(cls, plugin_file: Path) -> None:
        """Load a plugin module and register its actors."""
        # Dynamic import
        module_name = f"actor_plugins.{plugin_file.stem}"
        spec = importlib.util.spec_from_file_location(module_name, plugin_file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Look for Actor classes (inherit from base actor)
        for _, obj in inspect.getmembers(module, inspect.isclass):
            if cls._is_valid_plugin_actor(obj):
                cls._register_plugin_actor(obj)

    @classmethod
    def _is_valid_plugin_actor(cls, actor_class: Type) -> bool:
        """Check if a class is a valid plugin actor."""
        # Must have required attributes (support both __step_name__ and __plugin_name__)
        has_name = hasattr(actor_class, '__step_name__') or hasattr(actor_class, '__plugin_name__')
        return (
            has_name and
            hasattr(actor_class, '__param_class__') and
            hasattr(actor_class, '__call__')
        )

    @classmethod
    def _register_plugin_actor(cls, actor_class: Type) -> None:
        """Register a discovered plugin actor."""
        # Support both __step_name__ (new) and __plugin_name__ (backward compat)
        step_name = getattr(actor_class, '__step_name__', None) or actor_class.__plugin_name__

        cls._plugin_actors[step_name] = {
            'class': actor_class,
            'param_class': actor_class.__param_class__,
            'is_plugin': True,
            'dependencies': getattr(actor_class, '__dependencies__', None),
        }

    @classmethod
    def register_plugin(cls,
                       step_name: str,
                       actor_class: Type,
                       param_class: Type,
                       dependencies: Optional[List[str]] = None) -> None:
        """
        Manually register a plugin actor.

        Use this for runtime registration without file-based plugins.

        Args:
            step_name: Name of the pipeline step
            actor_class: Actor class
            param_class: Parameter class for the actor
            dependencies: List of required step names (not attr_names)
        """
        cls._plugin_actors[step_name] = {
            'class': actor_class,
            'param_class': param_class,
            'is_plugin': True,
            'dependencies': dependencies,
        }

    @classmethod
    def _ensure_plugins_loaded(cls) -> None:
        """Ensure plugins are discovered before accessing registry."""
        if not cls._plugins_loaded:
            cls.discover_plugins()

    @classmethod
    def get_actor_info(cls, step_name: str) -> Optional[Dict[str, Any]]:
        """Get actor configuration for a given step."""
        cls._ensure_plugins_loaded()

        # Check core actors first
        if step_name in cls._core_actors:
            return cls._core_actors[step_name]

        # Then check plugins
        return cls._plugin_actors.get(step_name)

    @classmethod
    def get_all_step_names(cls) -> List[str]:
        """Get all registered step names (core + plugins)."""
        cls._ensure_plugins_loaded()
        return list(cls._core_actors.keys()) + list(cls._plugin_actors.keys())

    @classmethod
    def get_core_steps(cls) -> List[str]:
        """Get core actor step names only."""
        return list(cls._core_actors.keys())

    @classmethod
    def get_plugin_steps(cls) -> List[str]:
        """Get plugin actor step names only."""
        cls._ensure_plugins_loaded()
        return list(cls._plugin_actors.keys())

    @classmethod
    def is_plugin(cls, step_name: str) -> bool:
        """Check if a step is a plugin."""
        info = cls.get_actor_info(step_name)
        return info.get('is_plugin', False) if info else False

    @classmethod
    def get_param_attr(cls, step_name: str) -> str:
        """
        Get the param attribute name for a step.

        Auto-derived from step name: 'confs' -> 'confs_params'
        """
        return f"{step_name}_params"

    @classmethod
    def get_param_class(cls, step_name: str) -> Optional[Type]:
        """Get the param class for a step."""
        info = cls.get_actor_info(step_name)
        return info['param_class'] if info else None

    @classmethod
    def get_dependencies(cls, step_name: str) -> Optional[List[str]]:
        """
        Get the dependencies (step names) for a step.

        Returns list of step names that must appear earlier in pipeline.
        """
        info = cls.get_actor_info(step_name)
        return info.get('dependencies') if info else None

    @classmethod
    def has_backends(cls, step_name: str) -> bool:
        """
        Check if a step has registered backends.

        Args:
            step_name: Step name (e.g., 'confs', 'source')

        Returns:
            True if step has backends registered in BackendRegistry
        """
        backends = BackendRegistry.list_backends(step_name)
        return len(backends.get(step_name, [])) > 0

    @classmethod
    def get_actor_instance(cls, owner: Any, step_name: str) -> Optional[Any]:
        """
        Get actor instance from pipeline by step name.

        Args:
            owner: ConstructPipe or MolForge instance that owns actors
            step_name: Step name (e.g., 'confs', 'chembl')

        Returns:
            Actor instance or None if not found

        Example:
            >>> actor = ActorRegistry.get_actor_instance(pipeline, Steps.CONFS)
        """
        cls._ensure_plugins_loaded()

        # Actors are stored by step name on pipeline
        if hasattr(owner, step_name):
            return getattr(owner, step_name)

        return None

    @classmethod
    def validate_steps(cls, steps: List[str]) -> Tuple[bool, Optional[str]]:
        """
        Validate pipeline steps for existence and dependency satisfaction.

        Args:
            steps: List of step names in the pipeline

        Returns:
            (is_valid, error_message) tuple
        """
        cls._ensure_plugins_loaded()

        # Check all steps exist
        for step in steps:
            if not cls.get_actor_info(step):
                return False, f"Unknown step: '{step}'"

        # Check dependencies (must exist and appear earlier)
        available_steps = set()
        for step in steps:
            dependencies = cls.get_dependencies(step)

            if dependencies:
                for dep_step in dependencies:
                    if dep_step not in available_steps:
                        return False, (
                            f"Dependency error: Step '{step}' requires '{dep_step}' "
                            f"which is not available or appears later in the pipeline."
                        )

            # Add current step to available set
            available_steps.add(step)

        return True, None

    @classmethod
    def get_molecule_provider_step(cls, steps: List[str]) -> Optional[str]:
        """
        Get the step name that provides molecules for boxing.

        Checks for conformer generator (confs), then falls back to curated molecules (curate).

        Args:
            steps: List of step names in the pipeline

        Returns:
            Step name that provides molecules ('confs' or 'curate'), or None
        """
        cls._ensure_plugins_loaded()

        # Priority: conformer generator > curated molecules
        if Steps.CONFS in steps:
            return Steps.CONFS
        elif Steps.CURATE in steps:
            return Steps.CURATE

        return None
