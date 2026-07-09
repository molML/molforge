from typing import Any, List, Dict, Optional
import importlib

from .registry import ActorRegistry
from .params import PipeParams
from ..configuration.logger import PipelineLogger


class ActorFactory:
    """Factory for creating and configuring pipeline actors."""

    def __init__(self, params: PipeParams, logger: Optional[PipelineLogger] = None):
        self.params = params
        self.logger = logger

    def create_actor(self, step_name: str) -> Optional[Any]:
        """
        Create a single actor for the given step.

        Args:
            step_name: Step identifier (e.g., 'confs', 'source')

        Returns:
            Actor instance or None if creation failed
        """
        # Auto-import backends for this step if they exist
        try:
            importlib.import_module(f'molforge.backends.{step_name}')
        except ModuleNotFoundError:
            pass  # No backends for this step

        actor_info = ActorRegistry.get_actor_info(step_name)
        if not actor_info:
            return None

        actor_class = actor_info['class']
        is_plugin = actor_info.get('is_plugin', False)

        # Get params from appropriate location
        if is_plugin:
            actor_params = self.params.plugin_params.get(step_name)
        else:
            param_attr = ActorRegistry.get_param_attr(step_name)
            actor_params = getattr(self.params, param_attr, None)

        if actor_params is None:
            return None

        # Create actor (step name comes from class attribute)
        actor = actor_class(actor_params, self.logger)

        return actor

    def create_actors(self, step_names: List[str]) -> Dict[str, Any]:
        """
        Create actors for all specified steps.

        Args:
            step_names: List of step identifiers

        Returns:
            Dictionary mapping step_name to actor instance
        """
        actors = {}

        for step_name in step_names:
            actor = self.create_actor(step_name)
            if actor is not None:
                actors[step_name] = actor

        return actors