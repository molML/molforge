"""
Pipeline context for actor communication and shared state.

Uses step names as the single identifier throughout.
"""

from typing import Any, Optional, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from ..actors.base import BaseActor
    from ..actors.protocol import ActorOutput


class PipelineContext:
    """
    Shared context for actor communication and pipeline state.

    Provides three communication channels:
    1. Actor Registry - Access to upstream actor instances for method calls
    2. Actor Results - Access to upstream outputs (metadata, endpoints)
    3. Shared State - General key-value storage for cross-cutting concerns

    Communication Patterns:
        - Use get_actor() for accessing upstream actor methods (e.g., extract_molecules())
        - Use get_actor_result() for accessing metadata and endpoints
        - Use get()/set() sparingly for shared configuration

    All methods use step names (e.g., 'confs', 'source') as identifiers.
    """

    def __init__(self,
                 run_id: str,
                 output_dir: str,
                 input_id: Optional[str] = None,
                 input_type: Optional[str] = None,
                 input_source: Optional[str] = None):
        """
        Initialize pipeline context.

        Args:
            run_id: Unique identifier for this pipeline run
            output_dir: Directory for pipeline outputs
            input_id: Input identifier
            input_type: Type of input ('ChEMBL ID', 'CSV File', 'DataFrame')
            input_source: Description of input source
        """
        self.run_id = run_id
        self.output_dir = output_dir
        self.input_id = input_id
        self.input_type = input_type
        self.input_source = input_source

        # Actor registry - stores actor instances for downstream method access
        self._actors: Dict[str, 'BaseActor'] = {}

        # Actor results - stores outputs for metadata/endpoint access
        self._actor_results: Dict[str, 'ActorOutput'] = {}

        # Shared state - general key-value storage (include input metadata)
        self._shared_state: Dict[str, Any] = {}

        # Populate shared state with input metadata (for actor access via get_context)
        if input_id is not None:
            self._shared_state['input_id'] = input_id
        if input_type is not None:
            self._shared_state['input_type'] = input_type
        if input_source is not None:
            self._shared_state['input_source'] = input_source

    # ==================== Actor Registry ====================

    def register_actor(self, step_name: str, actor: 'BaseActor') -> None:
        """
        Register actor instance for downstream access.

        Args:
            step_name: Step identifier (e.g., 'confs', 'source')
            actor: Actor instance to register
        """
        self._actors[step_name] = actor

    def get_actor(self, step_name: str) -> Optional['BaseActor']:
        """
        Get upstream actor instance by step name.

        Enables downstream actors to call methods on upstream actors.

        Args:
            step_name: Step identifier (e.g., 'confs', 'source')

        Returns:
            Actor instance if found, None otherwise

        Example:
            >>> from molforge.configuration.steps import Steps
            >>> gc_actor = context.get_actor(Steps.CONFS)
            >>> molecules = gc_actor.extract_molecules()
        """
        return self._actors.get(step_name)

    # ==================== Actor Results ====================

    def set_actor_result(self, step_name: str, result: 'ActorOutput') -> None:
        """
        Store actor output for downstream access.

        Args:
            step_name: Step identifier (e.g., 'confs', 'source')
            result: Actor output containing data, metadata, endpoint
        """
        self._actor_results[step_name] = result

    def get_actor_result(self, step_name: str) -> Optional['ActorOutput']:
        """
        Get upstream actor output by step name.

        Enables access to metadata and endpoints from upstream actors.

        Args:
            step_name: Step identifier (e.g., 'confs', 'source')

        Returns:
            ActorOutput if found, None otherwise

        Example:
            >>> from molforge.configuration.steps import Steps
            >>> gc_result = context.get_actor_result(Steps.CONFS)
            >>> backend = gc_result.metadata['backend']
            >>> endpoint = gc_result.endpoint
        """
        return self._actor_results.get(step_name)

    # ==================== Shared State ====================

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get value from shared state.

        Use sparingly for cross-cutting concerns.

        Args:
            key: State key
            default: Default value if not found

        Returns:
            Value if found, default otherwise
        """
        return self._shared_state.get(key, default)

    def set(self, key: str, value: Any) -> None:
        """
        Store value in shared state.

        Use sparingly for cross-cutting concerns.

        Args:
            key: State key
            value: Value to store
        """
        self._shared_state[key] = value
