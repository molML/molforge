"""Standardized abstract base class for all MolForge actors."""

from typing import Optional, List, Any, Union
from abc import ABC, abstractmethod
from pathlib import Path
import pandas as pd

from .protocol import ActorInput, ActorOutput
from .params.base import BaseParams

class BaseActor(ABC):
    """
    Standardized base class for all actors.

    All actors should inherit from this class and implement:
    - process() method: Core processing logic
    - __step_name__ (class attr): Actor identifier for logging
    - required_columns (optional): Input column requirements
    - output_columns (optional): Output column specification
    - __initial__ (optional): Set True to declare the step initial (first)
    - __terminal__ (optional): Set True to declare the step terminal (last)
    """

    __initial__: bool = False
    """Set True in a subclass to declare the step initial: it is expected to be
    the first step in the pipeline. The pipeline configuration warns when an
    initial step is not placed first."""

    __terminal__: bool = False
    """Set True in a subclass to declare the step terminal: it operates on the
    final dataset and is expected to be the last step in the pipeline. The
    pipeline configuration warns when a terminal step is not placed last."""

    def __init__(self, params: BaseParams, logger: Optional[Any] = None, suppress_init_logs: bool = False):
        """
        Standard initialization for all actors.

        Args:
            params: Actor-specific parameters
            logger: PipelineLogger for logging (None for standalone/print to terminal)
            suppress_init_logs: If True, suppress logging during __post_init__.
                               Useful for multiprocessing workers to reduce noise.
        """
        self._params = params
        self.logger = logger
        self._context = None  # Set by pipeline at execution time
        self._suppress_init_logs = suppress_init_logs
        self._in_init_phase = True  # Track if we're in __post_init__

        self._init_from_params(params)

        # Call post-init hook if subclass defines it
        if hasattr(self, '__post_init__'):
            self.__post_init__()

        self._in_init_phase = False  # Init complete

    def _init_from_params(self, params: 'BaseParams') -> None:
        """Initialize actor attributes from params (standardized pattern)."""
        for key, value in params.to_dict().items():
            setattr(self, key, value)

    @property
    def _nametag(self) -> str:
        """
        Actor identifier for logging with consistent width.

        Auto-derived from class __step_name__ attribute. If actor has a backend,
        automatically appends backend name (e.g., 'CONFS-RDK', 'SOURCE-SQL').

        Returns:
            Formatted logging prefix (e.g., '[CONFS-RDK]', '[  SOURCE  ]')
        """
        from ..configuration.steps import Steps
        width = max(Steps.max_length(), 7)  # min 7 for backend tags
        step_name = getattr(self.__class__, '__step_name__', None)

        # Check if actor has backend instance with backend name
        if hasattr(self, 'backend_instance') and self.backend_instance:
            backend_name = getattr(self.backend_instance.__class__, '__backend_name__', None)
            if backend_name and step_name:
                tag = f"{step_name}-{backend_name}"
                return f"[{tag.upper():^{width}s}]"

        # Standard tag without backend
        if step_name:
            return f"[{step_name.upper():^{width}s}]"
        return f"[{self.__class__.__name__[:width].upper():^{width}s}]"

    @property
    def required_columns(self) -> List[str]:
        """
        Columns required in input DataFrame.

        Override in subclass to specify requirements.
        """
        return []

    @property
    def output_columns(self) -> List[str]:
        """
        Columns added/modified by this actor.

        Override in subclass to document output.
        """
        return []

    @property
    def required_context_keys(self) -> List[str]:
        """
        Context keys required by this actor.

        Override in subclass to specify requirements.
        Default: empty list (standard DataFrame processing)

        Common values:
            [] - Standard DataFrame processing (default)
            ['input_id'] - Requires target_id for data source actors

        Returns:
            List of context key names
        """
        return []

    def validate_input(self, data: pd.DataFrame) -> bool:
        """
        Validate input DataFrame has required columns.

        Args:
            data: Input DataFrame

        Returns:
            True if valid, False otherwise
        """
        if not self.required_columns:
            return True

        missing = set(self.required_columns) - set(data.columns)
        if missing:
            self.log(f"Missing required columns: {missing}", level='ERROR')
            return False
        return True

    @abstractmethod
    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Core processing logic (to be implemented by subclass).

        This is the main method actors implement.

        Args:
            data: Input DataFrame

        Returns:
            Processed DataFrame
        """
        ...

    def __call__(self, input_data: ActorInput) -> ActorOutput:
        """
        Standard execution wrapper.

        Handles validation, logging, error handling.
        Actors should NOT override this - override process() instead.

        Args:
            input_data: ActorInput object containing data and context

        Returns:
            ActorOutput with processed data, success status, and metadata
        """
        # Extract data and set context
        data = input_data.data
        self._context = input_data.context

        # Update backend context if backend exists (for actors with backends)
        if hasattr(self, 'backend_instance') and self.backend_instance:
            self.backend_instance.context = self._context

        # Validation
        if not self.validate_input(data):
            return ActorOutput(
                data=data,
                success=False,
                metadata={'error': 'Input validation failed'}
            )

        # Process
        try:
            processed_data = self.process(data)
            return self._create_output(processed_data)

        except Exception as e:
            self.log(f"Processing failed: {e}", level='ERROR')
            return ActorOutput(
                data=data,
                success=False,
                metadata={'error': str(e), 'exception_type': type(e).__name__}
            )

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """
        Create standard output.

        Override to add custom metadata or endpoint.

        Args:
            data: Processed DataFrame

        Returns:
            Actor output with metadata
        """
        return ActorOutput(
            data=data,
            success=True,
            metadata={'n_rows': len(data)}
        )

    # ==================== Logging ====================

    def log(self, message: str, level: str = 'INFO') -> None:
        """
        Standard logging interface.

        Args:
            message: Log message
            level: Log level (INFO, WARNING, ERROR, DEBUG)
        """
        # Suppress logs during init phase if flag is set
        if self._in_init_phase and self._suppress_init_logs:
            return

        if self.logger:
            self.logger.log(level.lower(), self._nametag, message)
        elif getattr(self, 'verbose', True):
            print(f"{self._nametag} | {level} | {message}")

    # ==================== Context Access ====================

    def get_actor(self, attr_name: str) -> Optional['BaseActor']:
        """
        Get upstream actor instance from pipeline context.

        Enables calling methods on upstream actors (e.g., extract_molecules()).
        This is the preferred way to access complex artifacts like molecules or models.

        Args:
            attr_name: Actor's attr_name (e.g., 'GC' for GenerateConfs)

        Returns:
            Actor instance if available, None otherwise

        Example:
            >>> # Access GenerateConfs actor to extract molecules
            >>> gc_actor = self.get_actor('GC')
            >>> if gc_actor:
            ...     molecules = gc_actor.extract_molecules()
        """
        if self._context:
            return self._context.get_actor(attr_name)
        return None

    def get_actor_result(self, attr_name: str) -> Optional[ActorOutput]:
        """
        Get result (metadata, endpoint) from upstream actor.

        Use this to access metadata and endpoints from upstream actors.
        For method calls on actors, use get_actor() instead.

        Args:
            attr_name: Actor's attr_name (e.g., 'GC' for GenerateConfs)

        Returns:
            Actor output if available, None otherwise

        Example:
            >>> # Access GenerateConfs metadata
            >>> gc_result = self.get_actor_result('GC')
            >>> if gc_result:
            ...     backend = gc_result.metadata['backend']
            ...     success_count = gc_result.metadata['success_count']
        """
        if self._context:
            return self._context.get_actor_result(attr_name)
        return None

    def get_context(self, key: str, default: Any = None) -> Any:
        """
        Get value from shared context state.

        Use sparingly for cross-cutting concerns.
        For actor-specific data, use get_actor() or get_actor_result().

        Args:
            key: Context key
            default: Default value if not found

        Returns:
            Context value or default
        """
        if self._context:
            return self._context.get(key, default)
        return default

    def set_context(self, key: str, value: Any) -> None:
        """
        Store value in shared context state.

        Use sparingly for cross-cutting concerns.

        Args:
            key: Context key
            value: Value to store
        """
        if self._context:
            self._context.set(key, value)

    # ==================== File Utilities ====================

    def _get_run_path(self, filename: str) -> str:
        """
        Get full path for a file in current run directory.

        Uses context.output_dir as the run directory.

        Args:
            filename: File name

        Returns:
            Full path to file
        """
        if not self._context or not hasattr(self._context, 'output_dir'):
            raise RuntimeError(
                "Cannot get run path: context not set or output_dir not available. "
                "This method should only be called during actor execution."
            )
        return str(Path(self._context.output_dir) / filename)

    # ==================== Utilities ====================

    def _format_duration(self, seconds: float) -> str:
        """Format duration in human-readable format."""
        if seconds < 60:
            return f"{seconds:.2f}s"
        elif seconds < 3600:
            minutes, secs = divmod(seconds, 60)
            return f"{int(minutes)}m {secs:.1f}s"
        else:
            hours, remainder = divmod(seconds, 3600)
            minutes, secs = divmod(remainder, 60)
            return f"{int(hours)}h {int(minutes)}m {secs:.0f}s"

    @property
    def forge_endpoint(self) -> Any:
        """
        Endpoint for MolForge integration.

        Override in subclass to provide meaningful endpoint
        (file path, column name, or in-memory objects).

        Returns:
            Endpoint for downstream processing, or None
        """
        return None
