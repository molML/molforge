from typing import Any, List, Type, Union, Optional
import pandas as pd

import json
import hashlib
import traceback
import copy
import time

from datetime import datetime
from pathlib import Path

from .factory import *
from ..actors.protocol import ActorInput, ActorOutput
from .context import PipelineContext
from .logger import PipelineLogger
from .input import InputParser


class ConstructPipe:
    """Main pipeline class using the actor registry and factory for configration of steps."""

    def __init__(self, params: PipeParams, **kwds) -> None:
        self._params = params

        # also set all params as attributes
        for key, value in params.to_dict().items():
            setattr(self, key, value)

        self.logger = None
        self._ensure_dir(self.output_root)

        # initialize logger for construction phase (console only)
        # overwrite existing logs.
        self.log_path = Path(self.output_root) / f"INIT_ID{params._get_config_hash()}.log"
        if self.log_path.exists():
            self.log_path.unlink()
    
        self.logger = PipelineLogger(
            name="pipeline_construction",
            log_file=self.log_path,  # No file logging during construction
            console_level=self.console_log_level,
            file_level=self.file_log_level,
            use_colors=True,
        )
        
        self._initialize_actors()
        self.print(f"Pipe order configured as: {' > '.join(self.steps)}")


    @property
    def method(self) -> str:
        """
        pipeline identifier for logging with consistent width.
        """
        from ..configuration.steps import Steps        
        width = max(Steps.max_length(), 7)  # min 7 for backend tags
        step_name = self._params._get_config_hash() if hasattr(self, 'run_id') else 'PIPELINE'
        
        return f"[{step_name.upper():^{width}s}]"
    
    
    def _initialize_actors(self) -> None:
        """Initialize all actors using the actor factory."""
        self.print('Configuring Actors.')

        factory = ActorFactory(self._params, self.logger)
        actors = factory.create_actors(self.steps)

        for step_name, actor in actors.items():
            setattr(self, step_name, actor)
            self.print(f"  Initialized {step_name}: {type(actor).__name__}")
    
    def __call__(self,
                 input_data: Union[str, pd.DataFrame] = None,
                 input_id: Optional[str] = None,
                 return_info: bool = False,
                 **kwds) -> pd.DataFrame:
        """
        Execute the pipeline with flexible input options.

        Args:
            input_data: Either a ChEMBL ID string, Path to a CSV file, or a pandas DataFrame to use as input
            input_id: Required when input_data is a DataFrame, optional name/identifier
            return_info: If True, return (output, success, step) tuple

            Ensure that SMILES_column is correctly configured for pipe-steps if deviating from default configuration.

        Returns:
            Processed DataFrame or None if failed (or tuple with step info if return_info=True)

        Examples:
            # Option 1: ChEMBL ID - simple string evaluation
            pipeline("CHEMBL123")
            # Output: <output_root>/CHEMBL123_IDabc/...

            # Option 2: CSV file - simple string evaluation
            pipeline("./data/compounds.csv", input_id=optional)
            # Output: <output_root>/compounds_IDabc/...

            # Option 3: DataFrame with name - input type evaluation
            pipeline(df, input_id="my_compounds")
            # Output: <output_root>/my_compounds_IDabc/...
        """
        # Get initialized actors
        actor_list = [getattr(self, step) for step in self.steps]

        # Parse and validate input
        initial_data, parsed_input_id, input_type, input_source = InputParser.parse(
            actor_list, input_data, input_id
        )

        # Store parsed values
        self.input_id = parsed_input_id
        self.input_type = input_type
        self.input_source = input_source
        self.run_id = f"{self.input_id}_ID{self._params._get_config_hash()}"

        # Setup output directory (always under output_root)
        self.output_dir = Path(self.output_root) / self.run_id

        self._ensure_dir(self.output_dir)

        # Create pipeline context (logger is separate)
        self.context = PipelineContext(
            run_id=self.run_id,
            output_dir=str(self.output_dir),
            input_id=self.input_id,
            input_type=self.input_type,
            input_source=self.input_source
        )

        # Save config
        self._save_config()

        # Setup logging
        self.log_path = Path(self.output_dir) / f"{self.run_id}.log"
        if self.log_path.exists():
            self.log_path.unlink()
        
        if self.logger.update_logger(log_file=str(self.log_path)):
            self.print(f"Logging to {self.log_path}")
        else:
            self.print("Log path was updated to None. Logging to console only.", level='WARNING')

        # Execute pipeline
        self.print('Starting pipe.')
        self.print(f"Input ({self.input_id}): {self.input_type} from {self.input_source} ({len(initial_data) if isinstance(initial_data, pd.DataFrame) else 0} rows)")

        output, success, step = self._execute_pipeline(initial_data, self.context)

        # handle final output.
        if not success and self.return_none_on_fail:
            self.print(f'Pipeline failed at step: {step}. Returning None.')
            output = None

        if output is not None and self.write_output:
            self._save_output(output, f"{self.run_id}.csv")

        return (output, success, step) if return_info else output
    
    def _execute_pipeline(self, initial_data: Any, context: PipelineContext) -> tuple[Any, bool, str]:
        """Execute all pipeline steps with timing and context."""
        output, success, last_step, n_steps = initial_data, True, 'start', len(self.steps)
        pipeline_start_time = time.time()
        
        for idx, step in enumerate(self.steps):
            prog = idx * 2
            if not success:
                break

            # Start timing
            self.print(f"{'█'*prog + '░'*(n_steps*2 - prog)} [{idx+1}/{n_steps}] Starting {step}.")
            step_start_time = time.time()

            # Execute step with context
            (output, success), last_step = self._execute_step(output, step, context), step

            # Record information for logging
            step_duration = self._format_duration(time.time() - step_start_time)
            status = 'success' if success else 'failed'
            length = len(output) if hasattr(output, '__len__') else 0

            # Optionally write step checkpoint
            if success and self.write_checkpoints:
                self._save_output(output, f"{self.run_id}_{step}.csv")

            # Log step info
            self.print(f"{'█'*(prog+1) + '░'*(n_steps*2 - (prog+1))} [{idx+1}/{n_steps}] Finished {step} | {status} | {length} rows | {step_duration}")

        # Total timing
        self.print(f"{'█'*(n_steps*2)} [{int(prog/2) + 1}/{n_steps}] Pipeline completed | Total time: {self._format_duration(time.time() - pipeline_start_time)}")
        return output, success, last_step

    def _execute_step(self, input_data: Any, step: str, context: PipelineContext) -> tuple[Any, bool]:
        """Execute a single pipeline step with context."""
        if not hasattr(self, step):
            self.print(f"ERROR: Unknown step '{step}'. Check registry.")
            return None, False

        actor = getattr(self, step)

        try:
            # Execute with ActorInput
            actor_input = ActorInput(data=input_data, context=context)
            actor_output = actor(actor_input)

            # Register actor for downstream access
            context.register_actor(step, actor)

            # Store result in context
            context.set_actor_result(step, actor_output)

            # Return data and success status
            success = actor_output.success and isinstance(actor_output.data, pd.DataFrame) and len(actor_output.data) > 0
            return actor_output.data, success

        except Exception as e:
            self._log_error(actor, e)
            return None, False
        
    def _log_error(self, actor: Any, error: Exception) -> None:
        """Log error with traceback."""
        error_msg = f"Error in {actor.method}: {type(error).__name__}: {error}"
        self.print(error_msg, level="ERROR")
    
        if self.verbose:
            tb_str = traceback.format_exc()
            self.print(f"Full traceback:\n{tb_str}", level="DEBUG")
            
    def _ensure_dir(self, dir_path: Union[str, Path]) -> None:
        """Ensure directory exists."""
        path = Path(dir_path)
        
        if not path.is_dir():
            try:
                path.mkdir(parents=True, exist_ok=True)
                self.print(f'Created output dir: {path}')
            except Exception as e:
                self.print(f'Failed to create dir {path}: {e}')
                raise
        else:
            self.print(f'Using existing dir: {path}')
    
    def _save_output(self, df: pd.DataFrame, filename: str) -> None:
        """Save DataFrame to output directory."""
        output_path = Path(self.output_dir) / filename
        df.to_csv(output_path, index=False)
    
    def _save_config(self) -> None:
        """Save configuration to JSON file."""
        config = self._get_params_dict()
        config.update({
            # TODO: add coarse timing (month/year) to config hash for chembl api versioning.
            'log_time': datetime.now().strftime("%Y-%m-%d %H:%M:%S"), 
            'input_type': self.input_type,
            'input_id': self.input_id,
            'input_source': self.input_source,
            'run_id': self.run_id
        })
        
        config_path = Path(self.output_dir) / f"{self.run_id}_config.json"
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=4, default=str)

        return str(config_path)

    def _get_params_dict(self) -> dict:
        """Convert params to nested dictionary."""
        def to_dict(obj):
            if hasattr(obj, 'to_dict'):
                result = obj.to_dict()
            elif isinstance(obj, dict):
                result = copy.deepcopy(obj)
            else:
                return obj
            
            # rcursively convert nested objects.
            return {k: to_dict(v) for k, v in result.items()}
        
        return to_dict(self._params)
    
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

    def print(self, *values, level: str = None, **kwargs):
        """
        Print/log method for pipeline messages.

        Args:
            *values: Values to print/log
            level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
            **kwargs: Additional arguments
        """
        message = ' '.join(str(v) for v in values)

        # Auto-detect level if not specified
        if level is None:
            level = 'INFO'

        if self.logger and self.verbose:
            self.logger.log(level.lower(), self.method, message)
        elif self.verbose:
            print(self.method, message, **kwargs)

    def get_actor(self, identifier: str) -> Optional[Any]:
        """
        Get an initialized actor by step name or attr_name.

        Args:
            identifier: Step name (e.g., 'confs', 'curate') or attr_name (e.g., 'GC', 'CM')

        Returns:
            Actor instance or None if not found

        Note:
            This method now delegates to ActorRegistry.get_actor_instance()
            for centralized actor retrieval logic.
        """
        return ActorRegistry.get_actor_instance(self, identifier)
