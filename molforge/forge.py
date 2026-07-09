"""
MolForge - Clean API for the molecular processing pipeline.

Forges raw molecular data through multi-step processing (ChEMBL retrieval,
curation, tokenization, conformer generation) into refined molecular datasets.
"""

from typing import Union, Optional, Tuple, Dict, Any
from pathlib import Path
import pandas as pd

from .configuration.pipe import ConstructPipe
from .configuration.factory import PipeParams
from .configuration.registry import ActorRegistry


class MolForge:
    """
    Clean API for the molecular processing pipeline.

    MolForge processes raw molecular data through configurable steps:
    - ChEMBL data retrieval
    - Data curation and filtering
    - Molecule curation and standardization
    - Tokenization
    - Property computation and distribution filtering
    - Conformer generation

    Results are returned as a pandas DataFrame and can be written to CSV.
    """

    def __init__(self,
                 params: Union[PipeParams, str, Dict] = None,
                 config_path: str = None,
                 output_root: Optional[str] = None,
                 **config_overrides):
        """
        Initialize MolForge pipeline.

        Args:
            params: PipeParams object, or dict of parameters
            config_path: Path to configuration file (YAML/JSON)
            output_root: Root output directory for INIT log and run subdirectories
            **config_overrides: Override specific config parameters

        Examples:
            >>> # From params object
            >>> forge = MolForge(params)

            >>> # From config file
            >>> forge = MolForge(config_path="pipeline_config.yaml")

            >>> # With custom root output directory
            >>> forge = MolForge(params, output_root="projects/kinases")

            >>> # With overrides
            >>> forge = MolForge(params, verbose=True)
        """
        # Handle different initialization methods
        if params is None and config_path is None:
            raise ValueError("Must provide either params or config_path")

        if config_path is not None:
            params = self._load_config(config_path)

        if isinstance(params, dict):
            params = PipeParams(**params)

        # Apply output_root override if provided
        if output_root is not None:
            config_overrides['output_root'] = output_root

        # Apply overrides
        if config_overrides:
            params_dict = params.to_dict()
            params_dict.update(config_overrides)
            params = PipeParams(**params_dict)

        # Initialize internal pipeline
        self._pipe = ConstructPipe(params)
        self._params = params

    def forge(self,
              input_data: Union[str, pd.DataFrame],
              input_id: Optional[str] = None,
              return_info: bool = False) -> Union[pd.DataFrame, Tuple[pd.DataFrame, bool, str]]:
        """
        Forge molecules through the pipeline.

        Args:
            input_data: ChEMBL target ID (str), CSV file path (str), or DataFrame
            input_id: Optional identifier (auto-detected for ChEMBL/CSV)
            return_info: If True, return (dataframe, success, failed_step) tuple

        Returns:
            Processed DataFrame, or tuple if return_info=True

        Examples:
            >>> # From ChEMBL (output to configured output_root)
            >>> df = forge.forge("CHEMBL123")
            >>> # Output: <output_root>/CHEMBL123_IDabc/...

            >>> # From DataFrame
            >>> df = forge.forge(input_df, input_id="my_compounds")
            >>> # Output: <output_root>/my_compounds_IDabc/...

            >>> # With step info
            >>> df, success, step = forge.forge("CHEMBL123", return_info=True)
        """
        return self._pipe(input_data, input_id, return_info=return_info)

    def forge_to_csv(self,
                     input_data: Union[str, pd.DataFrame],
                     output_path: str,
                     input_id: Optional[str] = None) -> Tuple[pd.DataFrame, bool]:
        """
        Forge molecules and save to CSV.

        Args:
            input_data: ChEMBL target ID, CSV path, or DataFrame
            output_path: Path for output CSV file
            input_id: Optional identifier

        Returns:
            Tuple of (dataframe, success)

        Note:
            Pipeline intermediate files are written to output_root (set at MolForge initialization).
            Only the final CSV is written to output_path.

        Example:
            >>> df, success = forge.forge_to_csv("CHEMBL123", "output.csv")
        """
        df, success, _ = self._pipe(input_data, input_id, return_info=True)

        if success and df is not None:
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_path, index=False)

        return df, success

    def batch_forge(self,
                    input_list: list,
                    output_dir: str) -> Dict[str, Tuple[bool, str]]:
        """
        Forge multiple inputs in batch, writing one CSV per input.

        Args:
            input_list: List of inputs (ChEMBL IDs, CSV paths, or DataFrames)
            output_dir: Directory for output CSV files

        Returns:
            Dict mapping input identifiers to (success, output_path) tuples

        Example:
            >>> results = forge.batch_forge(
            ...     ["CHEMBL123", "CHEMBL456"],
            ...     "output_dir"
            ... )
        """
        results = {}
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        for input_data in input_list:
            # Generate identifier
            if isinstance(input_data, str):
                input_id = Path(input_data).stem if input_data.endswith('.csv') else input_data
            else:
                input_id = f"batch_{len(results)}"

            try:
                out_file = output_path / f"{input_id}.csv"
                df, success = self.forge_to_csv(input_data, str(out_file), input_id)
                results[input_id] = (success, str(out_file))

            except Exception as e:
                self._pipe.print(f"Failed to process {input_id}: {e}", level="ERROR")
                results[input_id] = (False, str(e))

        return results

    def _load_config(self, config_path: str) -> PipeParams:
        """Load configuration from a JSON or YAML file via the shared codec."""
        from .configuration.serialization import ConfigCodec

        if not Path(config_path).exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        return ConfigCodec.from_file(config_path, PipeParams)

    @property
    def config(self) -> PipeParams:
        """Access current configuration."""
        return self._params

    @property
    def steps(self) -> list:
        """Get configured pipeline steps."""
        return self._pipe.steps

    def get_actor(self, identifier: str) -> Any:
        """
        Retrieve an initialized actor from the pipe.

        Args:
            identifier: Step name (e.g., 'confs', 'curate')

        Returns:
            Actor instance or None if not found

        Note:
            Delegates to ActorRegistry.get_actor_instance().
        """
        return ActorRegistry.get_actor_instance(self._pipe, identifier)

    def __repr__(self) -> str:
        steps_str = ' → '.join(self.steps)
        return f"MolForge({steps_str})"
