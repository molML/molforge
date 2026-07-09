from typing import List, Dict, Any
from dataclasses import dataclass
from abc import ABC, abstractmethod


@dataclass
class BaseParams(ABC):
    """Base class for all parameter dataclasses."""
    verbose: bool = True
    """Emit informational (verbose) logging for this component."""

    def __post_init__(self):
        """Validate parameters after initialization."""
        self._validate_params()
        self._post_init_hook()
    
    @abstractmethod
    def _validate_params(self) -> None:
        pass
    
    def _post_init_hook(self) -> None:
        pass
    
    def to_dict(self) -> Dict[str, Any]:
        return {k: v for k, v in self.__dict__.items() if not k.startswith('_')}
    
    def _validate_policy(self, param_name: str, value: str, allowed_values: List[str|None]) -> None:
        if value not in allowed_values:
            raise ValueError(f"Invalid {param_name}='{value}'. Allowed values: {allowed_values}")

