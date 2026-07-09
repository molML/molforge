"""
Actor protocol and data structures for standardized pipeline execution.

Defines the contracts that all actors must follow.
"""

from typing import Any, Optional, Dict, TYPE_CHECKING
from dataclasses import dataclass, field
import pandas as pd

if TYPE_CHECKING:
    from ..configuration.context import PipelineContext


@dataclass
class ActorInput:
    """
    Standard input wrapper for actors.

    Context flows with data as a unit. BaseActor extracts and sets context internally.
    """
    data: pd.DataFrame
    context: 'PipelineContext'


@dataclass
class ActorOutput:
    """
    Standard output wrapper for actors.

    Attributes:
        data: Processed DataFrame
        success: Whether processing succeeded
        metadata: Additional information about processing
        endpoint: Optional endpoint for downstream usage (file path, objects, column name)
    """
    data: pd.DataFrame
    success: bool = True
    metadata: Dict[str, Any] = field(default_factory=dict)
    endpoint: Optional[Any] = None
