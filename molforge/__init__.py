"""MolForge - Molecular processing pipeline."""

from .forge import MolForge
from .configuration.params import ForgeParams
from .configuration.pipe import ConstructPipe

# Re-export common parameter classes for convenience
from .actors.params import (
    ChEMBLSourceParams,
    ChEMBLCuratorParams,
    CurateMolParams,
    TokenizeDataParams,
    CurateDistributionParams,
    GenerateConfsParams,
)
from .actors.params.distributions import PropertyThreshold

__all__ = [
    'MolForge',
    'ForgeParams',
    'ConstructPipe',
    'ChEMBLSourceParams',
    'ChEMBLCuratorParams',
    'CurateMolParams',
    'TokenizeDataParams',
    'CurateDistributionParams',
    'GenerateConfsParams',
    'PropertyThreshold',
]

from ._version import __version__