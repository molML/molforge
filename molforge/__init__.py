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
]

__version__ = '2.0.0'