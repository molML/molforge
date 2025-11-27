"""
MolForge constants: Toolkit support detection and package-wide settings.

Actor-specific parameters should be defined in their respective param classes,
not here.
"""

from multiprocessing import cpu_count

# =============================================================================
# Toolkit Support Detection
# =============================================================================

# OpenEye support (optional dependency)
try:
    from openeye import oechem
    HAS_OPENEYE = True
except ImportError:
    HAS_OPENEYE = False
    oechem = None

# Export for centralized import management
__all__ = ['HAS_OPENEYE', 'oechem', 'DEFAULT_N_JOBS', 'DEFAULT_MP_THRESHOLD', 'MAX_CHUNK_SIZE', 'DEFAULT_LOG_LEVEL']

# =============================================================================
# Multiprocessing Configuration
# =============================================================================

# Default number of parallel jobs
# Set to all available CPUs minus 1 to leave headroom for system
DEFAULT_N_JOBS = max(1, cpu_count() - 1)

# Threshold for triggering multiprocessing in actors
# If dataset size >= this value, actors will use multiprocessing
# Used by: CurateMol, CurateDistribution
DEFAULT_MP_THRESHOLD = 1_000

# Maximum items per chunk to avoid IPC overhead
# Large chunks can cause excessive pickling/serialization overhead
# Used by: CurateMol, CurateDistribution
MAX_CHUNK_SIZE = 50_000

# =============================================================================
# Logging Configuration
# =============================================================================

# Default logging level for pipeline operations
DEFAULT_LOG_LEVEL = 'INFO'
