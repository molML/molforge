"""
Unified ChEMBL data source parameters.

Supports multiple backends (sql, api) with shared and backend-specific parameters.
"""

from typing import Literal, Optional
from dataclasses import dataclass, field

import os
from ...utils.actortools.chembl_downloader import get_latest_chembl_version, ensure_chembl_db
from .base import BaseParams


@dataclass
class ChEMBLSourceParams(BaseParams):
    """
    Unified ChEMBL data source parameters.

    Provides a single parameter class that supports multiple backends,
    with both shared parameters (work across all backends) and
    backend-specific parameters (only used by relevant backend).
    """

    # ==================== Backend Selection ====================

    backend: Literal['sql', 'api'] = 'sql'
    """ChEMBL data source backend: 'sql' or 'api'"""

    # ==================== Shared Parameters ====================
    # These parameters work across all backends

    n: int = 1000
    """Maximum number of activity entries to fetch, applied only when search_all is False (ignored when search_all=True; in the API backend it also bounds the per-request page size)."""

    search_all: bool = True
    """Whether to fetch all entries instead of limiting to n"""

    # ==================== SQL-Specific Parameters ====================
    # These are only used when backend='sql'

    db_path: Optional[str] = None
    """Path to the local ChEMBL SQLite database (SQL backend)"""

    version: str | int = 'latest'
    """ChEMBL version to use (SQL backend)"""

    auto_download: bool = True
    """If True, downloads the database if not found (SQL backend)"""

    download_dir: Optional[str] = field(default="./data/chembl")
    """Directory to store the downloaded database (SQL backend)"""

    # ==================== API-Specific Parameters ====================
    # These are only used when backend='api'
    # (Currently API backend uses shared parameters only)

    def _validate_params(self) -> None:
        """Validate parameter values based on selected backend."""
        # Shared validation
        if self.n <= 0:
            raise ValueError("n must be positive")

        # SQL-specific validation
        if self.backend == 'sql':
            self._validate_sql_params()

    def _validate_sql_params(self) -> None:
        """
        Resolve and validate the SQL database, handling auto-download.

        The user-facing ``db_path`` stays as provided. The path used at run time
        is resolved into the private ``_resolved_db_path`` (machine-specific, so
        it stays out of the serialized config). ``version`` is pinned to the
        resolved ChEMBL version when a path is derived automatically, or to a
        relative reference to the file when a path is supplied manually.
        """
        # Create target directory if missing
        os.makedirs(self.download_dir, exist_ok=True)

        if self.db_path is None:
            # Auto: resolve the version and derive the standard path.
            self.version = get_latest_chembl_version() if self.version == 'latest' else self.version
            self._resolved_db_path = os.path.join(self.download_dir, f"{self.version}.db")
        else:
            # Manual: pin the user's file; record a relative reference as version.
            self._resolved_db_path = self.db_path
            self.version = os.path.relpath(self.db_path)

        db_path = self._resolved_db_path

        # Check if database exists and is valid
        if os.path.exists(db_path):
            if not self._check_database_integrity(db_path):
                print(f"Database file corrupted: {db_path}")
                if self.auto_download:
                    print("Removing corrupted file and re-downloading.")
                    os.remove(db_path)
                else:
                    raise RuntimeError(
                        f"ChEMBL database is corrupted at {db_path}. "
                        "Enable `auto_download=True` to re-download automatically."
                    )

        # Download if missing or was corrupted
        if not os.path.exists(db_path):
            if not self.auto_download:
                raise FileNotFoundError(
                    f"ChEMBL database not found at {db_path}. "
                    "Enable `auto_download=True` to download automatically."
                )
            print(f"ChEMBL SQL Database not found. Auto-downloading to {db_path} ...")
            ensure_chembl_db(db_path, version=self.version)

            if not self._check_database_integrity(db_path):
                raise RuntimeError(
                    f"Downloaded database failed integrity check: {db_path}"
                )
    
    def _check_database_integrity(self, db_path: str) -> bool:
        """
        Check if a SQLite database file is valid and not corrupted.
        
        Uses a tiered approach:
        1. Quick header check (fast, catches obvious corruption)
        2. SQLite connection test (catches structural issues)
        3. Schema query (verifies database is usable)
        
        Args:
            db_path: Path to the SQLite database file
        
        Returns:
            bool: True if database is valid, False if corrupted
        """
        import sqlite3
        
        # Level 1: Quick header check (< 1ms)
        try:
            with open(db_path, 'rb') as f:
                header = f.read(16)
            
            if not header.startswith(b'SQLite format 3'):
                return False
        except Exception:
            return False
        
        # Level 2: Connection test (~ 10ms)
        try:
            conn = sqlite3.connect(db_path)
        except sqlite3.DatabaseError:
            return False
        
        # Level 3: Schema validation (~ 50ms)
        try:
            cursor = conn.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1;")
            cursor.fetchone()
            conn.close()
            return True
        except sqlite3.DatabaseError:
            try:
                conn.close()
            except:
                pass
            return False
