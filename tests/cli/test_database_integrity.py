"""
Tests for CLI database handling and auto-download functionality.

These tests validate:
1. CLI works with missing database (auto-download enabled by default)
2. Error handling for corrupted databases
3. Database integrity checking
4. Auto-recovery from corruption

Test Categories:
- Unit tests (fast): Test individual components with mocking
- Integration tests (slow): Test full CLI with real downloads
"""

import pytest
import os
import sqlite3
import tempfile
import shutil
from unittest.mock import patch, MagicMock
from pathlib import Path

from molforge.cli.main import main
from molforge.actors.params.source import ChEMBLSourceParams


@pytest.mark.database
class TestDatabaseAutoDownload:
    """Test automatic database download functionality"""

    def test_auto_download_enabled_by_default(self):
        """Verify auto_download is enabled by default in params"""
        params = ChEMBLSourceParams(backend='sql')
        assert params.auto_download is True, "auto_download should be True by default"

    @pytest.mark.slow
    @pytest.mark.integration
    def test_cli_with_missing_database(self, tmp_path, capsys, monkeypatch):
        """Test CLI handles missing database with auto-download (SLOW: real download)"""
        # This test actually downloads the full ChEMBL database
        # Skip in fast test runs with: pytest -m "not slow"
        
        # Create temporary directory for database
        db_dir = tmp_path / "chembl_test"
        db_dir.mkdir()
        db_path = db_dir / "36.db"
        
        # Mock sys.argv to simulate CLI call
        test_args = ['molforge', 'run', 'CHEMBL234', '-q']
        monkeypatch.setattr('sys.argv', test_args)
        
        # Mock ChEMBLSourceParams to use our temp directory
        original_init = ChEMBLSourceParams.__init__
        
        def mock_init(self, *args, **kwargs):
            kwargs['download_dir'] = str(db_dir)
            original_init(self, *args, **kwargs)
        
        with patch.object(ChEMBLSourceParams, '__init__', mock_init):
            # Run CLI - should trigger auto-download
            try:
                exit_code = main()
            except SystemExit as e:
                exit_code = e.code
            
            # Verify database was downloaded
            assert db_path.exists(), "Database should be downloaded"
            assert db_path.stat().st_size > 1_000_000, "Database should be substantial size"
            
            # Verify no corruption errors
            captured = capsys.readouterr()
            assert 'database disk image is malformed' not in captured.err.lower()
    
    def test_auto_download_disabled(self, tmp_path):
        """Test that auto_download=False raises error for missing database"""
        db_dir = tmp_path / "chembl_test"
        db_dir.mkdir()
        db_path = db_dir / "36.db"
        
        # Create params with auto_download disabled
        with pytest.raises(FileNotFoundError, match="ChEMBL database not found"):
            params = ChEMBLSourceParams(
                backend='sql',
                db_path=str(db_path),
                auto_download=False
            )
            # Trigger validation
            params._validate_params()
    
    @staticmethod
    def _create_minimal_db(path):
        """Create a minimal valid SQLite database for testing"""
        conn = sqlite3.connect(path)
        cursor = conn.cursor()
        
        # Create minimal schema
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS activities (
                activity_id INTEGER PRIMARY KEY,
                molregno INTEGER,
                standard_value REAL,
                standard_type TEXT,
                standard_units TEXT
            )
        ''')
        
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecule_dictionary (
                molregno INTEGER PRIMARY KEY,
                chembl_id TEXT,
                pref_name TEXT
            )
        ''')
        
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS compound_structures (
                molregno INTEGER PRIMARY KEY,
                canonical_smiles TEXT
            )
        ''')
        
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS target_dictionary (
                tid INTEGER PRIMARY KEY,
                chembl_id TEXT,
                pref_name TEXT
            )
        ''')
        
        conn.commit()
        conn.close()


@pytest.mark.database
class TestDatabaseIntegrity:
    """Test database corruption detection"""
    
    def test_detect_corrupted_database(self, tmp_path):
        """Test detection of corrupted database file"""
        # Create a corrupted database file
        corrupted_db = tmp_path / "corrupted.db"
        with open(corrupted_db, 'wb') as f:
            f.write(b"This is not a valid SQLite database file")
        
        # Try to open it with sqlite3
        with pytest.raises(sqlite3.DatabaseError):
            conn = sqlite3.connect(str(corrupted_db))
            conn.execute("PRAGMA integrity_check;")
            conn.close()
    
    def test_integrity_check_on_valid_database(self, tmp_path):
        """Test integrity check passes on valid database"""
        # Create a valid database
        valid_db = tmp_path / "valid.db"
        conn = sqlite3.connect(str(valid_db))
        conn.execute("CREATE TABLE test (id INTEGER PRIMARY KEY)")
        conn.commit()
        
        # Run integrity check
        cursor = conn.execute("PRAGMA integrity_check;")
        result = cursor.fetchone()[0]
        conn.close()
        
        assert result == "ok", f"Integrity check should return 'ok', got: {result}"
    
    def test_quick_header_check(self, tmp_path):
        """Test quick SQLite header validation"""
        # Valid database should start with "SQLite format 3\x00"
        valid_db = tmp_path / "valid.db"
        conn = sqlite3.connect(str(valid_db))
        conn.execute("CREATE TABLE test (id INTEGER)")
        conn.commit()
        conn.close()
        
        # Check header
        with open(valid_db, 'rb') as f:
            header = f.read(16)
        
        assert header.startswith(b'SQLite format 3'), "Valid database should have SQLite header"
        
        # Create corrupted file
        corrupted_db = tmp_path / "corrupted.db"
        with open(corrupted_db, 'wb') as f:
            f.write(b"Invalid header content")
        
        # Check corrupted header
        with open(corrupted_db, 'rb') as f:
            bad_header = f.read(16)
        
        assert not bad_header.startswith(b'SQLite format 3'), "Corrupted file should not have SQLite header"
    
    def test_check_database_integrity_method(self, tmp_path):
        """Test ChEMBLSourceParams._check_database_integrity() method"""
        params = ChEMBLSourceParams(backend='sql')
        
        # Test with valid database
        valid_db = tmp_path / "valid.db"
        TestDatabaseAutoDownload._create_minimal_db(str(valid_db))
        assert params._check_database_integrity(str(valid_db)) is True
        
        # Test with corrupted database (bad header)
        corrupted_db = tmp_path / "corrupted.db"
        with open(corrupted_db, 'wb') as f:
            f.write(b"CORRUPTED DATA NOT SQLITE")
        assert params._check_database_integrity(str(corrupted_db)) is False
        
        # Test with non-existent file
        missing_db = tmp_path / "missing.db"
        assert params._check_database_integrity(str(missing_db)) is False


@pytest.mark.database
class TestCorruptionRecovery:
    """Test automatic recovery from corrupted databases"""
    
    def test_corrupted_db_with_auto_download_disabled(self, tmp_path):
        """Test that corrupted DB with auto_download=False raises error"""
        # Create corrupted database
        corrupted_db = tmp_path / "corrupted.db"
        with open(corrupted_db, 'wb') as f:
            f.write(b"NOT A VALID SQLITE DATABASE FILE")
        
        # Should raise RuntimeError about corruption
        with pytest.raises(RuntimeError, match="corrupted"):
            params = ChEMBLSourceParams(
                backend='sql',
                db_path=str(corrupted_db),
                auto_download=False
            )
            params._validate_params()
    
    def test_corrupted_db_auto_recovery_with_mock(self, tmp_path, capsys):
        """Test that corrupted DB triggers auto-download when enabled (fast with mock)"""
        db_dir = tmp_path / "chembl"
        db_dir.mkdir()
        corrupted_db = db_dir / "36.db"
        
        # Create corrupted database
        with open(corrupted_db, 'wb') as f:
            f.write(b"CORRUPTED DATABASE FILE WITH RANDOM BYTES")
        
        # Mock the download function to create valid DB instead
        with patch('molforge.actors.params.source.ensure_chembl_db') as mock_download:
            mock_download.side_effect = lambda path, version: TestDatabaseAutoDownload._create_minimal_db(path)
            
            # Create params with auto_download enabled
            params = ChEMBLSourceParams(
                backend='sql',
                download_dir=str(db_dir),
                auto_download=True
            )
            
            # Validate - should detect corruption and re-download
            params._validate_params()
            
            # Verify download was called
            assert mock_download.called, "ensure_chembl_db should be called for corrupted database"
            
            # Verify warning messages
            captured = capsys.readouterr()
            assert "corrupted" in captured.out.lower() or "corrupted" in captured.err.lower()
            
            # Verify database is now valid
            assert params._check_database_integrity(params.db_path) is True
    
    @pytest.mark.slow
    @pytest.mark.integration
    def test_corrupted_db_full_recovery(self, tmp_path, capsys):
        """Test full corruption recovery with real download (SLOW)"""
        db_dir = tmp_path / "chembl"
        db_dir.mkdir()
        corrupted_db = db_dir / "36.db"
        
        # Create corrupted database with random bytes
        with open(corrupted_db, 'wb') as f:
            # Write random data that looks like it could be a file
            import random
            f.write(bytes([random.randint(0, 255) for _ in range(10000)]))
        
        # Verify it's corrupted
        params_check = ChEMBLSourceParams(backend='sql')
        assert params_check._check_database_integrity(str(corrupted_db)) is False
        
        # Create params with auto_download - should recover
        params = ChEMBLSourceParams(
            backend='sql',
            download_dir=str(db_dir),
            auto_download=True,
            version=36
        )
        
        # Validate - this will download real database
        params._validate_params()
        
        # Verify recovery
        assert params._check_database_integrity(params.db_path) is True
        assert Path(params.db_path).stat().st_size > 1_000_000
        
        # Verify messages
        captured = capsys.readouterr()
        assert "corrupted" in captured.out.lower()
        assert "downloading" in captured.out.lower() or "auto-downloading" in captured.out.lower()


@pytest.mark.database
class TestDatabasePath:
    """Test database path handling"""
    
    def test_default_download_directory(self):
        """Test default download directory is ./data/chembl"""
        params = ChEMBLSourceParams(backend='sql')
        assert params.download_dir == "./data/chembl"
    
    def test_custom_download_directory(self, tmp_path):
        """Test custom download directory is respected"""
        custom_dir = str(tmp_path / "custom_chembl")
        
        with patch('molforge.actors.params.source.ensure_chembl_db') as mock_download:
            # Mock needs to create actual valid database for integrity check
            mock_download.side_effect = lambda path, version: TestDatabaseAutoDownload._create_minimal_db(path)
            
            params = ChEMBLSourceParams(
                backend='sql',
                download_dir=custom_dir,
                auto_download=True
            )
            params._validate_params()
        
        assert params.download_dir == custom_dir
        assert os.path.exists(custom_dir), "Custom directory should be created"
    
    def test_db_path_override(self, tmp_path):
        """Test explicit db_path overrides default location"""
        custom_db = tmp_path / "my_chembl.db"
        
        # Create database file
        TestDatabaseAutoDownload._create_minimal_db(str(custom_db))
        
        params = ChEMBLSourceParams(
            backend='sql',
            db_path=str(custom_db),
            auto_download=False
        )
        params._validate_params()
        
        assert params.db_path == str(custom_db)
        assert params.version == 'manual path'


@pytest.fixture
def mock_chembl_db(tmp_path):
    """Fixture providing a mock ChEMBL database for testing"""
    db_path = tmp_path / "36.db"
    TestDatabaseAutoDownload._create_minimal_db(str(db_path))
    return str(db_path)


@pytest.mark.database
def test_cli_run_with_valid_database(mock_chembl_db, tmp_path, monkeypatch, capsys):
    """Integration test: CLI run with valid database"""
    # Extract directory from db path
    db_dir = os.path.dirname(mock_chembl_db)
    
    test_args = ['molforge', 'run', 'CHEMBL234', '-q']
    monkeypatch.setattr('sys.argv', test_args)
    
    # Mock ChEMBLSourceParams to use our test database
    original_init = ChEMBLSourceParams.__init__
    
    def mock_init(self, *args, **kwargs):
        kwargs['download_dir'] = db_dir
        kwargs['db_path'] = mock_chembl_db
        original_init(self, *args, **kwargs)
    
    with patch.object(ChEMBLSourceParams, '__init__', mock_init):
        try:
            exit_code = main()
        except SystemExit as e:
            exit_code = e.code
        
        # Should run without database errors (may fail on missing data, which is OK)
        captured = capsys.readouterr()
        assert 'database disk image is malformed' not in captured.err.lower()
