"""
OpenEye OMEGA conformer generation backend.

Uses OpenEye's OMEGA executable for high-quality conformer generation.
File-based workflow: SMILES → .smi → OMEGA → .oeb
"""

from typing import Dict, Any, Iterator
from pathlib import Path
import os
import subprocess
import threading
import time
import math
import pandas as pd

from rdkit import Chem

from .base import ConformerBackend
from ...utils.constants import HAS_OPENEYE, oechem
from ...utils.mols.convert import oe_to_rdkit


class OpenEyeBackend(ConformerBackend):
    """
    OpenEye OMEGA conformer generation backend.

    Features:
    - High-quality conformer generation
    - GPU acceleration support
    - Multiple generation modes (classic, macrocycle, etc.)
    - File-based workflow
    """
    __backend_name__ = 'OEO'

    # Backend metadata
    description = "Generate conformers using OpenEye OMEGA (requires license)"
    required_params = ['SMILES_column']
    optional_params = [
        'names_column', 'max_confs', 'rms_threshold', 'dropna', 'timeout',
        'mode', 'use_gpu', 'mpi_np', 'strict',
        'flipper', 'flipper_warts', 'flipper_maxcenters',
        'oeomega_path', 'convert_to_rdkit',
    ]

    def __init__(self, params, logger=None, context=None):
        """Initialize OpenEye backend."""
        super().__init__(params, logger, context)

        # Validate OpenEye availability
        if not HAS_OPENEYE:
            raise ImportError(
                "OpenEye Python toolkit not available. "
                "Install openeye-toolkits to use OpenEye backend."
            )

        # Validate OMEGA executable
        self._validate_omega_installation()

        # Track if generation has been run
        self._generated = False
        self._stop_monitoring = False
        self._expected_total = 0

        self.log(
            f"OpenEye OMEGA initialized: mode={params.mode}, "
            f"max_confs={params.max_confs}, gpu={params.use_gpu}, timeout={params.timeout}s"
        )

    def _validate_omega_installation(self):
        """Check that OMEGA executable is available."""
        import shutil

        omega_path = self.params.oeomega_path

        if not omega_path or not (Path(omega_path).exists() or shutil.which(omega_path)):
            raise RuntimeError(
                f"OMEGA executable not found at '{omega_path}'. "
                "Please set oeomega_path parameter to valid OMEGA location."
            )

    @property
    def _confs_dir(self) -> Path:
        """Conformers output directory."""
        if self.context and hasattr(self.context, 'output_dir'):
            return Path(self.context.output_dir) / "conformers"
        return Path("conformers")

    @property
    def _input_file(self) -> str:
        """Path to input SMILES file."""
        return str(self._confs_dir / "OEOMEGA_smiles.smi")

    @property
    def _output_file(self) -> str:
        """Path to output conformers file."""
        return str(self._confs_dir / "OEOMEGA_conformers.oeb")

    @property
    def _report_file(self) -> str:
        """Path to OMEGA report file."""
        return str(self._confs_dir / "OEOMEGA_rpt.csv")

    @property
    def _path_prefix(self) -> str:
        """Path prefix for OMEGA command (without file extension)."""
        return str(self._confs_dir / "OEOMEGA")

    def generate_conformers(self, smiles_list: list[str], names_list: list[str]) -> Dict[str, Dict[str, Any]]:
        """
        Generate conformers using OpenEye OMEGA with progress monitoring.

        Args:
            smiles_list: List of SMILES strings
            names_list: List of molecule identifiers

        Returns:
            Dict mapping names to basic generation results
        """
        self.log(f"Starting OMEGA generation for {len(smiles_list)} SMILES")

        # Store expected total for progress monitoring
        self._expected_total = len(smiles_list)

        # Write SMILES to input file
        self._write_smiles_file(smiles_list, names_list)

        # Run OMEGA executable with progress monitoring
        self._execute_omega()

        # Reorder output to match input order (OMEGA may process asynchronously)
        self._sort_output()

        self._generated = True

        # Get successful names and report
        success_mols = self.get_successful_names()

        # Return basic results dict
        results = {}
        for name in names_list:
            results[name] = {
                'success': name in success_mols,
                'n_conformers': 0,
                'status': '' if name in success_mols else 'Generation failed',
            }

        success_count = len(success_mols)
        self.log(f"OMEGA generation complete: {success_count}/{len(smiles_list)} succeeded")

        return results

    def _write_smiles_file(self, smiles_list: list[str], names_list: list[str]):
        """Write SMILES to input file."""
        Path(self._input_file).parent.mkdir(parents=True, exist_ok=True)

        with open(self._input_file, 'w') as f:
            for smiles, name in zip(smiles_list, names_list):
                if smiles and name:  # Basic validation
                    f.write(f"{smiles} {name}\n")

    def _sort_output(self):
        """
        Reorder OEB output file to match the input SMILES order.

        OMEGA may process molecules out of order (e.g. with MPI parallelism),
        so the output OEB is not guaranteed to match the input ordering.
        This method reorders the OEB so that molecules appear in the same
        order as the input SMILES file, which downstream actors rely on.

        Uses OEMolDatabase.Order() + Save() to reorder in-place, preserving
        OMEGA's compact internal encoding (no re-serialization overhead).
        """
        if not Path(self._output_file).exists():
            return

        db = oechem.OEMolDatabase(self._output_file)
        n_mols = db.NumMols()
        if n_mols == 0:
            return

        # Build name → database index mapping from OEB (title-only, fast)
        name_to_db_idx = {}
        for i in range(n_mols):
            name_to_db_idx[db.GetTitle(i)] = i

        # Read input order from SMILES file
        input_names = []
        with open(self._input_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    input_names.append(parts[-1])

        # Build reordered index vector (input order → db indices)
        indices = oechem.OEUIntVector()
        for name in input_names:
            if name in name_to_db_idx:
                indices.append(name_to_db_idx[name])

        # Check if already in order
        already_sorted = all(
            indices[i] == i for i in range(len(indices))
        ) if len(indices) == n_mols else False

        if already_sorted:
            self.log("OEB output already in input order, skipping sort")
            return

        # Reorder database and save to temp file, then replace original
        db.Order(indices)
        sorted_path = str(Path(self._output_file).parent / '_sorted_temp.oeb')
        db.Save(sorted_path)
        os.replace(sorted_path, self._output_file)
        self.log(f"Sorted OEB to input order ({len(indices)} molecules)")

    def get_successful_names(self) -> list[str]:
        """Read names of successfully generated molecules from output file (matches original)."""
        if not Path(self._output_file).exists():
            self.log("Output file does not exist", level='WARNING')
            return []

        if Path(self._output_file).stat().st_size == 0:
            self.log("Output file is empty", level='WARNING')
            return []

        # Read molecule names from OEB file
        try:
            ifs = oechem.oemolistream(self._output_file)
            if not ifs.IsValid():
                self.log(f"Unable to open output file: {self._output_file}", level='ERROR')
                return []

            names = []
            for oemol in ifs.GetOEMols():
                mol = oechem.OEMol(oemol)
                names.append(mol.GetTitle())

            return names

        except Exception as e:
            self.log(f"Failed to read output file: {e}", level='ERROR')
            return []

    def get_report_dataframe(self) -> pd.DataFrame:
        """Get OMEGA report as DataFrame.

        OMEGA native format: Molecule (SMILES), Title (name), Rotors, Conformers, ElapsedTime(s), Status
        """
        if not Path(self._report_file).exists():
            self.log("OMEGA report file not found", level='WARNING')
            return pd.DataFrame(columns=['Molecule', 'Title', 'Rotors', 'Conformers', 'ElapsedTime(s)', 'Status'])

        try:
            return pd.read_csv(self._report_file)
        except Exception as e:
            self.log(f"Failed to read OMEGA report: {e}", level='ERROR')
            return pd.DataFrame(columns=['Molecule', 'Title', 'Rotors', 'Conformers', 'ElapsedTime(s)', 'Status'])

    def _execute_omega(self):
        """Execute OMEGA conformer generation with progress monitoring."""
        # Start progress monitoring thread
        self._stop_monitoring = False
        progress_thread = threading.Thread(
            target=self._monitor_progress,
            daemon=True  # Dies when main thread exits
        )
        progress_thread.start()

        # Build OMEGA command
        cmd_args = [
            self.params.oeomega_path,
            self.params.mode,
            "-in", self._input_file,
            "-out", self._output_file,
            "-verbose", "false",
            "-maxconfs", str(self.params.max_confs),
            "-rms", str(self.params.rms_threshold),
            "-mpi_np", str(self.params._resolved_mpi_np),
            "-useGPU", str(self.params.use_gpu).lower(),
            "-prefix", self._path_prefix,
            "-flipper", str(self.params.flipper).lower(),
            "-flipper_warts", str(self.params.flipper_warts).lower(),
            "-flipper_maxcenters", str(self.params.flipper_maxcenters),
            "-strict", str(self.params.strict).lower(),
        ]

        # Execute
        try:
            result = subprocess.run(
                cmd_args,
                check=False,
                stdout=False,
                stderr=False,
                timeout=self.params.timeout  # 1 hour timeout
            )

            # Signal progress monitoring to stop
            self._stop_monitoring = True

            if result.returncode != 0:
                error_msg = f"OMEGA execution failed with return code {result.returncode}"
                self.log(error_msg, level='ERROR')
                raise RuntimeError(error_msg)

            self.log(f"OMEGA execution successful. Output: {self._output_file}")

        except subprocess.TimeoutExpired:
            self._stop_monitoring = True
            raise RuntimeError("OMEGA execution timeout (1 hour)")
        except Exception as e:
            self._stop_monitoring = True
            raise RuntimeError(f"OMEGA execution failed: {str(e)}")

    def _monitor_progress(self):
        """Monitor OMEGA progress by reading report file line count."""
        poll_interval = 0.1
        base_log_interval = 5
        start_time = time.time()
        
        last_count = 0
        last_log_time = start_time # don't log immediately
        has_logged = False
        
        while not self._stop_monitoring:
            try:
                time.sleep(poll_interval)
                
                if not Path(self._report_file).exists():
                    continue
                
                # Count lines in report file (skip header)
                with open(self._report_file, 'r') as f:
                    current_count = sum(1 for line in f) - 1
                
                # Skip if no new progress
                if current_count <= last_count:
                    continue
                
                now = time.time()
                elapsed = now - start_time
                
                # Calculate metrics
                rate = current_count / elapsed if elapsed > 0 else 0
                
                # Determine if we should log
                should_log = False
                log_interval = base_log_interval
                
                if self._expected_total > 0:
                    percent = (current_count / self._expected_total) * 100
                    eta = (self._expected_total - current_count) / rate if rate > 0 else 0
                    
                    # Adaptive log interval based on ETA
                    if eta > 0:
                        log_interval = min(max(base_log_interval, int(10 ** (math.log10(max(eta, 1)) - 1))), 1000)
                    
                    should_log = (now - last_log_time) >= log_interval
                    
                    if should_log:
                        self.log(
                            f"Progress: {current_count}/{self._expected_total} ({percent:.1f}%) | "
                            f"Rate: {rate:.2f} mol/s | ETA: {self._format_duration(eta)}"
                        )
                else:
                    should_log = (now - last_log_time) >= log_interval
                    
                    if should_log:
                        self.log(
                            f"Progress: {current_count} molecules processed | "
                            f"Rate: {rate:.2f} mol/s | Elapsed: {self._format_duration(elapsed)}"
                        )
                
                if should_log:
                    last_log_time = now
                    last_count = current_count
                    has_logged = True
                
            except Exception as e:
                self.log(f"Progress monitoring error: {e}", level='DEBUG')
        
        # Final report - only if we had progress and either:
        # 1. Never logged (short process), OR
        # 2. Logged before but final count differs from last logged count
        try:
            if Path(self._report_file).exists():
                with open(self._report_file, 'r') as f:
                    final_count = sum(1 for line in f) - 1
                
                # Only print final report if there's something to report and we haven't already reported completion
                if final_count > 0 and (not has_logged or final_count > last_count):
                    elapsed = time.time() - start_time
                    rate = final_count / elapsed if elapsed > 0 else 0
                    
                    if self._expected_total > 0:
                        percent = (final_count / self._expected_total) * 100
                        self.log(
                            f"Progress: {final_count}/{self._expected_total} ({percent:.1f}%) | "
                            f"Rate: {rate:.2f} mol/s | Elapsed: {self._format_duration(elapsed)}"
                        )
                    else:
                        self.log(
                            f"Completed: {final_count} molecules | "
                            f"Rate: {rate:.2f} mol/s | Elapsed: {self._format_duration(elapsed)}"
                        )
        except Exception:
            pass

    def _format_duration(self, seconds: float) -> str:
        """Format duration in human-readable format."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            return f"{seconds / 60:.1f}m"
        else:
            hours = seconds // 3600
            minutes = (seconds % 3600) / 60
            return f"{int(hours)}h{int(minutes)}m"

    def get_endpoint(self) -> str:
        """
        Return path to conformer output file.

        For file-based backends, endpoint is the file path.
        """
        return self._output_file

    def extract_molecules(self) -> Iterator[Chem.Mol]:
        """
        Extract molecules from OEB file, optionally converting to RDKit.

        Yields:
            RDKit molecules (if convert_to_rdkit=True) or OpenEye molecules (if False)
        """
        if not Path(self._output_file).exists():
            self.log("Output file not found", level='ERROR')
            return

        # Read OEB file
        try:
            ifs = oechem.oemolistream(self._output_file)

            for oemol in ifs.GetOEMols():
                if self.params.convert_to_rdkit:
                    # Convert OpenEye mol to RDKit
                    try:
                        rdkit_mol = oe_to_rdkit(oemol)
                        if rdkit_mol:
                            yield rdkit_mol
                    except Exception as e:
                        self.log(f"Failed to convert molecule: {e}", level='WARNING')
                        continue
                else:
                    # Yield OpenEye molecule directly (no conversion)
                    yield oechem.OEMol(oemol)

        except Exception as e:
            self.log(f"Failed to read output file: {e}", level='ERROR')