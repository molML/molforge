"""
OpenEye OMEGA conformer generation backend.

Uses OpenEye's OMEGA executable for high-quality conformer generation.
File-based workflow: SMILES → .smi → OMEGA → .oeb
"""

from typing import Dict, Any, Iterator
from pathlib import Path
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
    required_params = ['SMILES_column', 'output_file']
    optional_params = [
        'names_column', 'max_confs', 'rms_threshold', 'mode',
        'verbose', 'dropna'
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
            f"max_confs={params.max_confs}, gpu={params.use_gpu}"
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
    def _input_file(self) -> str:
        """Path to input SMILES file."""
        if self.context and hasattr(self.context, 'output_dir'):
            return str(Path(self.context.output_dir) / "OEOMEGA_smiles.smi")
        return "OEOMEGA_smiles.smi"

    @property
    def _output_file(self) -> str:
        """Path to output conformers file."""
        if self.context and hasattr(self.context, 'output_dir'):
            return str(Path(self.context.output_dir) / "OEOMEGA_conformers.oeb")
        return "OEOMEGA_conformers.oeb"

    @property
    def _report_file(self) -> str:
        """Path to OMEGA report file."""
        if self.context and hasattr(self.context, 'output_dir'):
            return str(Path(self.context.output_dir) / "OEOMEGA_rpt.csv")
        return "OEOMEGA_rpt.csv"

    @property
    def _path_prefix(self) -> str:
        """Path prefix for OMEGA command (without file extension)."""
        if self.context and hasattr(self.context, 'output_dir'):
            base = str(Path(self.context.output_dir) / "OEOMEGA")
            return base.replace(".OEOMEGA", "")
        return "OEOMEGA"

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

        self._generated = True

        # Get successful names and report
        success_mols = self.get_successful_names()

        # Return basic results dict
        results = {}
        for name in names_list:
            results[name] = {
                'success': name in success_mols,
                'n_conformers': 0,  # Will be filled from report
                'status': 'SUCCESS' if name in success_mols else 'FAILED',
                'error': None if name in success_mols else 'Generation failed'
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
        """Get OMEGA report as DataFrame (matches original)."""
        if not Path(self._report_file).exists():
            self.log("OMEGA report file not found", level='WARNING')
            return pd.DataFrame(columns=['Molecule', 'SMILES', 'Rotors', 'Conformers', 'ElapsedTime(s)', 'Status'])

        try:
            return pd.read_csv(self._report_file)
        except Exception as e:
            self.log(f"Failed to read OMEGA report: {e}", level='ERROR')
            return pd.DataFrame(columns=['Molecule', 'SMILES', 'Rotors', 'Conformers', 'ElapsedTime(s)', 'Status'])

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
            "-mpi_np", str(self.params.mpi_np),
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
                timeout=3600  # 1 hour timeout
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
        update_interval = 10  # Check every 10 seconds
        last_count = 0
        start_time = time.time()

        while not self._stop_monitoring:
            try:
                time.sleep(update_interval)

                # Check if report file exists yet
                if not Path(self._report_file).exists():
                    continue

                # Count lines in report file (skip header)
                with open(self._report_file, 'r') as f:
                    current_count = sum(1 for line in f) - 1

                # Only report if progress was made
                if current_count <= last_count:
                    continue

                # Calculate metrics
                elapsed = time.time() - start_time
                rate = current_count / elapsed if elapsed > 0 else 0

                if self._expected_total > 0:
                    percent = (current_count / self._expected_total) * 100
                    eta = (self._expected_total - current_count) / rate if rate > 0 else 0
                    self.log(
                        f"Progress: {current_count}/{self._expected_total} ({percent:.1f}%) | "
                        f"Rate: {rate:.2f} mol/s | ETA: {self._format_duration(eta)}"
                    )
                else:
                    self.log(
                        f"Progress: {current_count} molecules processed | "
                        f"Rate: {rate:.2f} mol/s | Elapsed: {self._format_duration(elapsed)}"
                    )

                last_count = current_count

                # Adjust update interval based on ETA to avoid spam
                if eta > 1:
                    update_interval = min(max(10, int(10 ** (math.log10(eta) - 1))), 1000)

            except Exception as e:
                self.log(f"Progress monitoring error: {e}", level='DEBUG')
                time.sleep(30)  # Wait longer if error occurred

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

    def _parse_omega_report(self, names_list: list[str]) -> Dict[str, Dict[str, Any]]:
        """
        Parse OMEGA report CSV to get generation results.

        Args:
            names_list: List of molecule names that were processed

        Returns:
            Dict mapping molecule names to generation results
        """
        results = {}

        if not Path(self._report_file).exists():
            self.log("OMEGA report file not found", level='WARNING')
            # Return failure for all molecules
            for name in names_list:
                results[name] = {
                    'success': False,
                    'n_conformers': 0,
                    'status': 'REPORT_NOT_FOUND',
                    'error': 'Report file not found'
                }
            return results

        # Parse CSV report
        try:
            df = pd.read_csv(self._report_file)

            for name in names_list:
                # Find this molecule in report
                row = df[df['Molecule'] == name]

                if len(row) == 0:
                    results[name] = {
                        'success': False,
                        'n_conformers': 0,
                        'status': 'NOT_IN_REPORT',
                        'error': 'Not in OMEGA report'
                    }
                else:
                    row = row.iloc[0]
                    status = row.get('Status', 'UNKNOWN')
                    n_confs = int(row.get('conformers', 0))
                    success = (n_confs > 0 and status == 'SUCCESS')

                    results[name] = {
                        'success': success,
                        'n_conformers': n_confs,
                        'status': status,
                        'error': None if success else f'OMEGA status: {status}'
                    }

        except Exception as e:
            self.log(f"Failed to parse OMEGA report: {e}", level='ERROR')
            # Return failure for all
            for name in names_list:
                results[name] = {
                    'success': False,
                    'n_conformers': 0,
                    'status': 'PARSE_ERROR',
                    'error': f'Report parse error: {str(e)}'
                }

        return results

    def get_endpoint(self) -> str:
        """
        Return path to conformer output file.

        For file-based backends, endpoint is the file path.
        """
        return self._output_file

    def extract_molecules(self) -> Iterator[Chem.Mol]:
        """
        Extract molecules from OEB file.

        Yields:
            RDKit molecules with conformers
        """
        if not Path(self._output_file).exists():
            self.log("Output file not found", level='ERROR')
            return

        # Read OEB file and convert to RDKit
        try:
            ifs = oechem.oemolistream(self._output_file)

            for oemol in ifs.GetOEMols():
                # Convert OpenEye mol to RDKit
                try:
                    rdkit_mol = oe_to_rdkit(oemol)
                    if rdkit_mol:
                        yield rdkit_mol
                except Exception as e:
                    self.log(f"Failed to convert molecule: {e}", level='WARNING')
                    continue

        except Exception as e:
            self.log(f"Failed to read output file: {e}", level='ERROR')