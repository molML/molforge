"""
Train/validation/test data splitter plugin.

Partitions a curated dataset into train, validation, and test sets along two
orthogonal axes.

``unit`` — the atom of assignment
    ``"scaffold"``  group molecules by Bemis-Murcko scaffold and assign whole
                    clusters together, so no scaffold spans more than one split.
                    Acyclic molecules and rows without a computable scaffold are
                    always assigned to train.
    ``"molecule"``  treat each molecule as its own unit.  Scaffolds may then
                    span splits; this is the interpolation baseline.

``method`` — how units are ordered before filling test → val → train
    ``"isolation"``  rank units by size-weighted mean ECFP4 Tanimoto distance to
                     the rest of the dataset and fill the most isolated first,
                     yielding a chemically out-of-distribution test set for
                     prospective evaluation.  Deterministic.
    ``"random"``     order units by a seeded permutation.  Reproducible from
                     ``seed``.

The four combinations are scaffold/isolation (the default), scaffold/random,
molecule/random, and molecule/isolation.

This plugin requires scaffold assignment from the upstream ``scaffold``
actor (``scaffold_smiles`` column by default).

Design notes
------------
All fingerprints are represented as float32 NumPy arrays throughout.  This
keeps every joblib worker boundary free of RDKit C++ objects and ensures
complete pickling safety.  Pairwise Tanimoto matrices are computed via
vectorised NumPy matrix multiplication rather than row-by-row loops.

Output columns
--------------
split
    ``'train'``, ``'val'``, or ``'test'`` for every row.

Report card
-----------
A JSON report card is written alongside the pipeline output containing split
composition, scaffold statistics, structural separation (ECFP4 nearest-
neighbour Tanimoto distances, all pairwise split combinations), and optional activity
distribution diagnostics.

Example
-------
::

    from molforge import MolForge, ForgeParams
    from molforge.actor_plugins.split import SplitParams

    params = ForgeParams(
        steps=["source", "chembl", "curate", "scaffold", "split"],
        plugin_params={
            "scaffold": ComputeScaffoldsParams(...),
            "split":    SplitParams(test_ratio=0.10, val_ratio=0.10),
        },
        ...
    )
    forge  = MolForge(params)
    result = forge.forge("CHEMBL203")
"""

import json
import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from molforge.actors.base import BaseActor
from molforge.actors.params.base import BaseParams
from molforge.actors.protocol import ActorOutput
from molforge.utils.constants import DEFAULT_N_JOBS, DEFAULT_MP_THRESHOLD


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

@dataclass
class SplitParams(BaseParams):
    """Configuration for two-axis train/val/test splitting.

    See the per-field docstrings below for details on each parameter.
    """

    unit: str = "scaffold"
    """Atom of assignment.

    - ``"scaffold"`` — group molecules by ``scaffold_column`` (Bemis-Murcko) and
      assign whole clusters together, so no scaffold spans more than one split.
      Acyclic molecules and rows without a computable scaffold are always
      assigned to train.
    - ``"molecule"`` — treat each molecule as its own unit; scaffolds may span
      splits (the interpolation baseline)."""

    method: str = "isolation"
    """Ordering used to fill test → val → train.

    - ``"isolation"`` — rank units by size-weighted mean ECFP4 Tanimoto distance
      to the rest of the dataset and fill the most isolated first.
      Deterministic; ignores ``seed``.
    - ``"random"`` — order units by a seeded permutation, reproducible from
      ``seed``."""

    SMILES_column: str = "curated_smiles"
    """Column containing curated SMILES strings used for fingerprint
    computation.  Should be the output of the ``curate`` actor.  Also the
    representative structure fingerprinted per unit when ``unit='molecule'``."""

    split_column: str = "split"
    """Name of the output column that receives the ``train``/``val``/``test``
    label for every row."""

    train_label: str = "train"
    """Value written to ``split_column`` for training-set rows."""

    val_label: str = "val"
    """Value written to ``split_column`` for validation-set rows."""

    test_label: str = "test"
    """Value written to ``split_column`` for test-set rows."""

    scaffold_column: str = "scaffold_smiles"
    """Column containing Bemis-Murcko scaffold SMILES from an upstream
    ``scaffold`` actor.  This column is required.  It defines the clusters when
    ``unit='scaffold'`` and provides the scaffold statistics in the report card
    for both units."""

    test_ratio: float = 0.10
    """Fraction of molecules to assign to the test set."""

    val_ratio: float = 0.10
    """Fraction of molecules to assign to the validation set.  The training set
    receives the remaining ``1 − test_ratio − val_ratio``."""

    max_units_for_isolation: int = 10_000
    """Compute bound for ``method='isolation'``: the maximum number of units for
    which the full pairwise Tanimoto distance matrix is built to rank units by
    structural isolation.  The matrix is O(n²) in both time and memory:

    - ≤  1,000 units → ~  1 s,   ~  16 MB
    - ≤  5,000 units → ~ 10 s,   ~ 100 MB
    - ≤ 10,000 units → ~ 40 s,   ~ 400 MB
    - ≤ 20,000 units → ~200 s,   ~ 1.6 GB

    When the number of units exceeds this limit the run raises ``ValueError``.
    Raise this bound to permit the larger matrix, or set ``method='random'`` for
    an order-independent split."""

    low_n_threshold: int = 5_000
    """Advisory threshold on dataset size.  When the dataset has fewer rows than
    this, the actor logs a warning that the split — isolation-based in
    particular — may be unreliable or leave a set underpopulated.  Advisory
    only; it does not stop the run."""

    seed: int = 42
    """Random seed for the permutation used by ``method='random'``.  Ignored by
    ``method='isolation'``, which is deterministic."""

    activity_column: Optional[str] = None
    """Column name for the continuous activity label used in the report card's
    activity-distribution statistics and KS test.  When left as ``None`` it is
    resolved automatically from the upstream curator's declared endpoint (the
    standardized activity column, e.g. pIC50/pEC50).  Activity diagnostics are
    silently skipped when no such column can be resolved or found."""

    report_card_filename: str = "split_report.json"
    """Filename for the JSON report card saved in the pipeline output
    directory."""

    ecfp4_radius: int = 2
    """Morgan fingerprint radius for ECFP4 computation.  Used for both scaffold
    isolation scoring and molecule-level NN distance statistics."""

    ecfp4_n_bits: int = 2048
    """Fingerprint bit vector length.  Larger values reduce hash collisions at
    the cost of memory."""

    max_nn_dist_n: int = 10_000
    """Maximum dataset size (number of molecules) for which molecule-level
    ECFP4 nearest-neighbour Tanimoto distance statistics are computed.  The
    full pairwise distance matrix is O(n²); for 10,000 molecules this is
    ~400 MB and takes a few seconds.  Skipped with an INFO note in the report
    card when the dataset exceeds this limit."""

    def _validate_params(self) -> None:
        self._validate_policy('unit', self.unit, ['scaffold', 'molecule'])
        self._validate_policy('method', self.method, ['isolation', 'random'])
        if not self.SMILES_column:
            raise ValueError("SMILES_column cannot be empty.")
        if not self.split_column:
            raise ValueError("split_column cannot be empty.")
        if not self.train_label or not self.val_label or not self.test_label:
            raise ValueError("train_label, val_label, and test_label cannot be empty.")
        if len({self.train_label, self.val_label, self.test_label}) != 3:
            raise ValueError("train_label, val_label, and test_label must be distinct.")
        if not self.scaffold_column:
            raise ValueError("scaffold_column cannot be empty.")
        if not (0.0 < self.test_ratio < 1.0):
            raise ValueError(f"test_ratio must be in (0, 1), got {self.test_ratio}.")
        if not (0.0 < self.val_ratio < 1.0):
            raise ValueError(f"val_ratio must be in (0, 1), got {self.val_ratio}.")
        if self.test_ratio + self.val_ratio >= 1.0:
            raise ValueError(
                f"test_ratio + val_ratio must be < 1.0, "
                f"got {self.test_ratio + self.val_ratio:.2f}."

            )
        if self.max_units_for_isolation < 1:
            raise ValueError(
                f"max_units_for_isolation must be >= 1, "
                f"got {self.max_units_for_isolation}."
            )
        if self.low_n_threshold < 0:
            raise ValueError(
                f"low_n_threshold must be >= 0, got {self.low_n_threshold}."
            )
        if self.ecfp4_radius < 1:
            raise ValueError(
                f"ecfp4_radius must be >= 1, got {self.ecfp4_radius}."
            )
        if self.ecfp4_n_bits < 64:
            raise ValueError(
                f"ecfp4_n_bits must be >= 64, got {self.ecfp4_n_bits}."
            )
        if self.max_nn_dist_n < 1:
            raise ValueError(
                f"max_nn_dist_n must be >= 1, got {self.max_nn_dist_n}."
            )


# ---------------------------------------------------------------------------
# Actor
# ---------------------------------------------------------------------------

class Splitter(BaseActor):
    """Assign train/val/test labels to every row in the pipeline along two axes.

    Adds a single ``split`` column with values ``'train'``, ``'val'``, or
    ``'test'``.  Rows are grouped into units by ``unit`` and the units are
    ordered by ``method`` before filling test → val → train:

    - ``unit='scaffold'`` groups molecules by Bemis-Murcko scaffold and assigns
      whole clusters together, so no scaffold spans more than one split.
      Acyclic molecules (empty scaffold) and rows without a computable scaffold
      are always assigned to train.
    - ``unit='molecule'`` makes each molecule its own unit; scaffolds may span
      splits (the interpolation baseline).
    - ``method='isolation'`` assigns the most structurally isolated units — those
      whose ECFP4 representative fingerprint is most distant from the rest of the
      dataset (weighted-mean Tanimoto) — to test first, then val, then train.
      Deterministic.
    - ``method='random'`` orders units by a seeded permutation.

    A JSON report card is written to the pipeline output directory summarising
    split composition, scaffold statistics, structural separation, and
    (optionally) activity distribution diagnostics.

    Notes
    -----
    Under ``method='isolation'`` a run with more than
    ``max_units_for_isolation`` units raises ``ValueError`` rather than building
    the O(n²) distance matrix; raise the bound or use ``method='random'``.

    Molecule-level ECFP4 NN-distance statistics are skipped for datasets
    larger than ``max_nn_dist_n`` rows (default 10,000).
    """

    __step_name__   = 'split'
    __dependencies__ = ['scaffold']
    __param_class__ = SplitParams
    __terminal__    = True   # terminal step: expected to be the last step in the pipeline

    # Report-card figure palette (self-contained; no external style dependency).
    _TRAIN_COLOR   = '#f46a67'   # burnt sienna
    _VAL_COLOR     = '#769ff5'   # vista blue
    _TEST_COLOR    = '#c157ba'   # steel pink

    _HM_CMAP_COLORS = ["white", "#219492"]   # white → teal; used for NN heatmap
    _OVERLAP_COLOR = '#c0c0c0'   # grey: scaffold shared with train
    _HIST_ALPHA    = 0.5
    _FIG_SIZE      = (16, 8)
    _FIG_DPI       = 120
    _FONT_SMALL    = 8
    _FONT_MEDIUM   = 9
    _FONT_TITLE    = 13

    @property
    def required_columns(self) -> List[str]:
        return [self.SMILES_column, self.scaffold_column]

    @property
    def output_columns(self) -> List[str]:
        return [self.split_column]

    @property
    def _split_colors(self) -> Dict[str, str]:
        return {
            self.train_label: self._TRAIN_COLOR,
            self.val_label: self._VAL_COLOR,
            self.test_label: self._TEST_COLOR,
        }

    # ------------------------------------------------------------------
    # ECFP4 / scaffold helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _ecfp4_batch(
        smiles_list: List[str],
        radius: int,
        n_bits: int,
    ) -> List[Optional[np.ndarray]]:
        """Compute ECFP4 fingerprints for a batch of SMILES."""
        from rdkit import Chem
        from rdkit.Chem import DataStructs, rdFingerprintGenerator

        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        results: List[Optional[np.ndarray]] = []
        for smiles in smiles_list:
            if not smiles or not isinstance(smiles, str):
                results.append(None)
                continue
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    results.append(None)
                    continue
                fp  = gen.GetFingerprint(mol)
                arr = np.zeros(n_bits, dtype=np.float32)
                DataStructs.ConvertToNumpyArray(fp, arr)
                results.append(arr)
            except Exception:
                results.append(None)
        return results

    @staticmethod
    def _ecfp4_parallel(
        smiles_list: List[str],
        radius: int,
        n_bits: int,
        n_jobs: int = DEFAULT_N_JOBS,
    ) -> List[Optional[np.ndarray]]:
        """Compute ECFP4 fingerprints, using parallel workers for large inputs."""
        n = len(smiles_list)
        if n < DEFAULT_MP_THRESHOLD:
            return Splitter._ecfp4_batch(smiles_list, radius, n_bits)

        batch_size = max(1, (n + n_jobs - 1) // n_jobs)
        batches    = [smiles_list[i:i + batch_size] for i in range(0, n, batch_size)]
        results    = Parallel(n_jobs=n_jobs)(
            delayed(Splitter._ecfp4_batch)(batch, radius, n_bits)
            for batch in batches
        )
        return [fp for batch in results for fp in batch]

    @staticmethod
    def _tanimoto_dist_matrix(fps: List[np.ndarray]) -> np.ndarray:
        """Pairwise Tanimoto distance matrix for a list of ECFP4 NumPy arrays."""
        M            = np.vstack(fps)
        intersections = M @ M.T
        row_sums     = M.sum(axis=1)
        unions       = row_sums[:, None] + row_sums[None, :] - intersections
        sim          = np.where(unions > 0, intersections / unions, 0.0)
        dist         = (1.0 - sim).astype(np.float32)
        np.fill_diagonal(dist, 0.0)
        return dist

    @staticmethod
    def _isolation_scores(
        fps: List[np.ndarray],
        unit_sizes: np.ndarray,
    ) -> np.ndarray:
        """Weighted-mean Tanimoto isolation score for each unit."""
        dist_mat   = Splitter._tanimoto_dist_matrix(fps)
        total_size = unit_sizes.sum()
        weighted = dist_mat @ unit_sizes
        other_wt = total_size - unit_sizes
        return np.where(other_wt > 0, weighted / other_wt, 0.0)

    @staticmethod
    def _greedy_assign(
        sorted_units: List,
        sorted_sizes: List[int],
        n_total: int,
        test_ratio: float,
        val_ratio: float,
        train_label: str,
        val_label: str,
        test_label: str,
    ) -> Dict:
        """Greedily fill test → val → train from an ordered list of units.

        Keys of the returned mapping are the unit identifiers passed in
        ``sorted_units``; each value is the assigned split label.
        """
        n_target_test = max(1, round(n_total * test_ratio))
        n_target_val  = max(1, round(n_total * val_ratio))
        assignments: Dict = {}
        n_test = n_val = 0

        for unit, size in zip(sorted_units, sorted_sizes):
            if n_test < n_target_test:
                assignments[unit] = test_label
                n_test += size
            elif n_val < n_target_val:
                assignments[unit] = val_label
                n_val += size
            else:
                assignments[unit] = train_label
        return assignments

    @staticmethod
    def _random_order_assign(
        units: List,
        sizes: List[int],
        n_total: int,
        test_ratio: float,
        val_ratio: float,
        seed: int,
        train_label: str,
        val_label: str,
        test_label: str,
    ) -> Dict:
        """Assign units to splits in a seeded random order."""
        rng   = np.random.default_rng(seed)
        order = rng.permutation(len(units)).tolist()
        return Splitter._greedy_assign(
            [units[i] for i in order],
            [sizes[i] for i in order],
            n_total,
            test_ratio,
            val_ratio,
            train_label,
            val_label,
            test_label,
        )

    @staticmethod
    def _ks_test(
        a: np.ndarray,
        b: np.ndarray,
    ) -> Tuple[Optional[float], Optional[float]]:
        """Two-sample Kolmogorov-Smirnov test; returns ``(stat, pval)`` or ``(None, None)``."""
        if len(a) < 2 or len(b) < 2:
            return None, None
        try:
            from scipy.stats import ks_2samp
            stat, pval = ks_2samp(a, b)
            return round(float(stat), 4), round(float(pval), 4)
        except Exception:
            return None, None

    # ------------------------------------------------------------------
    # Scaffold diversity
    # ------------------------------------------------------------------

    @staticmethod
    def _scaffold_entropy(scaffold_smiles: pd.Series) -> Optional[float]:
        """Scaled Shannon entropy of the scaffold frequency distribution.

        Based on Medina-Franco et al. (2009): *Scaffold Diversity Analysis of
        Compound Data Sets Using an Entropy-Based Measure.*  The entropy is
        normalised by ``log2(n_unique)`` so that values are in [0, 1]:
        1.0 = perfectly uniform (all scaffolds appear equally often);
        0.0 = all molecules share a single scaffold.

        Acyclic molecules (empty scaffold SMILES) are excluded.  Returns
        ``None`` when fewer than 2 unique scaffolds are present.

        Parameters
        ----------
        scaffold_smiles:
            Series of Bemis-Murcko scaffold SMILES, empty string for acyclic.

        Returns
        -------
        float in [0, 1] or ``None``.
        """
        non_empty = scaffold_smiles[scaffold_smiles != ""]
        counts    = non_empty.value_counts().values   # sorted descending
        n         = len(counts)
        if n < 2:
            return None
        p = counts / counts.sum()
        h = -float(np.sum(p * np.log2(p)))
        return round(h / np.log2(n), 4)

    # ------------------------------------------------------------------
    # Pipeline entry point
    # ------------------------------------------------------------------

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        split_col = self.split_column
        train_label = self.train_label
        val_label = self.val_label
        test_label = self.test_label

        n = len(data)

        # Advisory: small datasets can yield unreliable or underpopulated splits.
        if n < self.low_n_threshold:
            self.log(
                f"Dataset has {n:,} molecules, below low_n_threshold="
                f"{self.low_n_threshold:,}. The split may be unreliable and can "
                f"leave a set underpopulated, especially method='isolation', "
                f"whose ranking is sensitive to sparse chemical coverage. Treat "
                f"the resulting partition with caution.",
                level="WARNING",
            )

        # Resolve the activity column from the upstream curator when it was not
        # pinned explicitly.  The curator advertises its standardized activity
        # column (e.g. pIC50/pEC50) via the ``output_label`` metadata key.
        if self.activity_column is None:
            res = self.get_actor_result('chembl')
            if res is not None:
                self.activity_column = res.metadata.get('output_label')

        self.log(
            f"Splitter (unit={self.unit}, method={self.method}): {n:,} molecules → "
            f"{self.test_ratio:.0%} test / {self.val_ratio:.0%} val / "
            f"{1 - self.test_ratio - self.val_ratio:.0%} train."
        )

        # Scaffold SMILES power the report card for both units.
        scaffold_series = data[self.scaffold_column].fillna("").astype(str)

        # 1. Group rows into units and identify rows forced to train.
        unit_reps, unit_indices, forced_train_idx = self._group_units(data)
        unit_sizes = [len(ix) for ix in unit_indices]
        self.log(
            f"  {len(unit_reps):,} {self.unit} units "
            f"({len(forced_train_idx):,} rows → train unconditionally)."
        )

        # 2. Order and assign units per method (no fallback).
        assignments, n_unit_fp_failed = self._assign_units(unit_reps, unit_sizes, n)

        # 3. Map unit assignments back to rows; forced-train rows stay train.
        split_labels = [train_label] * n
        for pos, label in assignments.items():
            for idx in unit_indices[pos]:
                split_labels[idx] = label
        df = data.copy()
        df[split_col] = split_labels

        counts = df[split_col].value_counts()
        self.log(
            f"  Split: train={counts.get(train_label, 0):,}  "
            f"val={counts.get(val_label, 0):,}  "
            f"test={counts.get(test_label, 0):,}"
        )

        # 4. Write report card
        report      = self._build_report_card(
            df,
            scaffold_series,
            f"{self.unit}/{self.method}",
            n_scaffold_fp_failed=n_unit_fp_failed,
        )
        report_path = self._get_run_path(self.report_card_filename)
        with open(report_path, "w") as fh:
            json.dump(report, fh, indent=2)
        self.log(f"  Report card → {report_path}")
        self._log_summary(report)

        # 6. Save figure alongside report card
        fig_path = self._save_report_figure(df, scaffold_series, report)
        if fig_path:
            self.log(f"  Report figure → {fig_path}")

        return df

    # ------------------------------------------------------------------
    # ActorOutput
    # ------------------------------------------------------------------

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        split_col = self.split_column
        counts = data[split_col].value_counts() if split_col in data.columns else {}
        return ActorOutput(
            data=data,
            success=True,
            metadata={
                f"n_{self.train_label}": int(counts.get(self.train_label, 0)),
                f"n_{self.val_label}":   int(counts.get(self.val_label,   0)),
                f"n_{self.test_label}":  int(counts.get(self.test_label,  0)),
                "split_column": split_col,
            },
            endpoint=self._get_run_path(self.report_card_filename),
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _log_summary(self, report: dict) -> None:
        """Log a formatted summary table after a split run."""
        splits = [self.train_label, self.val_label, self.test_label]

        # ── Split composition ─────────────────────────────────────────────
        rows = {
            "n_molecules": {s: report.get(f"n_{s}") for s in splits},
            "fraction":    {s: report.get(f"{s}_ratio") for s in splits},
            "n_scaffolds": {
                s: (report.get(f"scaffolds_{s}") or {}).get("n_unique")
                for s in splits
            },
            "entropy":     {
                s: (report.get(f"scaffolds_{s}") or {}).get("scaffold_entropy")
                for s in splits
            },
        }
        comp_df = pd.DataFrame(rows).T
        comp_df.columns = splits

        # ── NN distance matrix (median) ───────────────────────────────────
        nn_d    = report.get("nn_distances", {})
        nn_rows = {}
        for query in splits:
            nn_rows[query] = {}
            for ref in splits:
                v = (nn_d.get(f"{query}_{ref}") or {}).get("nn_dist_median")
                nn_rows[query][ref] = f"{v:.3f}" if v is not None else "n/a"
        nn_df = pd.DataFrame(nn_rows).T
        nn_df.columns = splits
        nn_df.index.name = "query \\ ref"

        try:
            self.log(
                f"\n  ── Split Summary ──\n"
                f"{comp_df.to_markdown(floatfmt='.3f')}\n\n"
                f"  ── NN Tanimoto Distance (mol-level, median) ──\n"
                f"{nn_df.to_markdown()}"
            )
        except Exception:
            # tabulate not installed — fall back to plain repr
            self.log(f"\n{comp_df}\n\n{nn_df}")

    def _group_units(
        self,
        data: pd.DataFrame,
    ) -> Tuple[List[str], List[List[int]], List[int]]:
        """Group rows into assignment units per ``unit``.

        Returns ``(unit_reps, unit_indices, forced_train_idx)`` where
        ``unit_reps[i]`` is the representative SMILES fingerprinted for unit
        ``i`` under ``method='isolation'``, ``unit_indices[i]`` are the row
        positions belonging to that unit, and ``forced_train_idx`` are the row
        positions assigned to train unconditionally.

        For ``unit='scaffold'`` the units are the distinct non-empty scaffold
        SMILES (the representative is the scaffold itself); acyclic and
        failed-scaffold rows are forced to train.  For ``unit='molecule'`` each
        row is its own unit (the representative is its ``SMILES_column`` value)
        and no rows are forced to train.
        """
        if self.unit == "scaffold":
            scaffold_series = data[self.scaffold_column].fillna("").astype(str)
            groups: Dict[str, List[int]] = {}
            forced_train_idx: List[int] = []
            for idx, scaf in enumerate(scaffold_series):
                if scaf == "":
                    forced_train_idx.append(idx)
                else:
                    groups.setdefault(scaf, []).append(idx)
            unit_reps    = list(groups.keys())
            unit_indices = list(groups.values())
            return unit_reps, unit_indices, forced_train_idx

        # unit == "molecule": every row is its own unit.
        smiles_series = data[self.SMILES_column].fillna("").astype(str)
        unit_reps    = smiles_series.tolist()
        unit_indices = [[i] for i in range(len(data))]
        return unit_reps, unit_indices, []

    def _assign_units(
        self,
        unit_reps: List[str],
        unit_sizes: List[int],
        n_total: int,
    ) -> Tuple[Dict[int, str], int]:
        """Order units per ``method`` and greedily fill test → val → train.

        Returns ``(assignments, n_fp_failed)`` where ``assignments`` maps a unit
        position (index into ``unit_reps``) to its split label and
        ``n_fp_failed`` is the number of units whose representative fingerprint
        could not be computed (assigned to train; always 0 for ``random``).
        """
        n_units   = len(unit_reps)
        positions = list(range(n_units))

        if self.method == "random":
            assignments = self._random_order_assign(
                positions, unit_sizes, n_total,
                self.test_ratio, self.val_ratio, self.seed,
                self.train_label, self.val_label, self.test_label,
            )
            return assignments, 0

        # ── Isolation: guard the O(n²) distance matrix behind a compute bound ──
        if n_units > self.max_units_for_isolation:
            self.log(
                f"{n_units:,} {self.unit} units exceed the isolation compute "
                f"bound max_units_for_isolation={self.max_units_for_isolation:,}. "
                f"Ranking units by structural isolation builds a pairwise "
                f"Tanimoto matrix that is O(n²) in time and memory; this limit is "
                f"a deliberate guard against unbounded computation. Raise "
                f"max_units_for_isolation to permit the larger matrix, or set "
                f"method='random' for an order-independent split.",
                level="ERROR",
            )
            raise ValueError(
                f"Splitter aborted: {n_units:,} {self.unit} units exceed "
                f"max_units_for_isolation={self.max_units_for_isolation:,}. "
                f"Raise max_units_for_isolation or use method='random'."
            )

        t0  = time.time()
        fps = self._ecfp4_parallel(unit_reps, self.ecfp4_radius, self.ecfp4_n_bits)

        assignments: Dict[int, str] = {}
        valid_pos  = [i for i, fp in enumerate(fps) if fp is not None]
        failed_pos = [i for i, fp in enumerate(fps) if fp is None]

        for i in failed_pos:
            assignments[i] = self.train_label   # unresolved representative → train

        if valid_pos:
            valid_fps = [fps[i] for i in valid_pos]
            sizes     = np.array([unit_sizes[i] for i in valid_pos], dtype=float)
            scores    = self._isolation_scores(valid_fps, sizes)
            order     = np.argsort(scores)[::-1]        # most isolated first
            assignments.update(
                self._greedy_assign(
                    [valid_pos[k]      for k in order],
                    [int(sizes[k])     for k in order],
                    n_total, self.test_ratio, self.val_ratio,
                    self.train_label, self.val_label, self.test_label,
                )
            )

        self.log(
            f"  Isolation scores computed for {n_units:,} units "
            f"in {time.time() - t0:.2f}s."
        )
        return assignments, len(failed_pos)

    def _build_report_card(
        self,
        df: pd.DataFrame,
        scaffold_series: pd.Series,
        assignment_method: str,
        n_scaffold_fp_failed: int = 0,
    ) -> dict:
        """Assemble and return the JSON report-card dictionary."""
        split_col = self.split_column
        train_label = self.train_label
        val_label = self.val_label
        test_label = self.test_label

        n_total    = len(df)
        train_mask = df[split_col] == train_label
        val_mask   = df[split_col] == val_label
        test_mask  = df[split_col] == test_label
        n_train    = int(train_mask.sum())
        n_val      = int(val_mask.sum())
        n_test     = int(test_mask.sum())

        def _pair_key(left: str, right: str) -> str:
            return f"{left}_{right}"

        warnings: List[str] = []

        # ── Scaffold stats ────────────────────────────────────────────────────
        scaf_col = scaffold_series.fillna("").astype(str)
        n_acyclic_or_missing_scaffold_rows = int((scaf_col == "").sum())

        def _scaffold_stats(mask: pd.Series) -> dict:
            valid = scaf_col[mask]
            valid = valid[valid != ""]
            cnts  = valid.value_counts()
            return {
                "n_unique":            int(valid.nunique()),
                "singleton_count":     int((cnts == 1).sum()),
                "singleton_fraction":  round(float((cnts == 1).mean()), 4) if len(cnts) else 0.0,
                "scaffold_entropy":    self._scaffold_entropy(valid),
            }

        train_scaffolds = set(scaf_col[train_mask]) - {""}
        test_scaffolds  = set(scaf_col[test_mask])  - {""}
        val_scaffolds   = set(scaf_col[val_mask])   - {""}
        n_overlap_test_train = len(test_scaffolds & train_scaffolds)
        n_overlap_val_train  = len(val_scaffolds  & train_scaffolds)
        n_overlap_test_val   = len(test_scaffolds & val_scaffolds)

        # Under unit='scaffold' whole clusters move together, so any scaffold
        # shared between test and train signals a pipeline fault.  Under
        # unit='molecule' scaffolds may span splits by design, so the overlap is
        # recorded but carries no warning.
        if self.unit == "scaffold" and n_overlap_test_train > 0:
            warnings.append(
                f"WARN: {n_overlap_test_train} {test_label} scaffold(s) also appear "
                f"in {train_label} "
                "(unexpected — check pipeline)."
            )

        # ── ECFP4 molecule-level NN distances — all 6 directed pairs ────────────
        # FPs for all molecules are computed once; the full pairwise Tanimoto
        # distance matrix is built once; all 6 directed pair stats are extracted
        # by slicing.  Skipped when the dataset exceeds max_nn_dist_n rows.
        _null_nn     = dict(nn_dist_mean=None, nn_dist_median=None,
                            nn_dist_p10=None,  nn_dist_p90=None)
        nn_distances = {k: _null_nn for k in (
            _pair_key(test_label, train_label),
            _pair_key(val_label, train_label),
            _pair_key(test_label, val_label),
            _pair_key(train_label, test_label),
            _pair_key(train_label, val_label),
            _pair_key(val_label, test_label),
            _pair_key(train_label, train_label),
            _pair_key(val_label, val_label),
            _pair_key(test_label, test_label),
        )}

        molecule_fp_qc: Dict[str, object] = {
            "computed": False,
            "reason": f"dataset size {n_total:,} > max_nn_dist_n={self.max_nn_dist_n:,}",
            "n_valid_total": None,
            "n_failed_total": None,
            "by_split": {},
        }

        if n_total <= self.max_nn_dist_n:
            smiles_col = df[self.SMILES_column]
            all_fps    = self._ecfp4_parallel(
                smiles_col.tolist(), self.ecfp4_radius, self.ecfp4_n_bits,
            )

            # Track which rows have a valid (non-None) fingerprint
            valid_arr = np.array([fp is not None for fp in all_fps])
            valid_idx = np.where(valid_arr)[0]          # global positions → dist_mat rows
            valid_fps = [all_fps[i] for i in valid_idx]

            n_valid_total = int(valid_arr.sum())
            n_failed_total = int(n_total - n_valid_total)
            molecule_fp_qc = {
                "computed": True,
                "reason": None,
                "n_valid_total": n_valid_total,
                "n_failed_total": n_failed_total,
                "by_split": {
                    train_label: {
                        "n_total": n_train,
                        "n_valid": int((train_mask.values & valid_arr).sum()),
                        "n_failed": int(n_train - (train_mask.values & valid_arr).sum()),
                    },
                    val_label: {
                        "n_total": n_val,
                        "n_valid": int((val_mask.values & valid_arr).sum()),
                        "n_failed": int(n_val - (val_mask.values & valid_arr).sum()),
                    },
                    test_label: {
                        "n_total": n_test,
                        "n_valid": int((test_mask.values & valid_arr).sum()),
                        "n_failed": int(n_test - (test_mask.values & valid_arr).sum()),
                    },
                },
            }
            if n_failed_total > 0:
                warnings.append(
                    f"WARN: {n_failed_total} molecule(s) failed ECFP4 generation and were excluded from NN distance stats."
                )

            if len(valid_fps) > 1:
                dist_mat = self._tanimoto_dist_matrix(valid_fps)  # (n_valid × n_valid)

                def _pos(mask: pd.Series) -> np.ndarray:
                    """Global split indices → rows/cols in dist_mat."""
                    global_idx = np.where(mask.values & valid_arr)[0]
                    return np.searchsorted(valid_idx, global_idx)

                def _nn_stats(
                    q: np.ndarray, r: np.ndarray, exclude_self: bool = False
                ) -> dict:
                    """Nearest-neighbour Tanimoto-distance summary from ``q`` to ``r``.

                    With ``exclude_self=True`` the query and reference sets are the
                    same split, so the diagonal (each molecule matching itself) is
                    masked out to measure within-split diversity.
                    """
                    if exclude_self:
                        if len(q) < 2:
                            return _null_nn
                        sub = dist_mat[np.ix_(q, r)].copy()
                        np.fill_diagonal(sub, np.inf)
                        nn = sub.min(axis=1)
                    else:
                        if not len(q) or not len(r):
                            return _null_nn
                        nn = dist_mat[np.ix_(q, r)].min(axis=1)
                    return dict(
                        nn_dist_mean   = round(float(nn.mean()),             4),
                        nn_dist_median = round(float(np.median(nn)),         4),
                        nn_dist_p10    = round(float(np.percentile(nn, 10)), 4),
                        nn_dist_p90    = round(float(np.percentile(nn, 90)), 4),
                    )

                tr = _pos(train_mask)
                va = _pos(val_mask)
                te = _pos(test_mask)
                nn_distances = {
                    # Inter-split (directed)
                    _pair_key(test_label, train_label):  _nn_stats(te, tr),
                    _pair_key(val_label, train_label):   _nn_stats(va, tr),
                    _pair_key(test_label, val_label):    _nn_stats(te, va),
                    _pair_key(train_label, test_label):  _nn_stats(tr, te),
                    _pair_key(train_label, val_label):   _nn_stats(tr, va),
                    _pair_key(val_label, test_label):    _nn_stats(va, te),
                    # Intra-split (within-split diversity, self excluded)
                    _pair_key(train_label, train_label): _nn_stats(tr, tr, exclude_self=True),
                    _pair_key(val_label, val_label):     _nn_stats(va, va, exclude_self=True),
                    _pair_key(test_label, test_label):   _nn_stats(te, te, exclude_self=True),
                }
        else:
            warnings.append(
                f"INFO: Molecule-level NN distance stats skipped "
                f"(dataset size {n_total:,} > max_nn_dist_n={self.max_nn_dist_n:,})."
            )

        # ── Activity distribution stats (optional) ────────────────────────────
        activity_stats: Dict = {}
        if self.activity_column and self.activity_column in df.columns:
            acts    = pd.to_numeric(df[self.activity_column], errors="coerce")
            tr_acts = acts[train_mask].dropna().values
            va_acts = acts[val_mask].dropna().values
            te_acts = acts[test_mask].dropna().values

            def _summary(arr: np.ndarray) -> dict:
                if len(arr) == 0:
                    return dict(mean=None, std=None, n=0)
                return dict(
                    mean=round(float(arr.mean()), 4),
                    std =round(float(arr.std()),  4),
                    n   =int(len(arr)),
                )

            ks_stat, ks_pval = self._ks_test(tr_acts, te_acts)
            activity_stats = {
                "activity_column":    self.activity_column,
                train_label:          _summary(tr_acts),
                val_label:            _summary(va_acts),
                test_label:           _summary(te_acts),
                f"{train_label}_{test_label}_ks_stat": ks_stat,
                f"{train_label}_{test_label}_ks_pval": ks_pval,
            }

        return {
            # Method
            "assignment_method":         assignment_method,
            "max_units_for_isolation":   self.max_units_for_isolation,
            # Split composition
            "n_total":     n_total,
            f"n_{train_label}": n_train,
            f"n_{val_label}":   n_val,
            f"n_{test_label}":  n_test,
            f"{train_label}_ratio": round(n_train / n_total, 4) if n_total else 0.0,
            f"{val_label}_ratio":   round(n_val   / n_total, 4) if n_total else 0.0,
            f"{test_label}_ratio":  round(n_test  / n_total, 4) if n_total else 0.0,
            # Scaffold stats
            "n_unique_scaffolds_total":    int(scaf_col[scaf_col != ""].nunique()),
            "n_acyclic_or_missing_scaffold_rows": n_acyclic_or_missing_scaffold_rows,
            "n_scaffold_fp_failed_clusters": int(n_scaffold_fp_failed),
            "scaffold_entropy_total":      self._scaffold_entropy(scaf_col),
            f"scaffolds_{train_label}": _scaffold_stats(train_mask),
            f"scaffolds_{val_label}":   _scaffold_stats(val_mask),
            f"scaffolds_{test_label}":  _scaffold_stats(test_mask),
            f"scaffold_overlap_{test_label}_{train_label}": n_overlap_test_train,
            f"scaffold_overlap_{val_label}_{train_label}":  n_overlap_val_train,
            f"scaffold_overlap_{test_label}_{val_label}":   n_overlap_test_val,
            # Fingerprint QC / exclusions
            "molecule_fp_qc": molecule_fp_qc,
            # Structural separation — molecule-level ECFP4 NN distances
            # All 6 directed pairs; None values when dataset > max_nn_dist_n.
            "nn_distances":  nn_distances,
            # Activity stats (empty dict when no activity column)
            **activity_stats,
            # Diagnostics
            "warnings": warnings,
        }

    def _save_report_figure(
        self,
        df: pd.DataFrame,
        scaffold_series: pd.Series,
        report: dict,
    ) -> str:
        """Generate and save a split overview figure alongside the report card.

        Layout (GridSpec 2 × 4):

        +──────────────────+──────────────────+
        | (a) Activity     | (b) Mol count    |   ← 8-wide × 4-tall each
        +──────────────────+────────+─────────+
        | (c) Scaffold     | (d1)   | (d2)    |   ← (c) 8-wide; (d1/d2) 4×4 square
        |     leakage      | Retrie-| NN dist |
        |                  | val    | heatmap |
        +──────────────────+────────+---------+

        Panels:
        (a) activity_column stepfilled density histograms per split.
        (b) Molecule count bar chart with count labels.
        (c) Scaffold leakage vs train — stacked bar (unique / shared-with-train)
            with leakage percentage for val and test.
        (d1) Scaffold retrieval curves — cumulative fraction of molecules covered
             vs fraction of unique scaffolds examined (most frequent first), one
             line per split + total dataset.  A diagonal reference shows a
             perfectly uniform distribution.
        (d2) Pairwise ECFP4 nearest-neighbour Tanimoto distance heatmap (3 × 3).
             Rows = query split, columns = reference split.  Off-diagonal cells
             show the median distance from each molecule in the query split to its
             nearest neighbour in the reference split.  Diagonal cells show the
             within-split NN distance (self excluded), reflecting intra-split
             chemical diversity.  Rendered only when the dataset is within
             ``max_nn_dist_n``; otherwise all cells show ``n/a``.

        If matplotlib is not installed the method logs a warning and returns ``""``.

        Returns
        -------
        str
            Absolute path of the saved PNG, or ``""`` if skipped.
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.gridspec as gridspec
            from matplotlib.colors import LinearSegmentedColormap
        except ImportError:
            self.log(
                "matplotlib not installed — skipping report card figure.",
                level="WARNING",
            )
            return ""

        split_col = self.split_column
        split_names = [self.train_label, self.val_label, self.test_label]
        colors      = [self._split_colors[n] for n in split_names]
        masks       = {
            self.train_label: df[split_col] == self.train_label,
            self.val_label:   df[split_col] == self.val_label,
            self.test_label:  df[split_col] == self.test_label,
        }
        scaf_col = scaffold_series.fillna("").astype(str)

        # ── Figure / GridSpec ─────────────────────────────────────────────
        fig = plt.figure(figsize=self._FIG_SIZE)
        gs  = gridspec.GridSpec(2, 4, figure=fig, wspace=0.42, hspace=0.45)
        ax_a  = fig.add_subplot(gs[0, :2])   # activity        — 8 wide
        ax_b  = fig.add_subplot(gs[0, 2:])   # mol count       — 8 wide
        ax_c  = fig.add_subplot(gs[1, :2])   # leakage         — 8 wide
        ax_d1 = fig.add_subplot(gs[1,  2])   # retrieval curve — 4 wide (square)
        ax_d2 = fig.add_subplot(gs[1,  3])   # NN-distance heatmap — 4 wide (square)
        fig.suptitle("Scaffold Split Overview", fontsize=self._FONT_TITLE, fontweight="bold")

        # ── Panel (a): activity distribution ─────────────────────────────
        has_activity = bool(
            self.activity_column and self.activity_column in df.columns
        )
        if has_activity:
            acts = pd.to_numeric(df[self.activity_column], errors="coerce")
            for color, name in zip(colors, split_names):
                vals = acts[masks[name]].dropna().values
                if vals.size:
                    ax_a.hist(
                        vals, bins=30, density=True,
                        histtype="stepfilled", alpha=self._HIST_ALPHA,
                        linewidth=1.5, color=color,
                        label=f"{name} (n={vals.size:,})",
                    )
            ax_a.set_xlabel(self.activity_column)
            ax_a.set_ylabel("Density")
            ax_a.legend(fontsize=self._FONT_SMALL)
        else:
            ax_a.text(0.5, 0.5, "No activity column",
                      ha="center", va="center", transform=ax_a.transAxes,
                      fontsize=self._FONT_MEDIUM, color="grey")
        ax_a.set_title("Activity Distribution by Split")

        # ── Panel (b): molecule count per split ───────────────────────────
        n_mols   = [int(masks[n].sum()) for n in split_names]
        fp_qc    = report.get("molecule_fp_qc") or {}
        by_split = fp_qc.get("by_split") if isinstance(fp_qc, dict) else None

        if fp_qc.get("computed") and isinstance(by_split, dict):
            n_valid = [int((by_split.get(n) or {}).get("n_valid", 0)) for n in split_names]
            n_failed = [int((by_split.get(n) or {}).get("n_failed", 0)) for n in split_names]
            max_mols = max((v + f for v, f in zip(n_valid, n_failed)), default=1)
            bars_valid = ax_b.bar(
                split_names,
                n_valid,
                color=colors,
                edgecolor="white",
                label="valid for NN",
            )
            ax_b.bar(
                split_names,
                n_failed,
                bottom=n_valid,
                color="#c0c0c0",
                edgecolor="white",
                label="failed/removed for NN",
            )
            for bar, valid_count, failed_count in zip(bars_valid, n_valid, n_failed):
                total = valid_count + failed_count
                ax_b.text(
                    bar.get_x() + bar.get_width() / 2,
                    total + 0.01 * max_mols,
                    f"{total:,}",
                    ha="center",
                    va="bottom",
                    fontsize=self._FONT_SMALL,
                )
                if failed_count > 0:
                    ax_b.text(
                        bar.get_x() + bar.get_width() / 2,
                        valid_count + failed_count / 2,
                        f"-{failed_count:,}",
                        ha="center",
                        va="center",
                        fontsize=self._FONT_SMALL,
                        color="dimgrey",
                    )
            ax_b.legend(fontsize=self._FONT_SMALL, frameon=False)
            ax_b.set_title("Molecule Count per Split (NN QC)")
            ax_b.set_ylabel("Count")
            ax_b.set_ylim(0, max_mols * 1.18)
        else:
            max_mols = max(n_mols, default=1)
            bars = ax_b.bar(split_names, n_mols, color=colors, edgecolor="white")
            for bar, count in zip(bars, n_mols):
                ax_b.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01 * max_mols,
                    f"{count:,}", ha="center", va="bottom", fontsize=self._FONT_SMALL,
                )
            ax_b.set_title("Molecule Count per Split")
            ax_b.set_ylabel("Count")
            ax_b.set_ylim(0, max_mols * 1.15)

        # ── Panel (c): scaffold leakage vs train ──────────────────────────
        train_scaffolds = set(scaf_col[masks[self.train_label]].values) - {""}
        n_unique: List[int] = []
        n_shared: List[int] = []
        for name in split_names:
            total = int(masks[name].sum())
            if name == self.train_label or not train_scaffolds:
                n_unique.append(total)
                n_shared.append(0)
            else:
                non_empty    = scaf_col[masks[name]][scaf_col[masks[name]] != ""]
                shared_valid = int(non_empty.isin(train_scaffolds).sum())
                n_unique.append(total - shared_valid)
                n_shared.append(shared_valid)

        ax_c.bar(split_names, n_unique, color=colors, edgecolor="white",
                 label="unique scaffold")
        ax_c.bar(split_names, n_shared, bottom=n_unique,
                 color=self._OVERLAP_COLOR, edgecolor="white", label="scaffold in train")

        for i, (name, nu, ns) in enumerate(zip(split_names, n_unique, n_shared)):
            total = nu + ns
            if total > 0 and name != self.train_label:
                pct = 100.0 * ns / total
                ax_c.text(
                    i, total + 0.01 * max(n_mols, default=1),
                    f"{pct:.0f}% leak", ha="center", va="bottom",
                    fontsize=self._FONT_SMALL, color="dimgrey",
                )

        max_stack = max((u + s for u, s in zip(n_unique, n_shared)), default=1)
        ax_c.set_title("Scaffold Leakage vs. Train")
        ax_c.set_ylabel("Molecule Count")
        ax_c.set_ylim(0, max_stack * 1.18)
        ax_c.legend(fontsize=self._FONT_SMALL, frameon=False)

        # ── Panel (d1): scaffold retrieval curves ─────────────────────────
        # x = fraction of unique scaffolds examined (sorted by descending frequency)
        # y = cumulative fraction of molecules covered by those scaffolds
        # A high-diversity set stays near the diagonal; a scaffold-dominated set
        # rises steeply at the left.
        def _retrieval_curve(scaf: pd.Series):
            ne = scaf[scaf != ""]
            if ne.empty:
                return np.array([0., 1.]), np.array([0., 1.])
            counts = ne.value_counts().values         # descending
            n_u    = len(counts)
            x      = np.arange(1, n_u + 1) / n_u
            y      = np.cumsum(counts) / counts.sum()
            return x, y

        # Diagonal reference — x = y means every molecule has a unique scaffold
        ax_d1.plot([0, 1], [0, 1], color="grey", linestyle="--",
                   linewidth=1.0, label="uniform")
        ax_d1.text(
            0.5, 0.5, "perfectly uniform",
            color="grey", fontsize=self._FONT_SMALL,
            ha="center", va="bottom",
            rotation=45, rotation_mode="anchor",
            transform=ax_d1.transData,
        )

        # Total dataset curve first, then per-split
        x_tot, y_tot = _retrieval_curve(scaf_col)
        ax_d1.plot(x_tot, y_tot, color="black", linewidth=1.5,
                   label=f"total (n={len(scaf_col[scaf_col != '']):,})")

        for color, name in zip(colors, split_names):
            x_s, y_s = _retrieval_curve(scaf_col[masks[name]])
            n_s      = int((scaf_col[masks[name]] != "").sum())
            ax_d1.plot(x_s, y_s, color=color, linewidth=1.2,
                       label=f"{name} (n={n_s:,})")

        ax_d1.set_xlim(0, 1)
        ax_d1.set_ylim(0, 1)
        ax_d1.set_xlabel("Fraction of unique scaffolds")
        ax_d1.set_ylabel("Cumulative fraction of molecules")
        ax_d1.set_title("Scaffold Retrieval")
        ax_d1.legend(fontsize=self._FONT_SMALL, frameon=False)

        # ── Panel (d2): pairwise NN distance heatmap (3 × 3) ─────────────
        # Rows = query split, cols = reference split.
        # Off-diagonal: median ECFP4 NN Tanimoto distance from each molecule
        #   in the query split to its nearest neighbour in the reference split.
        # Diagonal: within-split NN distance (self excluded) — a measure of
        #   intra-split chemical diversity.

        tr, va, te = self.train_label, self.val_label, self.test_label
        _key_map = {
            (0, 0): f"{tr}_{tr}", (0, 1): f"{tr}_{va}", (0, 2): f"{tr}_{te}",
            (1, 0): f"{va}_{tr}", (1, 1): f"{va}_{va}", (1, 2): f"{va}_{te}",
            (2, 0): f"{te}_{tr}", (2, 1): f"{te}_{va}", (2, 2): f"{te}_{te}",
        }
        nn_d   = report.get("nn_distances", {})
        matrix = np.full((3, 3), np.nan)
        for (i, j), key in _key_map.items():
            v = (nn_d.get(key) or {}).get("nn_dist_median")
            if v is not None:
                matrix[i, j] = v
        
        hm_cmap = LinearSegmentedColormap.from_list("wh_teal", self._HM_CMAP_COLORS)
        im = ax_d2.imshow(matrix, cmap=hm_cmap, vmin=0.0, vmax=1.0, aspect="auto")

        ax_d2.set_xticks(range(3))
        ax_d2.set_yticks(range(3))
        ax_d2.set_xticklabels(split_names, fontsize=self._FONT_SMALL)
        ax_d2.set_yticklabels(split_names, fontsize=self._FONT_SMALL)
        ax_d2.set_xlabel("Reference", fontsize=self._FONT_SMALL)
        ax_d2.set_ylabel("Query",     fontsize=self._FONT_SMALL)
        ax_d2.set_title("Median Nearest Neighbor Distance")

        for i in range(3):
            for j in range(3):
                v = matrix[i, j]
                txt   = f"{v:.3f}" if not np.isnan(v) else "n/a"
                color = "white" if (not np.isnan(v) and v > 0.6) else "black"
                ax_d2.text(j, i, txt, ha="center", va="center",
                           fontsize=self._FONT_SMALL, color=color)

        plt.colorbar(im, ax=ax_d2, fraction=0.046, pad=0.04,
                     label="Tanimoto distance")
        
        ax_d2.set_aspect(1)
        ax_d2.set_box_aspect(1)
        # ── Shared formatting ─────────────────────────────────────────────
        for ax in [ax_a, ax_b, ax_c, ax_d1, ax_d2]:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=self._FONT_SMALL)

        fig_filename = self.report_card_filename.replace(".json", ".png")
        fig_path     = self._get_run_path(fig_filename)
        plt.savefig(fig_path, dpi=self._FIG_DPI, bbox_inches="tight")
        plt.close()
        return fig_path
