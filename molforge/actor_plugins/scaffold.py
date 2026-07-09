"""
Bemis-Murcko scaffold analysis actor plugin.

Computes Bemis-Murcko scaffolds (and optionally generic carbon-skeleton
scaffolds) for every molecule in the pipeline.

Parallel processing uses ``joblib`` with a module-level worker function so
that only plain, serialisable data crosses the process boundary — no actor
class reference is ever pickled.

Only RDKit is required.

Output columns
--------------
scaffold_smiles
    Bemis-Murcko scaffold SMILES (atom-typed: ring systems + linkers,
    side-chains stripped).
scaffold_generic_smiles
    Generic (carbon-skeleton) scaffold SMILES — all atoms replaced with C,
    all bonds made single.  Always added; contains the generic scaffold SMILES
    when ``include_generic=True``, otherwise ``None``.
scaffold_success
    ``True`` when scaffold computation succeeded, ``False`` when the input
    SMILES could not be parsed or an RDKit error occurred.

Example
-------
::

    from molforge import MolForge, ForgeParams
    from molforge.actor_plugins.scaffold import ComputeScaffoldsParams

    params = ForgeParams(
        steps=["source", "chembl", "curate", "scaffold"],
        plugin_params={"scaffold": ComputeScaffoldsParams(include_generic=True)},
        ...
    )
    forge = MolForge(params)
    result = forge.forge("CHEMBL203")
"""

import gc
import time
from dataclasses import dataclass
from typing import List, Optional, Tuple

import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric

from molforge.actors.base import BaseActor
from molforge.actors.params.base import BaseParams
from molforge.actors.protocol import ActorOutput
from molforge.utils.constants import DEFAULT_MP_THRESHOLD, DEFAULT_N_JOBS


# ---------------------------------------------------------------------------
# Module-level worker function
# ---------------------------------------------------------------------------

def _compute_scaffold_batch(
    smiles_list: List[str],
    include_generic: bool,
    include_chirality: bool,
) -> List[Tuple[Optional[str], Optional[str], bool]]:
    """
    Compute Bemis-Murcko scaffolds for a batch of SMILES strings.

    Defined at module level so ``joblib`` can serialise it without
    pickling the actor class.  All arguments are plain Python primitives.

    Parameters
    ----------
    smiles_list:
        Input SMILES strings for this batch.
    include_generic:
        Whether to also compute the generic (carbon-skeleton) scaffold.
    include_chirality:
        Whether to preserve stereo annotations in the scaffold SMILES.

    Returns
    -------
    List of ``(murcko_smiles, generic_smiles, success)`` tuples, one per
    input SMILES.  On failure the first two elements are ``None`` and
    ``success`` is ``False``.
    """
    results = []
    for smiles in smiles_list:
        if not smiles or not isinstance(smiles, str):
            results.append((None, None, False))
            continue
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                results.append((None, None, False))
                continue

            scaffold_mol = GetScaffoldForMol(mol)

            # Empty scaffold (no rings) is a valid result — stored as ""
            murcko_smiles = Chem.MolToSmiles(
                scaffold_mol,
                isomericSmiles=include_chirality,
            )

            generic_smiles = None
            if include_generic:
                generic_mol = MakeScaffoldGeneric(scaffold_mol)
                generic_smiles = Chem.MolToSmiles(generic_mol)

            results.append((murcko_smiles, generic_smiles, True))

        except Exception:
            results.append((None, None, False))

    return results


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

@dataclass
class ComputeScaffoldsParams(BaseParams):
    """
    Configuration for Bemis-Murcko scaffold computation.

    See the per-field docstrings below for details on each parameter.
    """

    SMILES_column:     str  = 'curated_smiles'
    """Input column containing SMILES strings.  Should be the output of the
    ``curate`` actor (default ``"curated_smiles"``)."""

    include_generic:   bool = True
    """Also compute the generic (carbon-skeleton) scaffold and store it in
    ``scaffold_generic_smiles``.  The column is always created; it holds
    ``None`` when this is ``False``."""

    include_chirality: bool = False
    """Preserve stereo annotations in the Murcko scaffold SMILES.  Has no
    effect on the generic scaffold, which always strips annotations."""

    dropna:            bool = False
    """Remove rows where ``scaffold_success`` is ``False`` — i.e. the input
    SMILES could not be parsed or the scaffold algorithm raised an error."""

    acyclic_policy:    str  = 'keep'
    """Controls how acyclic molecules are handled.  Acyclic molecules have no
    ring system; their Murcko scaffold is an empty string and
    ``scaffold_success`` is ``True``.

    - ``'keep'``   — retain acyclic molecules in the output (default).
    - ``'remove'`` — drop acyclic molecules from the output."""

    def _validate_params(self) -> None:
        if not self.SMILES_column:
            raise ValueError("SMILES_column cannot be empty.")
        self._validate_policy('acyclic_policy', self.acyclic_policy, ['keep', 'remove'])


# ---------------------------------------------------------------------------
# Actor
# ---------------------------------------------------------------------------

class ComputeScaffolds(BaseActor):
    """
    Compute Bemis-Murcko scaffolds for each molecule in the pipeline.

    For datasets below ``DEFAULT_MP_THRESHOLD`` rows the computation runs
    in a single process.  Larger datasets are split into batches and
    processed in parallel using ``joblib``.

    The actor adds three columns to the DataFrame:

    - ``scaffold_smiles`` — atom-typed Murcko scaffold SMILES (empty string
      for acyclic molecules, ``None`` on failure)
    - ``scaffold_generic_smiles`` — always added; contains the carbon-skeleton
      scaffold when ``include_generic=True``, otherwise ``None``
    - ``scaffold_success`` — ``True`` when the scaffold algorithm succeeded,
      including for acyclic molecules whose scaffold is an empty string

    Output filtering is controlled by two independent parameters:

    - ``dropna`` — removes rows where ``scaffold_success`` is ``False``
    - ``acyclic_policy`` — ``'remove'`` drops acyclic molecules (empty scaffold,
      ``scaffold_success=True``); ``'keep'`` retains them (default)
    """

    __step_name__    = 'scaffold'
    __dependencies__ = ['curate']
    __param_class__  = ComputeScaffoldsParams

    OUTPUT_COLUMNS = ['scaffold_smiles', 'scaffold_generic_smiles', 'scaffold_success']

    @property
    def required_columns(self) -> List[str]:
        return [self.SMILES_column]

    @property
    def output_columns(self) -> List[str]:
        return self.OUTPUT_COLUMNS

    # ------------------------------------------------------------------
    # Pipeline entry point
    # ------------------------------------------------------------------

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        smiles_list = data[self.SMILES_column].tolist()
        n = len(smiles_list)

        if n < DEFAULT_MP_THRESHOLD:
            self.log(f"Computing scaffolds for {n:,} molecules (single process).")
            results = _compute_scaffold_batch(
                smiles_list, self.include_generic, self.include_chirality
            )
        else:
            self.log(f"Computing scaffolds for {n:,} molecules (parallel).")
            results = self._compute_scaffolds_parallel(smiles_list)

        return self._combine_results(data, results)

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        n_ok = int(data['scaffold_success'].sum()) if 'scaffold_success' in data.columns else 0
        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'n_scaffolds_computed': n_ok,
                'n_unique_scaffolds': (
                    data['scaffold_smiles'].nunique()
                    if 'scaffold_smiles' in data.columns else 0
                ),
            },
        )

    # ------------------------------------------------------------------
    # Parallel processing
    # ------------------------------------------------------------------

    def _compute_scaffolds_parallel(self, smiles_list: List[str]) -> List[Tuple]:
        """Split SMILES into batches and compute scaffolds in parallel."""
        n = len(smiles_list)
        batch_size = max(1, (n + DEFAULT_N_JOBS - 1) // DEFAULT_N_JOBS)
        batches = [
            smiles_list[i:i + batch_size]
            for i in range(0, n, batch_size)
        ]

        self.log(
            f"Scaffold: {DEFAULT_N_JOBS} workers, "
            f"{len(batches)} batches of ~{batch_size:,}."
        )

        start = time.time()
        batch_results = Parallel(n_jobs=DEFAULT_N_JOBS)(
            delayed(_compute_scaffold_batch)(
                batch, self.include_generic, self.include_chirality
            )
            for batch in batches
        )
        self.log(f"Parallel scaffold complete in {time.time() - start:.1f}s.")

        gc.collect()
        return [item for batch in batch_results for item in batch]

    # ------------------------------------------------------------------
    # Result assembly
    # ------------------------------------------------------------------

    def _combine_results(
        self, df: pd.DataFrame, results: List[Tuple]
    ) -> pd.DataFrame:
        """Attach scaffold result tuples as new columns on the DataFrame."""
        murcko, generic, success = zip(*results) if results else ([], [], [])

        df = df.copy()
        df['scaffold_smiles']         = list(murcko)
        df['scaffold_generic_smiles'] = list(generic)
        df['scaffold_success']        = list(success)

        n_ok    = sum(success)
        n_empty = sum(1 for s in murcko if s == '')
        self.log(
            f"Scaffold computation: {n_ok}/{len(df)} succeeded "
            f"({n_empty} acyclic, {len(df) - n_ok} failed)."
        )

        if self.dropna:
            df = df[df['scaffold_success']]

        if self.acyclic_policy == 'remove':
            df = df[df['scaffold_smiles'].ne('')]

        return df
