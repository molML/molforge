"""
ChEMBL SQL backend implementation.

Wraps the ChEMBLSQL actor logic for use with unified ChEMBLSource actor.
"""

import sqlite3
import pandas as pd
from tqdm.auto import tqdm
from typing import Optional

from .base import ChEMBLSourceBackend


class ChEMBLSQLBackend(ChEMBLSourceBackend):
    """
    SQL backend for ChEMBL data retrieval.

    Queries a local ChEMBL SQLite database for activity data.
    """
    __backend_name__ = 'SQL'

    # Backend metadata
    description = "Query local ChEMBL SQLite database"
    required_params = []
    optional_params = ['db_path', 'version', 'auto_download', 'verbose', 'n_jobs']

    @property
    def source_description(self) -> str:
        """Data source description for metadata."""
        return 'ChEMBL SQL'

    def fetch_activities(self, target_chembl_id: str) -> pd.DataFrame:
        """
        Fetch activity data from local ChEMBL database.

        Args:
            target_chembl_id: ChEMBL identifier for the target protein

        Returns:
            DataFrame containing activity data

        Raises:
            sqlite3.Error: If database query fails
        """
        self.log(
            f"\n{'='*50}\n"
            f"|   Fetching activity data for {target_chembl_id}\t|\n"
            f"|   Database: {self.params._resolved_db_path}\n"
            f"{'='*50}"
        )

        # First, get the total count
        total_count = self._get_total_count(target_chembl_id)
        self.log(f"Total entries available: {total_count}")

        if total_count == 0:
            self.log("No activities found for this target.")
            return pd.DataFrame()

        # Fetch activities based on configuration
        if self.params.search_all:
            activities_df = self._fetch_all_entries(target_chembl_id, total_count)
        else:
            activities_df = self._fetch_n_entries(target_chembl_id, total_count)

        self.log(f"{len(activities_df)} entries retrieved.")
        return activities_df

    def _get_connection(self) -> sqlite3.Connection:
        """
        Create a database connection.

        Returns:
            SQLite database connection
        """
        conn = sqlite3.connect(self.params._resolved_db_path)
        conn.row_factory = sqlite3.Row  # Enable column access by name
        return conn

    def _get_total_count(self, target_chembl_id: str) -> int:
        """
        Get the total number of activities for a target.

        Args:
            target_chembl_id: ChEMBL identifier for the target protein

        Returns:
            Total count of activities
        """
        conn = self._get_connection()
        cursor = conn.cursor()

        query = """
        SELECT COUNT(*) as count
        FROM activities a
        JOIN assays ass ON a.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        WHERE td.chembl_id = ?
        """

        cursor.execute(query, (target_chembl_id,))
        result = cursor.fetchone()
        conn.close()

        return result['count'] if result else 0

    def _fetch_all_entries(self, target_chembl_id: str, total_count: int) -> pd.DataFrame:
        """
        Fetch all available entries for a target.

        Args:
            target_chembl_id: ChEMBL identifier for the target protein
            total_count: Total number of entries available

        Returns:
            DataFrame with all activities
        """
        conn = self._get_connection()

        query = self._build_activity_query()

        self.log(f"Fetching all {total_count} entries...")

        with tqdm(total=total_count, desc="Fetching entries") as pbar:
            # Use pandas read_sql_query which handles large result sets efficiently
            df = pd.read_sql_query(
                query,
                conn,
                params=(target_chembl_id,)
            )
            pbar.update(len(df))

        conn.close()
        return df

    def _fetch_n_entries(self, target_chembl_id: str, total_count: int) -> pd.DataFrame:
        """
        Fetch exactly `n` entries for a target.

        Args:
            target_chembl_id: ChEMBL identifier for the target protein
            total_count: Total number of entries available

        Returns:
            DataFrame with n activities
        """
        conn = self._get_connection()

        limit = min(self.params.n, total_count)
        query = self._build_activity_query(limit=limit)

        self.log(f"Fetching {limit} entries...")

        with tqdm(total=limit, desc="Fetching entries") as pbar:
            df = pd.read_sql_query(
                query,
                conn,
                params=(target_chembl_id,)
            )
            pbar.update(len(df))

        conn.close()
        return df

    def _build_activity_query(self, limit: Optional[int] = None) -> str:
        """
        Build the SQL query to fetch activity data.

        This query mimics the structure of data returned by the ChEMBL API,
        including all relevant fields from activities, assays, targets, and molecules.

        Args:
            limit: Optional limit on number of rows to return

        Returns:
            SQL query string
        """
        query = """
        SELECT
            a.activity_id,
            a.assay_id,
            a.doc_id,
            a.record_id,
            a.molregno,
            a.standard_relation,
            a.standard_value,
            a.standard_units,
            a.standard_flag,
            a.standard_type,
            a.activity_comment,
            a.data_validity_comment,
            a.potential_duplicate,
            a.pchembl_value,
            a.bao_endpoint,
            a.uo_units,
            a.qudt_units,
            a.toid,
            a.upper_value,
            a.standard_upper_value,
            a.src_id,
            a.type,
            a.relation,
            a.value,
            a.units,
            a.text_value,
            a.standard_text_value,

            -- Molecule/Compound information
            m.chembl_id as molecule_chembl_id,
            m.pref_name as molecule_pref_name,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,

            -- Assay information
            ass.assay_type,
            ass.description as assay_description,
            ass.assay_organism,
            ass.assay_tax_id,
            ass.assay_strain,
            ass.assay_tissue,
            ass.assay_cell_type,
            ass.assay_subcellular_fraction,
            ass.bao_format,

            -- Target information
            td.chembl_id as target_chembl_id,
            td.pref_name as target_pref_name,
            td.target_type,
            td.organism as target_organism,
            td.tax_id as target_tax_id,

            -- Document information
            d.chembl_id as document_chembl_id,
            d.journal,
            d.year as document_year,
            d.volume,
            d.issue,
            d.first_page,
            d.last_page,
            d.pubmed_id,
            d.doi

        FROM activities a
        JOIN assays ass ON a.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        LEFT JOIN molecule_dictionary m ON a.molregno = m.molregno
        LEFT JOIN compound_structures cs ON m.molregno = cs.molregno
        LEFT JOIN docs d ON a.doc_id = d.doc_id
        WHERE td.chembl_id = ?
        """

        if limit is not None:
            query += f"\nLIMIT {limit}"

        return query
