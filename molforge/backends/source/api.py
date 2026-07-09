"""
ChEMBL API backend implementation.

Wraps the ChEMBLAPI actor logic for use with unified ChEMBLSource actor.
"""

import time
import requests
import pandas as pd
from tqdm.auto import tqdm
from typing import Dict, Any

from .base import ChEMBLSourceBackend


class ChEMBLAPIBackend(ChEMBLSourceBackend):
    """
    API backend for ChEMBL data retrieval.

    Fetches activity data from ChEMBL web API with pagination and rate limiting.
    """
    __backend_name__ = 'API'

    # Backend metadata
    description = "Fetch data from ChEMBL REST API"
    required_params = []
    optional_params = ['n', 'verbose', 'n_jobs']

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json?"
    MAX_PAGE_SIZE = 1000
    RATE_LIMIT_DELAY = 1.1  # seconds between requests

    @property
    def source_description(self) -> str:
        """Data source description for metadata."""
        return 'ChEMBL API'

    def fetch_activities(self, target_chembl_id: str) -> pd.DataFrame:
        """
        Fetch activity data from ChEMBL web API.

        Args:
            target_chembl_id: ChEMBL identifier for the target protein

        Returns:
            DataFrame containing activity data

        Raises:
            requests.HTTPError: If API request fails
        """
        limit = min(self.params.n, self.MAX_PAGE_SIZE)
        activities_url = f"{self.BASE_URL}target_chembl_id={target_chembl_id}&limit={limit}"

        self.log(
            f"\n{'='*50}\n"
            f"|   Fetching activity data for {target_chembl_id}\t|\n"
            f"|   URL: {activities_url}\n"
            f"{'='*50}"
        )

        time.sleep(self.RATE_LIMIT_DELAY)
        response = self.get_response(activities_url)
        total_count = response['page_meta']['total_count']
        activities = response['activities']

        if self.params.search_all:
            activities = self._fetch_all_entries(activities, target_chembl_id, total_count, limit)
        else:
            activities = self._fetch_n_entries(activities, target_chembl_id, total_count, limit)

        self.log(f"{len(activities)} entries retrieved.")
        return pd.DataFrame(activities)

    def get_response(self, url: str) -> Dict[str, Any]:
        """
        Execute API request and return JSON response.

        Args:
            url: Complete API endpoint URL

        Returns:
            Parsed JSON response

        Raises:
            requests.HTTPError: If request fails
        """
        for retry in range(3):  # retry up to 3 times
            response = requests.get(url)
            if response.status_code == 500:
                self.log(f"API Server Error: Internal Server Error, retrying in 5s. ({retry +1}/3)", level="ERROR")
                time.sleep(5)
                continue

            response.raise_for_status()
            return response.json()
        response.raise_for_status()

    def _fetch_all_entries(self, activities, target_id: str, total_count: int, limit: int):
        """Fetch all available entries for a target, tracking entries not pages."""
        if total_count <= limit:
            return activities

        entries_fetched = len(activities)
        with tqdm(total=total_count, initial=entries_fetched, desc="Fetching entries") as pbar:
            pages_needed = (total_count - 1) // limit
            for i in range(1, pages_needed + 1):
                url = f"{self.BASE_URL}limit={limit}&offset={i*limit}&target_chembl_id={target_id}"
                time.sleep(self.RATE_LIMIT_DELAY)
                response = self.get_response(url)
                new_activities = response['activities']
                activities.extend(new_activities)

                pbar.update(len(new_activities))
                entries_fetched += len(new_activities)

                if len(new_activities) < limit:
                    break

        return activities

    def _fetch_n_entries(self, activities, target_id: str, total_count: int, limit: int):
        """Fetch exactly `n` entries for a target, tracking entries not pages."""
        if len(activities) >= self.params.n:
            return activities[:self.params.n]

        entries_needed = min(self.params.n, total_count)
        entries_fetched = len(activities)

        with tqdm(total=entries_needed, initial=entries_fetched, desc="Fetching entries") as pbar:
            while entries_fetched < entries_needed:
                offset = entries_fetched // limit * limit
                url = f"{self.BASE_URL}limit={limit}&offset={offset}&target_chembl_id={target_id}"
                time.sleep(self.RATE_LIMIT_DELAY)
                response = self.get_response(url)
                new_activities = response['activities']

                remaining_needed = entries_needed - entries_fetched
                to_take = new_activities[:remaining_needed]
                activities.extend(to_take)
                entries_fetched += len(to_take)
                pbar.update(len(to_take))

                if len(new_activities) < limit or entries_fetched >= entries_needed:
                    break

        return activities[:self.params.n]
