"""
Utility for downloading and verifying the latest ChEMBL SQLite database.

Used by ChEMBLSQL to ensure a local database is available for querying.
"""

import os
import re
import requests
from tqdm.auto import tqdm
from typing import Optional, Union
import tarfile
        
def get_latest_chembl_version() -> int:
    """Return the latest ChEMBL release version number from the EBI FTP index."""
    url = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    r = requests.get(url)
    r.raise_for_status()
    match = re.search(r"chembl_(\d+)_sqlite", r.text)
    if match:
        return int(match.group(1))
    raise RuntimeError("Could not find latest ChEMBL version.")

def ensure_chembl_db(db_path: str, version: Union[int, str] = 'latest', logger: Optional[object] = None) -> str:
    """
    Ensure that the ChEMBL SQLite database exists locally.
    If missing, downloads and decompresses it.

    Args:
        db_path: Destination path for the SQLite database (e.g., './chembl_35.db').
        version: Which version of chembl should be downloaded, if db_path does not exist.
        logger: Optional logging-compatible object with a `.print()` method.

    Returns:
        Path to the verified ChEMBL database.
    """

    if version == 'latest':
        page = version
        version = get_latest_chembl_version()
    else:
        page = f'releases/chembl_{version}'
        
    CHEMBL_SQLITE_URL = f"https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/{page}/chembl_{version}_sqlite.tar.gz"

    if os.path.exists(db_path):
        print(f"ChEMBL database already exists at {db_path}")
        return db_path

    if not str(version) in db_path:
        print(f"WARNING: ChEMBL version {version} was requested, but version is not mentioned in db_path: {db_path}")
    
    gz_path = db_path + ".gz"
    os.makedirs(os.path.dirname(db_path) or ".", exist_ok=True)

    print(f"Downloading ChEMBL SQLite database from:\n{CHEMBL_SQLITE_URL}")

    response = requests.get(CHEMBL_SQLITE_URL, stream=True)
    response.raise_for_status()

    total_size = int(response.headers.get("content-length", 0))
    block_size = 1024 * 1024  # 1 MB chunks

    with open(gz_path, "wb") as f, tqdm(
        total=total_size, unit="B", unit_scale=True, desc="Downloading ChEMBL DB"
    ) as pbar:
        for chunk in response.iter_content(block_size):
            f.write(chunk)
            pbar.update(len(chunk))

    # Decompress the gzipped database
    print("Decompressing database...")

    with tarfile.open(gz_path, "r:gz") as tar:
        # find the .db file in the tar
        db_member = next((m for m in tar.getmembers() if m.name.endswith(".db")), None)
        if db_member is None:
            raise RuntimeError("No .db file found in tarball")
        
        # set the extraction path to your desired db_path
        db_member.name = os.path.basename(db_path)  # override to avoid nested folders
        tar.extract(db_member, path=os.path.dirname(db_path))
    
    os.remove(gz_path)

    print(f"Database saved to {db_path}")

    return db_path
