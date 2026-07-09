from typing import Optional
from dataclasses import dataclass

import os

from .base import BaseParams


@dataclass
class TokenizeDataParams(BaseParams):
    """Configuration parameters for SMILES tokenization component."""
    
    vocab_file: Optional[str] = None
    """Path to a vocabulary JSON file: loaded if the file exists, otherwise a new vocabulary is built from the data and saved to this path. When None, uses (or creates) vocab.json in the run directory."""
    SMILES_column: str = 'curated_smiles'
    """DataFrame column containing SMILES strings to tokenize."""
    dynamically_update_vocab: bool = True
    """When a non-empty vocabulary is loaded from file, extend it with any new tokens found in the data and re-save it. Ignored when a fresh vocabulary is built (which always includes all tokens from the data)."""

    def _validate_params(self) -> None:
        if not self.dynamically_update_vocab and self.vocab_file is not None and not os.path.isfile(self.vocab_file):
            raise ValueError(f"Invalid vocab_file. Does not exist: {self.vocab_file}")