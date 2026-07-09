"""SMILES tokenization actor."""

import os
import pandas as pd

from typing import List

from .base import BaseActor
from .protocol import ActorOutput
from ..configuration.steps import Steps

from ..utils.actortools.vocabulary import (
    Vocabulary,
    SMILESTokenizer,
    load_vocabulary,
    save_vocabulary
)


class TokenizeData(BaseActor):
    """SMILES tokenization and vocabulary management."""
    __step_name__ = Steps.TOKENS

    def __post_init__(self):
        """Post-initialization setup for tokenizer and vocabulary."""
        self.TOKENIZER = SMILESTokenizer()
        self.VOCABULARY = self._load_or_create_vocabulary()

    @property
    def required_columns(self) -> List[str]:
        """Required input columns."""
        return [self.SMILES_column]

    @property
    def output_columns(self) -> List[str]:
        """Output columns from tokenization."""
        return ['tokens', 'seqlen']

    @property
    def output_vocab_file(self) -> str:
        """Path to output vocabulary file (writes to input vocab_file if specified, else run directory)."""
        if self.vocab_file:
            return self.vocab_file
        # Context may not be available during initialization
        if self._context and hasattr(self._context, 'output_dir'):
            return self._get_run_path("vocab.json")
        return None

    def _load_or_create_vocabulary(self):
        """Load vocabulary from file or create empty one."""
        vocab_path = self.output_vocab_file

        if vocab_path and os.path.isfile(vocab_path):
            self.log(f'Loading vocab from {vocab_path}')
            return load_vocabulary(vocab_path)
        else:
            self.log(f'Creating new vocabulary.')
            return Vocabulary()


    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Perform SMILES-based tokenization for a dataframe.

        Args:
            data: Input DataFrame with SMILES column

        Returns:
            DataFrame with tokens and seqlen columns added

        Raises:
            ValueError: If SMILES column not found
        """
        if self.SMILES_column not in data.columns:
            raise ValueError(
                f"'{self.SMILES_column}' not in dataframe. "
                f"Available columns: {list(data.columns)}"
            )

        data = data.copy()  # do not edit original

        SMILES = data[self.SMILES_column]
        TOKENS = SMILES.apply(self.TOKENIZER.tokenize)

        # build vocab if its a fresh one
        if len(self.VOCABULARY) == 0:
            self._build_new_vocabulary(TOKENS)
            vocab_tokens = sorted(self.VOCABULARY.tokens())
            self.log(f"Vocabulary constructed: {len(self.VOCABULARY)} tokens: {vocab_tokens}")

        elif self.dynamically_update_vocab:
            old_tokens = set(self.VOCABULARY.tokens())
            self.VOCABULARY.update(self._extract_token_set(TOKENS))
            new_tokens = set(self.VOCABULARY.tokens()) - old_tokens

            if new_tokens:
                new_tokens_sorted = sorted(new_tokens)
                self.log(f"Vocabulary expanded: {len(old_tokens)} -> {len(self.VOCABULARY)} tokens. Added: {new_tokens_sorted}")
                self._write_vocabulary()

        data['tokens'] = TOKENS
        data['seqlen'] = data.tokens.apply(len)

        self.log(f"Tokenized {len(TOKENS)} sequences.")
        return data

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with vocabulary metadata."""
        return ActorOutput(
            data=data,
            success=True,
            metadata={
                'vocab_size': len(self.VOCABULARY),
                'n_tokenized': len(data)
            },
            endpoint=self.output_vocab_file  # Vocab file path in endpoint for downstream actors
        )
    
    def _extract_token_set(self, TOKENS: List[List[str]]) -> set:
        return {TOK for TOKEN in TOKENS for TOK in TOKEN}
            
    def _build_new_vocabulary(self, TOKENS:List[List[str]]):
        """Build vocabulary from SMILES list."""
        TOKENSET = self._extract_token_set(TOKENS)
        self.VOCABULARY.update((["*", "^", "$"]) + sorted(TOKENSET))
        
        for ring in list(range(1, 7)):
            if str(ring) not in self.VOCABULARY.tokens():
                self.VOCABULARY.update([str(ring)])
            
        self._write_vocabulary()

    def _write_vocabulary(self) -> None:
        path = os.path.abspath(self.output_vocab_file)
        self.log(f"{'Saved' if not os.path.isfile(path) else 'Updated'} vocabulary to {path}")
        save_vocabulary(self.VOCABULARY, path)
        return None

    def get_vocabulary(self) -> Vocabulary:
        """
        Public API for downstream actors to access vocabulary instance.

        Returns:
            Vocabulary instance with all tokens from processed data

        Example:
            >>> td = self.get_actor(Steps.TOKENS)
            >>> vocab = td.get_vocabulary()
            >>> if 'Br' in vocab:
            >>>     print("Bromine token present")
        """
        return self.VOCABULARY

