"""
Vocabulary and tokenization utilities for SMILES processing.

This module provides reusable vocabulary management and SMILES tokenization
functionality that can be used across different actors.
"""

import json
import os
import re
from typing import Optional, Set, List
import numpy as np


class Vocabulary:
    """Stores the tokens and their conversion to integer indices."""

    def __init__(
        self, tokens=None, starting_id=0, pad_token=0, bos_token=1, eos_token=2, unk_token=None
    ):
        self._tokens = {}
        self._current_id = starting_id

        self.pad_token = pad_token
        self.bos_token = bos_token
        self.eos_token = eos_token
        self.unk_token = unk_token

        if tokens:
            for token, idx in tokens.items():
                self._add(token, idx)
                self._current_id = max(self._current_id, idx + 1)

    def __getitem__(self, token_or_id):
        return self._tokens[token_or_id]

    def add(self, token):
        """
        Adds a token to the vocabulary.
        :param token: Token to add.
        :return: The id assigned to the token. If the token was already there,
                 the id of that token is returned instead.
        """
        if not isinstance(token, str):
            raise TypeError("Token is not a string")
        if token in self:
            return self[token]
        self._add(token, self._current_id)
        self._current_id += 1
        return self._current_id - 1

    def update(self, tokens):
        """Adds many tokens."""
        return [self.add(token) for token in tokens]

    def __delitem__(self, token_or_id):
        """
        Deletes a (token, id) tuple, given a token or an id.
        :param token_or_id: A token or an id.
        :return:
        """
        other_val = self._tokens[token_or_id]
        del self._tokens[other_val]
        del self._tokens[token_or_id]

    def __contains__(self, token_or_id):
        """
        Checks whether a token is contained in the vocabulary.
        :param token_or_id: token or an id to check
        :return : True if it is contained, otherwise False.
        """
        return token_or_id in self._tokens

    def __eq__(self, other_vocabulary):
        """
        Compares two vocabularies.
        :param other_vocabulary: Other vocabulary to be checked.
        :return: True if they are the same.
        """
        return self._tokens == other_vocabulary._tokens

    def __len__(self):
        """
        Calculates the length (number of tokens) of the vocabulary.
        :return : The number of tokens.
        """
        return len(self._tokens) // 2

    def encode(self, tokens):
        """
        Encodes a list of tokens as integer indices.

        Converts tokens to their vocabulary indices, which can be used as input to
        embedding layers in neural networks (e.g., LSTM, transformers).

        :param tokens: Tokens to encode.
        :return: A numpy array of integer indices. Unknown tokens are either
                 mapped to unk_token (if set) or filtered out (if unk_token is None).
        """
        encoded_indices = np.zeros(len(tokens), dtype=np.float32)
        keep_mask = np.ones_like(tokens, dtype=bool)
        for i, token in enumerate(tokens):
            if token not in self._tokens:
                if hasattr(self, "unk_token") and (self.unk_token is not None):
                    unk_symbol = self[self.unk_token]
                    encoded_indices[i] = self._tokens[unk_symbol]
                else:
                    keep_mask[i] = False
            else:
                encoded_indices[i] = self._tokens[token]
        return encoded_indices[keep_mask]

    def decode(self, encoded_indices):
        """
        Decodes integer indices back to a list of tokens.

        Converts vocabulary indices back to their corresponding token strings.

        :param encoded_indices: A numpy array or list of integer indices.
        :return: A list of decoded token strings.
        """
        tokens = []
        for idx in encoded_indices:
            tokens.append(self[idx])
        return tokens

    def _add(self, token, idx):
        if idx not in self._tokens:
            self._tokens[token] = idx
            self._tokens[idx] = token
        else:
            raise ValueError(f"Index {idx} already present in vocabulary")

    def tokens(self):
        """Returns the tokens from the vocabulary"""
        return [t for t in self._tokens if isinstance(t, str)]

    def word2idx(self):
        return {k: self._tokens[k] for k in self._tokens if isinstance(k, str)}

    def get_dictionary(self):
        return {
            "tokens": self.word2idx(),
            "pad_token": getattr(self, "pad_token", 0),
            "bos_token": getattr(self, "bos_token", 1),
            "eos_token": getattr(self, "eos_token", 2),
            "unk_token": getattr(self, "unk_token", None),
        }

    @classmethod
    def load_from_dictionary(cls, dictionary: dict):
        vocabulary = cls()
        for k, i in dictionary["tokens"].items():
            vocabulary._add(str(k), int(i))
        vocabulary.pad_token = dictionary["pad_token"]
        vocabulary.bos_token = dictionary["bos_token"]
        vocabulary.eos_token = dictionary["eos_token"]
        vocabulary.unk_token = dictionary["unk_token"]
        return vocabulary


class SMILESTokenizer:
    """Deals with the tokenization and untokenization of SMILES."""

    REGEXPS = {
        "brackets": re.compile(r"(\[[^\]]*\])"),
        "2_ring_nums": re.compile(r"(%\d{2})"),
        "brcl": re.compile(r"(Br|Cl)"),
    }
    REGEXP_ORDER = ["brackets", "2_ring_nums", "brcl"]

    def tokenize(self, smiles: str, with_begin_and_end: bool = True) -> List[str]:
        """
        Tokenizes a SMILES string.
        :param smiles: A SMILES string.
        :param with_begin_and_end : Whether to add begin and end tags.
        :return: A list with the tokens.
        """
        if smiles is None or smiles == "":
            return []

        def split_by(smiles_string, regexps):
            if not regexps:
                return list(smiles_string)
            regexp = self.REGEXPS[regexps[0]]
            splitted = regexp.split(smiles_string)
            tokens = []
            for i, split in enumerate(splitted):
                if i % 2 == 0:
                    tokens += split_by(split, regexps[1:])
                else:
                    tokens.append(split)
            return tokens

        tokens = split_by(smiles, self.REGEXP_ORDER)
        if with_begin_and_end:
            tokens = ["^"] + tokens + ["$"]
        return tokens

    def untokenize(self, tokens):
        """
        Untokenizes a SMILES string.
        :param tokens: List of tokens.
        :return : A SMILES string.
        """
        smi = ""
        for token in tokens:
            if token == "$":
                break
            if token != "^":
                smi += token
        return smi


# ==================== Vocabulary Management Utilities ====================


def load_vocabulary(vocab_file: str) -> Vocabulary:
    """
    Load vocabulary from JSON file.

    Args:
        vocab_file: Path to vocabulary JSON file

    Returns:
        Loaded Vocabulary instance

    Raises:
        FileNotFoundError: If vocabulary file doesn't exist
        ValueError: If vocabulary file is invalid
    """
    if not os.path.isfile(vocab_file):
        raise FileNotFoundError(f"Vocabulary file not found: {vocab_file}")

    with open(vocab_file, 'r') as f:
        vocab_dict = json.load(f)

    return Vocabulary.load_from_dictionary(vocab_dict)


def save_vocabulary(vocabulary: Vocabulary, vocab_file: str) -> None:
    """
    Save vocabulary to JSON file.

    Args:
        vocabulary: Vocabulary instance to save
        vocab_file: Path to save vocabulary JSON file
    """
    os.makedirs(os.path.dirname(vocab_file), exist_ok=True)
    with open(vocab_file, 'w') as f:
        json.dump(vocabulary.get_dictionary(), f, indent=2)


def create_vocabulary_from_tokens(
    tokens: Set[str],
    include_special: bool = True,
    include_ring_nums: bool = True
) -> Vocabulary:
    """
    Create a new vocabulary from a set of tokens.

    Args:
        tokens: Set of token strings
        include_special: Whether to include special tokens (^, $, *)
        include_ring_nums: Whether to ensure ring numbers 1-6 are included

    Returns:
        New Vocabulary instance
    """
    vocab = Vocabulary()

    # Add special tokens if requested
    if include_special:
        vocab.update(["*", "^", "$"])

    # Add provided tokens
    vocab.update(sorted(tokens))

    # Ensure ring numbers are included if requested
    if include_ring_nums:
        for ring in range(1, 7):
            if str(ring) not in vocab:
                vocab.update([str(ring)])

    return vocab


def extract_tokens_from_dataframe(df, tokens_column: str = 'tokens') -> Set[str]:
    """
    Extract unique tokens from a DataFrame's tokens column.

    Args:
        df: DataFrame with tokens column
        tokens_column: Name of column containing token lists

    Returns:
        Set of unique tokens across all molecules
    """
    unique_tokens = set()
    for tokens in df[tokens_column]:
        if isinstance(tokens, list):
            unique_tokens.update(tokens)
    return unique_tokens
