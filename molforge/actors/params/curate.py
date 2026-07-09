from typing import List, Union, Dict, Any, Literal
from dataclasses import dataclass, field

import os, shutil

from .base import BaseParams

@dataclass
class CurateMolParams(BaseParams):
    """Configuration parameters for molecule curation component."""

    mol_steps: List[str] = field(default_factory=lambda: [
        'desalt',
        'removeIsotope',
        'removeHs',
        'tautomers',
        'neutralize',
        'sanitize',
        'handleStereo',
        'computeProps',
    ])
    """Ordered list of mol-level curation steps to apply."""
    # smiles steps are handled separately, and at the end of the pipe
    smiles_steps: List[str] = field(default_factory=lambda: [
        'canonical',
        'kekulize',
    ])
    """SMILES-level steps, handled separately and applied at the end of the pipe."""

    neutralize_policy: Literal['keep', 'remove'] = 'keep'
    """'keep' returns already-neutral molecules (net formal charge 0) unchanged without running neutralization reactions; 'remove' always runs them (also neutralizing e.g. zwitterions). Charged molecules are neutralized under both; those that cannot reach zero charge fail the step."""
    desalt_policy: Literal['keep', 'remove'] = 'keep'
    """'keep' strips salts from multi-fragment molecules and keeps the desalted molecule; 'remove' drops any multi-fragment molecule without attempting desalting. Single-fragment molecules pass unchanged."""
    brute_force_desalt: bool = False
    """Use brute-force desalting to strip all fragments but the largest."""
    max_tautomers: int = 512
    """Maximum number of tautomers to enumerate per molecule."""

    step_timeout: int = 60
    """Per-step timeout in seconds per molecule (mainly for very large/complex molecules and CIP assignment)."""

    SMILES_column: str = 'canonical_smiles'
    """DataFrame column containing SMILES strings."""
    dropna: bool = True
    """Drop molecules that failed curation, keeping only rows where curation_success is True."""
    check_duplicates: bool = True
    """Check for and remove duplicate curated SMILES (post-curation)."""
    duplicates_policy: Literal['first', 'last', False] = 'first'
    """Which duplicate to keep when dropping duplicate curated SMILES; passed as pandas drop_duplicates `keep` ('first', 'last', or False). Only applied in the fallback path when no ChEMBLCurator actor is available."""

    stereo_policy: Literal['keep', 'remove', 'assign', 'enumerate'] = 'assign'
    """How to handle stereochemistry: keep, remove, assign, or enumerate stereoisomers."""
    assign_policy: Literal['first', 'random', 'lowest'] = 'random'
    """Which stereoisomer to assign when stereo_policy is 'assign'."""

    random_seed: int = 42
    """Random seed for reproducible stereoisomer assignment/enumeration."""
    try_embedding: bool = False
    """Passed to RDKit StereoEnumerationOptions.tryEmbedding; when True, enumerated stereoisomers that cannot be 3D-embedded are discarded."""
    only_unassigned: bool = True
    """Only enumerate/assign stereochemistry at unassigned stereocenters."""
    only_unique: bool = True
    """Passed to RDKit's StereoEnumerationOptions(unique=...); when True, duplicate enumerated stereoisomers are filtered out (applies in stereo_policy='enumerate')."""
    max_isomers: int = 32
    """Maximum number of stereoisomers to enumerate per molecule."""
    
    _VALID_MOL_STEPS: List[str] = field(default_factory=lambda: [
        'desalt', 'removeIsotope', 'removeHs', 'tautomers', 'neutralize', 'sanitize',
        'computeProps', 'handleStereo',
    ])
    _VALID_SMILES_STEPS: List[str] = field(default_factory=lambda: [
        'canonical', 'kekulize',
    ])
    _BASIC_POLICIES: List[str] = field(default_factory=lambda: ['keep', 'remove'])
    _STEREO_POLICIES: List[str] = field(default_factory=lambda: ['assign', 'enumerate'])
    _ASSIGN_POLICIES: List[str] = field(default_factory=lambda: ['first', 'random', 'lowest'])
    
        
    def _validate_params(self) -> None:
        self._validate_policy('desalt_policy', self.desalt_policy, self._BASIC_POLICIES)
        self._validate_policy('neutralize_policy', self.neutralize_policy, self._BASIC_POLICIES)
        self._validate_policy('duplicates_policy', self.duplicates_policy, ['first', 'last', False])

        self._validate_policy('stereo_policy', self.stereo_policy, self._BASIC_POLICIES + self._STEREO_POLICIES)
        self._validate_policy('assign_policy', self.assign_policy, self._ASSIGN_POLICIES)       
        
        self._validate_steps()

    def _validate_steps(self) -> None:
        """Validate specified steps."""

        for step in self.mol_steps:
            if step not in self._VALID_MOL_STEPS:
                raise ValueError(f"Invalid mol_step '{step}'. Valid mol steps are: {self._VALID_MOL_STEPS}")
        
        for step in self.smiles_steps:
            if step not in self._VALID_SMILES_STEPS:
                raise ValueError(f"Invalid smiles_step '{step}'. Valid SMILES steps are: {self._VALID_SMILES_STEPS}")
        
        if len(self.mol_steps) != len(set(self.mol_steps)):
            duplicates = [step for step in set(self.mol_steps) if self.mol_steps.count(step) > 1]
            raise ValueError(f"Duplicate mol_steps found: {duplicates}")
            
        if len(self.smiles_steps) != len(set(self.smiles_steps)):
            duplicates = [step for step in set(self.smiles_steps) if self.smiles_steps.count(step) > 1]
            raise ValueError(f"Duplicate smiles_steps found: {duplicates}")
        
        if set(self.mol_steps) & set(self.smiles_steps):
            raise ValueError(f"Steps cannot be in both mol_steps and smiles_steps: {set(self.mol_steps) & set(self.smiles_steps)}")