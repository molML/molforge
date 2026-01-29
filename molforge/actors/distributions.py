"""
Molecular distribution curation module for filtering based on chemical properties.

This module provides comprehensive distribution-based curation including:
- Configurable property thresholds (absolute or statistical)
- Token-based vocabulary filtering
- Parallel processing for large datasets
"""

from typing import List, Optional, Dict, Any, Tuple
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import gc
import os

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

from .params.distributions import PropertyThreshold
from .base import BaseActor
from .protocol import ActorOutput
from ..configuration.steps import Steps

from ..utils.actortools.multiprocess import multiprocess_worker, calculate_chunk_params
from ..utils.actortools.vocabulary import (
    create_vocabulary_from_tokens,
    extract_tokens_from_dataframe,
    save_vocabulary
)
from ..utils.constants import DEFAULT_MP_THRESHOLD, MAX_CHUNK_SIZE


class CurateDistribution(BaseActor):
    """Distribution-based molecular curation."""
    __step_name__ = Steps.DISTRIBUTIONS
    """
    Distribution-based molecular curation.

    Filters molecules based on chemical property distributions using either:
    - Absolute thresholds (min/max values)
    - Statistical thresholds (mean ± n*std)
    - Quantile thresholds (percentiles)
    - Token vocabulary filtering
    """

    # Class constants (imported from utils.constants)
    DEFAULT_MP_THRESHOLD = DEFAULT_MP_THRESHOLD
    MAX_CHUNK_SIZE = MAX_CHUNK_SIZE
    DEFAULT_PROGRESS_INTERVAL = 10
    MAX_PLOTTING_SAMPLE_SIZE = 3_000_000

    # Property computation patterns
    SMARTS_CHAINS = [Chem.MolFromSmarts("-".join(["[CR0H2]"] * i)) for i in range(1, 11)]

    def __post_init__(self):
        """Post-initialization setup for distribution curation."""
        self._statistics_cache = {}
        self._property_mapping = self._build_property_mapping()
        self._properties_to_compute = self._get_properties_to_compute()

        self.color_primary = "#769ff5"    # Vista Blue
        self.color_threshold = "#db5881"  # Blush Red
        self.color_filtered = '#D3D3D3'   # Light gray

        self._log_configuration()

    @property
    def required_columns(self) -> List[str]:
        """Required input columns."""
        return [self.SMILES_column]

    @property
    def output_columns(self) -> List[str]:
        """Output columns depend on computed properties."""
        return self._properties_to_compute

    @property
    def distributions_dir(self) -> str:
        """Path to distributions directory in current run directory."""
        if not self.plot_distributions:
            return None
        dir_path = self._get_run_path("distributions")
        os.makedirs(dir_path, exist_ok=True)
        return dir_path

    @property
    def curated_vocab_file(self) -> str:
        """Path to curated vocabulary file in current run directory."""
        return self._get_run_path("vocab_curated.json")

    @property
    def curation_results_file(self) -> str:
        """Path to curation results JSON file in current run directory."""
        return self._get_run_path("curation_results.json")

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Perform distribution-based curation on a DataFrame.

        Args:
            data: Input DataFrame with SMILES/tokens

        Returns:
            Curated DataFrame

        Raises:
            ValueError: If SMILES column not found
        """
        if self.SMILES_column not in data.columns:
            raise ValueError(
                f"'{self.SMILES_column}' not in dataframe. "
                f"Available columns: {list(data.columns)}"
            )

        data = data.copy()

        # Handle tokenization
        data = self._ensure_tokenization(data)

        # Compute properties if needed
        if self.compute_properties:
            data = self._compute_all_properties(data)

        # Apply distribution filters
        data = self._apply_distribution_filters(data)

        # Apply token filters if applicable
        if self.curate_tokens and self.tokens_column in data.columns:
            data = self._apply_token_filters(data)

        self.log(f"Distribution curation complete: {len(data)} molecules retained.")

        # Save curation results
        self._save_curation_results()

        return data

    def _create_output(self, data: pd.DataFrame) -> ActorOutput:
        """Create output with metadata and curated vocabulary endpoint."""
        metadata = {
            'n_curated': len(data),
            'properties_computed': self._properties_to_compute,
            'distributions_dir': self.distributions_dir
        }

        # Create curated vocabulary if token filtering was performed
        endpoint = None
        if self.curate_tokens and self.tokens_column in data.columns:
            curated_vocab_file = self._create_curated_vocabulary(data)
            if curated_vocab_file:
                endpoint = curated_vocab_file
                metadata['curated_vocab_file'] = curated_vocab_file

        return ActorOutput(
            data=data,
            success=True,
            metadata=metadata,
            endpoint=endpoint
        )
    
    # ==================== Tokenization ====================
    
    def _ensure_tokenization(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure tokens column exists, either by reusing existing or skipping token filters."""

        if self.tokens_column in df.columns:
            self.log(f"Using existing '{self.tokens_column}' column.")
            if 'seqlen' not in df.columns:
                df['seqlen'] = df[self.tokens_column].apply(
                    lambda t: len(t) if isinstance(t, list) else 0
                )
            return df

        if self.curate_tokens or self.token_frequency_threshold is not None:
            self.log(
                f"No '{self.tokens_column}' column found. "
                f"Token filtering will be skipped. "
                f"Run TokenizeData actor before CurateDistribution to enable token filtering.",
                level='WARNING'
            )

        return df

    # ==================== Property Computation ====================
    
    def _compute_all_properties(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute all molecular properties."""
        df = df.copy()
        
        # Separate molecular properties from token-based properties
        molecular_props = [p for p in self._properties_to_compute 
                          if p in self._property_mapping]
        token_props = [p for p in self._properties_to_compute 
                      if p in ['num_tokens', 'tokens_atom_ratio']]
        
        # Compute molecular properties
        if molecular_props:
            if len(df) < self.DEFAULT_MP_THRESHOLD:
                self.log(f"Computing properties for {len(df):,} molecules (single process).")
                df = self._compute_properties_single_process(df, molecular_props)
            else:
                self.log(f"Computing properties for {len(df):,} molecules (multiprocessing).")
                df = self._compute_properties_multiprocess(df, molecular_props)
        
        # Compute token-based properties
        if token_props:
            df = self._compute_token_properties(df, token_props)
        
        return df
    
    def _compute_properties_single_process(self, df: pd.DataFrame, 
                                          properties: List[str]) -> pd.DataFrame:
        """Compute molecular properties using single process."""
        smiles_list = df[self.SMILES_column].tolist()
        results = [self._compute_molecular_properties(smi, properties) 
                  for smi in smiles_list]
        return self._merge_property_results(df, results)

    def _compute_properties_multiprocess(self, df: pd.DataFrame, 
                                        properties: List[str]) -> pd.DataFrame:
        """Compute molecular properties using multiprocessing."""
        import multiprocessing as mp
        
        n_processes, chunk_size, n_chunks = calculate_chunk_params(len(df), self.MAX_CHUNK_SIZE)
        
        self.log(
            f"Processing with {n_processes} processes "
            f"({n_chunks} chunks of {chunk_size:,})."
        )
        
        # Prepare chunks
        smiles_list = df[self.SMILES_column].tolist()
        smiles_chunks = [smiles_list[i:i + chunk_size] 
                        for i in range(0, len(smiles_list), chunk_size)]
        
        # Process chunks in parallel
        results = self._process_chunks_parallel(smiles_chunks, properties, n_processes)
        return self._merge_property_results(df, results)

    def _process_chunks_parallel(self, chunks: List[List[str]], properties: List[str],
                                n_processes: int) -> List[Dict]:
        """Process chunks in parallel with progress tracking."""
        total_chunks = len(chunks)
        chunk_results = [None] * total_chunks
        
        start_time = time.time()
        last_progress_time = start_time
        completed_chunks = 0
        
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            worker_func = partial(
                multiprocess_worker,
                actor_class=self.__class__,
                actor_method='_compute_molecular_properties',
                actor_params=self._params,
                logger=self.logger,
                method_kwargs={'properties': properties},
            )
            
            future_to_chunk = {executor.submit(worker_func, chunk): i 
                              for i, chunk in enumerate(chunks)}
            
            for future in as_completed(future_to_chunk, timeout=3600):
                chunk_idx = future_to_chunk[future]
                try:
                    chunk_results[chunk_idx] = future.result(timeout=300)
                    completed_chunks += 1
                except Exception as e:
                    self.log(f"Chunk {chunk_idx} failed: {e}", level="ERROR")
                    chunk_results[chunk_idx] = [{prop: None for prop in properties}] * len(chunks[chunk_idx])
                    
                # Report progress every 10 seconds or on completion
                current_time = time.time()
                if current_time - last_progress_time >= self.DEFAULT_PROGRESS_INTERVAL or completed_chunks == total_chunks:
                    self._report_progress(completed_chunks, total_chunks, start_time, len(chunks[0]))
                    last_progress_time = current_time
        
        # Flatten results
        all_results = [result for chunk_result in chunk_results for result in chunk_result]
        gc.collect()
        
        elapsed = time.time() - start_time
        rate = len(all_results) / elapsed if elapsed > 0 else 0
        self.log(f"Completed in {self._format_duration(elapsed)} | Rate: {rate:.0f} mol/s")
        
        return all_results

    def _compute_molecular_properties(self, smiles: str, 
                                     properties: Optional[List[str]] = None) -> Dict[str, Any]:
        """Compute specified molecular properties for a single SMILES."""
        if properties is None:
            properties = [p for p in self._properties_to_compute if p in self._property_mapping]
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {prop: None for prop in properties}
        
        result = {}
        for prop in properties:
            try:
                result[prop] = self._property_mapping[prop](mol)
            except Exception as e:
                self.log(f"Error computing {prop} for {smiles}: {e}", level='DEBUG')
                result[prop] = None
        
        return result

    def _compute_token_properties(self, df: pd.DataFrame, properties: List[str]) -> pd.DataFrame:
        """Compute token-based properties from the tokens column."""
        if self.tokens_column not in df.columns:
            self.log(
                f"Cannot compute token properties: '{self.tokens_column}' column not found.",
                level='WARNING'
            )
            for prop in properties:
                df[prop] = None
            return df
        
        self.log(f"Computing token-based properties from '{self.tokens_column}' column.")
        
        # Compute num_tokens
        if 'num_tokens' in properties:
            df['num_tokens'] = df.get('seqlen', df[self.tokens_column].apply(
                lambda t: len(t) if isinstance(t, list) else None
            ))
        
        # Compute tokens_atom_ratio
        if 'tokens_atom_ratio' in properties:
            if 'num_atoms' not in df.columns:
                self.log(
                    "Cannot compute tokens_atom_ratio: num_atoms not available.",
                    level='WARNING'
                )
                df['tokens_atom_ratio'] = None
            else:
                num_tokens = df.get('num_tokens', df.get('seqlen', 
                    df[self.tokens_column].apply(lambda t: len(t) if isinstance(t, list) else None)
                ))
                df['tokens_atom_ratio'] = num_tokens / df['num_atoms'].replace(0, np.nan)
        
        return df
    
    def _merge_property_results(self, df: pd.DataFrame, results: List[Dict]) -> pd.DataFrame:
        """Merge computed properties back into DataFrame."""
        results_df = pd.DataFrame(results)
        for col in results_df.columns:
            if col not in df.columns:
                df[col] = results_df[col].values
        return df
    
    def _build_property_mapping(self) -> Dict[str, callable]:
        """Build mapping of property names to computation functions."""
        return {
            'num_atoms': lambda mol: mol.GetNumHeavyAtoms(),
            'num_rings': lambda mol: len(list(Chem.GetSSSR(mol))),
            'size_largest_ring': lambda mol: max([0] + [len(ring) for ring in mol.GetRingInfo().AtomRings()]),
            'c_atom_ratio': lambda mol: self._compute_ratio(sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6), mol.GetNumHeavyAtoms()),
            'longest_aliph_carbon': lambda mol: self._compute_longest_chain(mol),
            'molecular_weight': lambda mol: Descriptors.MolWt(mol),
            'logp': lambda mol: Descriptors.MolLogP(mol),
            'tpsa': lambda mol: Descriptors.TPSA(mol),
            'num_rotatable_bonds': lambda mol: Descriptors.NumRotatableBonds(mol),
            'num_h_donors': lambda mol: Descriptors.NumHDonors(mol),
            'num_h_acceptors': lambda mol: Descriptors.NumHAcceptors(mol),
            'fsp3': lambda mol: Descriptors.FractionCSP3(mol),
            'num_aromatic_rings': lambda mol: Descriptors.NumAromaticRings(mol),
            'num_stereocenters': lambda mol: self._compute_stereocenters(mol),
            'num_heteroatoms': lambda mol: sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6]),
            'heteroatom_ratio': lambda mol: self._compute_ratio(sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6]), mol.GetNumHeavyAtoms()),
        }

    @staticmethod
    def _compute_ratio(numerator: int, denominator: int) -> float:
        """Safely compute ratio, returning 0 if denominator is 0."""
        return numerator / denominator if denominator > 0 else 0.0
    
    @staticmethod
    def _compute_stereocenters(mol: Chem.Mol) -> int:
        """Compute number of stereocenters."""
        try:
            return len(Chem.FindMolChiralCenters(
                mol, includeUnassigned=True, useLegacyImplementation=False, includeCIP=False
            ))
        except Exception:
            return 0

    def _compute_longest_chain(self, mol: Chem.Mol) -> int:
        """Compute longest aliphatic carbon chain."""
        for i, chain in enumerate(self.SMARTS_CHAINS, 1):
            if not mol.HasSubstructMatch(chain):
                return i - 1
        return len(self.SMARTS_CHAINS)
    
    def _get_properties_to_compute(self) -> List[str]:
        """Determine which properties should be computed based on configuration."""
        if self.properties == 'all':
            return self._params._VALID_PROPERTIES.copy()
        
        if isinstance(self.properties, list):
            return self.properties
        
        # None: compute only properties with defined thresholds
        properties_with_thresholds = []
        for prop in self._params._VALID_PROPERTIES:
            threshold = self.thresholds.get(prop)
            if threshold and self._has_active_threshold(threshold):
                properties_with_thresholds.append(prop)
        
        return properties_with_thresholds
    
    @staticmethod
    def _has_active_threshold(threshold: PropertyThreshold) -> bool:
        """Check if threshold has any active bounds."""
        return any([
            threshold.min_value is not None,
            threshold.max_value is not None,
            threshold.statistical_lower is not None,
            threshold.statistical_upper is not None,
            threshold.quantile_lower is not None,
            threshold.quantile_upper is not None
        ])
    
    # ==================== Distribution Filtering ====================
    
    def _apply_distribution_filters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply all distribution-based filters."""
        initial_count = len(df)
        df_original = df.copy()
        
        # Compute statistics for statistical thresholds
        self._compute_statistics(df)
        
        # Compute PCA if requested
        if self.perform_pca:
            df = self._compute_pca_coordinates(df)

        # Show initial distribution statistics
        self._show_distribution_statistics(df, "INITIAL DISTRIBUTION")
        
        # Apply filters independently and track results
        filter_results = []
        combined_mask = pd.Series([True] * len(df), index=df.index)

        for prop_name, threshold_config in self.thresholds.items():
            if prop_name in df.columns and self._has_active_threshold(threshold_config):
                prop_mask = self._get_property_mask(df_original, prop_name, threshold_config)
                combined_mask &= prop_mask

                passed = prop_mask.sum()
                removed = len(df_original) - passed

                filter_results.append({
                    'property': prop_name,
                    'passed': passed,
                    'removed': removed,
                    'pass_rate': 100 * passed / len(df_original) if len(df_original) > 0 else 0
                })
        
        # Apply combined mask
        df = df[combined_mask]
    
        # Show filter summary
        if filter_results:
            self._show_filter_summary(filter_results, initial_count)
        
        removed = initial_count - len(df)
        self.log(
            f"Distribution filters: {len(df)}/{initial_count} molecules retained "
            f"({removed} removed, {100*len(df)/initial_count:.1f}% overall pass rate)."
        )
        
        return df
        
    def _compute_statistics(self, df: pd.DataFrame) -> None:
        """Compute and cache statistical measures for properties."""
        from scipy.stats import skew, kurtosis
        
        for prop in self._properties_to_compute:
            if prop in df.columns:
                values = df[prop].dropna()
                if len(values) > 0:
                    n_digits = 2

                    stats = {
                        'count': len(values),
                        'mean': round(values.mean(), n_digits),
                        'std': round(values.std(), n_digits),
                        'min': round(values.min(), n_digits),
                        'max': round(values.max(), n_digits),
                        'median': round(values.median(), n_digits),
                        'skewness': round(skew(values), n_digits),
                        'kurtosis': round(kurtosis(values), n_digits),
                    }
                    
                    threshold = self.thresholds.get(prop)
                    threshold_str = self._format_threshold(threshold, stats.get('mean'), stats.get('std'), values)

                    stats['thresholds'] = threshold_str
                    self._statistics_cache[prop] = stats
                        
    def _get_property_mask(self, df: pd.DataFrame, property_name: str, 
                        threshold: PropertyThreshold) -> pd.Series:
        """Get boolean mask for a single property filter."""
        mask = pd.Series([True] * len(df), index=df.index)
        stats = self._statistics_cache.get(property_name, {})
        values = df[property_name]

        # Compute bounds
        lower_bound, upper_bound = self._get_threshold_bounds(property_name, threshold,
            stats.get('mean'), stats.get('std'), values
        )
        
        # Apply bounds
        if lower_bound is not None:
            mask &= values >= lower_bound
        if upper_bound is not None:
            mask &= values <= upper_bound
        
        return mask

    @staticmethod
    def _compute_bound(absolute_value: Optional[float], statistical_value: Optional[float],
                      quantile_value: Optional[float], mean: Optional[float],
                      std: Optional[float], values: pd.Series) -> Optional[float]:
        """Compute threshold bound from absolute, statistical, or quantile specification."""
        if absolute_value is not None:
            return absolute_value
        
        if statistical_value is not None and mean is not None and std is not None:
            return mean + (statistical_value * std)
        
        if quantile_value is not None:
            return values.quantile(quantile_value)
        
        return None
    
    # ==================== Token Filtering ====================
    
    def _apply_token_filters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply token-based filtering with detailed logging."""
        if self.tokens_column not in df.columns:
            self.log("Token column not found. Skipping token filtering.", level='WARNING')
            return df

        # Get vocabulary from TokenizeData actor instance via context
        td_actor = self.get_actor(Steps.TOKENS)
        if not td_actor:
            self.log(
                "No TokenizeData actor found in context. Skipping token filtering. "
                "Run TokenizeData actor before CurateDistribution to enable token filtering.",
                level='WARNING'
            )
            return df

        vocabulary = td_actor.get_vocabulary()
        if not vocabulary or len(vocabulary) == 0:
            self.log(
                "TokenizeData vocabulary is empty. Skipping token filtering.",
                level='WARNING'
            )
            return df

        initial_count = len(df)

        # Filter by token frequency (doesn't require vocabulary object)
        if self.token_frequency_threshold is not None:
            df, rare_tokens = self._filter_by_token_frequency(df)

        removed = initial_count - len(df)
        if removed > 0:
            self.log(
                f"Token filtering: {len(df)}/{initial_count} molecules retained "
                f"({removed} removed, {100*len(df)/initial_count:.1f}% pass rate)."
            )

        return df

    def _create_curated_vocabulary(self, df: pd.DataFrame) -> Optional[str]:
        """
        Create and save curated vocabulary from filtered molecules.

        Args:
            df: DataFrame with tokens column

        Returns:
            Path to curated vocabulary file, or None if creation failed
        """
        try:
            # Extract unique tokens from curated dataset
            curated_tokens = extract_tokens_from_dataframe(df, self.tokens_column)

            if not curated_tokens:
                self.log("No tokens found for curated vocabulary", level='WARNING')
                return None

            # Create vocabulary from tokens
            curated_vocab = create_vocabulary_from_tokens(
                curated_tokens,
                include_special=True,
                include_ring_nums=True
            )

            # Save curated vocabulary
            curated_vocab_file = self.curated_vocab_file
            save_vocabulary(curated_vocab, curated_vocab_file)

            # Log statistics
            td_actor = self.get_actor(Steps.TOKENS)
            if td_actor:
                original_vocab = td_actor.get_vocabulary()
                original_size = len(original_vocab)
                curated_size = len(curated_vocab)
                self.log(
                    f"Curated vocabulary saved to {curated_vocab_file}\n"
                    f"  Original: {original_size} tokens\n"
                    f"  Curated: {curated_size} tokens\n"
                    f"  Removed: {original_size - curated_size} tokens"
                )
            else:
                self.log(f"Curated vocabulary saved to {curated_vocab_file} ({len(curated_vocab)} tokens)")

            return curated_vocab_file

        except Exception as e:
            self.log(f"Failed to create curated vocabulary: {e}", level='ERROR')
            return None

    def _filter_by_token_frequency(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, set]:
        """Filter molecules with rare tokens and create distribution plot."""
        # Compute token frequencies - each molecule should contribute once per UNIQUE token
        all_tokens = []
        for tokens in df[self.tokens_column]:
            if isinstance(tokens, list):
                all_tokens.extend(set(tokens))  # Use unique tokens per molecule
        
        token_counts = pd.Series(all_tokens).value_counts()
        total_molecules = len(df)
        token_freq_df = pd.DataFrame({
            'token': token_counts.index,
            'count': token_counts.values,
            'percent': (token_counts.values / total_molecules) * 100
        }).sort_values('percent', ascending=False).reset_index(drop=True)
        
        # Identify rare tokens
        rare_tokens = set(token_freq_df[token_freq_df['percent'] < self.token_frequency_threshold]['token'])
        
        if rare_tokens:
            self.log(
                f"Found {len(rare_tokens)} rare tokens below {self.token_frequency_threshold}% threshold: "
                f"{sorted(rare_tokens)}"
            )
            
            def has_no_rare_tokens(tokens):
                if not isinstance(tokens, list):
                    return False
                return not any(token in rare_tokens for token in tokens)
            
            mask = df[self.tokens_column].apply(has_no_rare_tokens)
            removed = (~mask).sum()
            
            if removed > 0:
                self.log(f"Removed {removed} molecules with rare tokens.")
            
            # Create token distribution plot
            if self.plot_distributions and self.distributions_dir:
                self._plot_token_distribution(token_freq_df, self.token_frequency_threshold)
            
            return df[mask], rare_tokens
        
        return df, set()
    
    
    # ==================== Statistics and Reporting ====================
    
    def _show_distribution_statistics(self, df: pd.DataFrame, label: str) -> None:
        """Display distribution statistics for all computed properties with their thresholds."""
        from scipy.stats import normaltest
        
        stats_data = []
        for prop in self._properties_to_compute:
            if prop in df.columns:
                values = df[prop].dropna()
                if len(values) > 0:
                    stats = self._statistics_cache.get(prop, {})
                    
                    # Normality test
                    if len(values) >= 8:
                        try:
                            _, p_value = normaltest(values)
                            normality_str = "✓" if p_value >= 0.05 else "✗"
                            pvalue_str = f"{p_value:.2e}"
                        except Exception:
                            normality_str = "—"
                            pvalue_str = "N/A"
                    else:
                        normality_str = "—"
                        pvalue_str = "N/A"
                    
                    # Practical normality (skewness and kurtosis)
                    skew_val = stats.get('skewness')
                    kurt_val = stats.get('kurtosis')
                    is_approx_normal = (abs(skew_val) < 1 and abs(kurt_val) < 2) if (
                        skew_val is not None and kurt_val is not None) else False
                    approx_normal_str = "✓" if is_approx_normal else "✗"

                    stats_data.append({
                        'property': prop,
                        'count': stats.get('count', len(values)),
                        'mean': stats.get('mean'),
                        'std': stats.get('std'),
                        'min': stats.get('min'),
                        'max': stats.get('max'),
                        'median': stats.get('median'),
                        'skew': skew_val,
                        'kurt': kurt_val,
                        'symmetric [*]': approx_normal_str,
                        'normal [**]': normality_str,
                        'p_value': pvalue_str,
                    })
        
        if stats_data:
            stats_df = pd.DataFrame(stats_data)
            self.log(f"{label}\n{stats_df.to_markdown(index=False)}")
            self.log(
                "[*] Practical normality: Skewness and Kurtosis tests. "
                "✓ = approximately normal (|skew|<1, |kurt|<2); ✗ = notably non-normal."
            )
            self.log(
                "[**] Normality test: D'Agostino-Pearson test. "
                "p < 0.05 suggests non-normal distribution."
            )

        # constrain sample size for very large input.
        sample_size = min(self.MAX_PLOTTING_SAMPLE_SIZE, len(df))
        
        # Create visualizations if enabled
        if self.plot_distributions and self.distributions_dir:
            self.log(f"Generating distribution plots with {sample_size:,} molecules...")
            for prop in self._properties_to_compute:
                if prop in df.columns:
                    try:
                        threshold = self.thresholds.get(prop)
                        stats = self._statistics_cache.get(prop, {})
                        self._plot_property_distribution(df, prop, threshold, stats, sample_size)
                    except Exception as e:
                        self.log(f"Failed to create plot for {prop}: {e}", level='WARNING')
        
        # Create PCA plot if enabled
        if self.perform_pca and self.distributions_dir:
            try:
                self._plot_pca(df, sample_size)
            except Exception as e:
                self.log(f"Failed to create PCA plot: {e}", level='WARNING')
                
    def _format_threshold(self, threshold: Optional[PropertyThreshold], 
                          mean: float, std: float, values: pd.Series = None) -> str:
        """Format threshold configuration as a readable string."""
        if threshold is None or not self._has_active_threshold(threshold):
            return "—"
        
        lower_parts = []
        upper_parts = []
        
        # Lower bound (priority: absolute > statistical > quantile)
        if threshold.min_value is not None:
            lower_parts.append(f"≥{threshold.min_value:.2f}")
        elif threshold.statistical_lower is not None and mean is not None and std is not None:
            computed_lower = mean + (threshold.statistical_lower * std)
            lower_parts.append(f"≥μ{threshold.statistical_lower:+.1f}σ ({computed_lower:.2f})")
        elif threshold.quantile_lower is not None and values is not None:
            computed_lower = values.quantile(threshold.quantile_lower)
            lower_parts.append(f"≥Q{threshold.quantile_lower:.3f} ({computed_lower:.2f})")
        
        # Upper bound
        if threshold.max_value is not None:
            upper_parts.append(f"≤{threshold.max_value:.2f}")
        elif threshold.statistical_upper is not None and mean is not None and std is not None:
            computed_upper = mean + (threshold.statistical_upper * std)
            upper_parts.append(f"≤μ{threshold.statistical_upper:+.1f}σ ({computed_upper:.2f})")
        elif threshold.quantile_upper is not None and values is not None:
            computed_upper = values.quantile(threshold.quantile_upper)
            upper_parts.append(f"≤Q{threshold.quantile_upper:.3f} ({computed_upper:.2f})")
        
        # Format with alignment
        if not lower_parts and not upper_parts:
            return "—"
        
        lower_str = ", ".join(lower_parts) if lower_parts else ""
        upper_str = ", ".join(upper_parts) if upper_parts else ""
        
        if lower_str and upper_str:
            return f"{lower_str:<20} , {upper_str:<20}"
        elif lower_str:
            return f"{lower_str:<20}"
        else:
            return f"{upper_str:<42}"
        
    def _show_filter_summary(self, filter_results: List[Dict], initial_count: int) -> None:
        """Display summary of filtering results per property."""
        for result in filter_results:
            stat = self._statistics_cache.get(result['property'], {})
            result['thresholds'] = stat.get('thresholds', '—')
        
        summary_df = pd.DataFrame(filter_results)
        summary_df['passed'] = summary_df['passed'].astype(str) + f"/{initial_count}"
        summary_df['pass_rate'] = summary_df['pass_rate'].round(1).astype(str) + "%"
        summary_df = summary_df.sort_values('removed', ascending=False)
        
        self.log(f"FILTER RESULTS\n{summary_df.to_markdown(index=False)}")
    
    # ==================== Plotting ====================
    
    def _plot_property_distribution(self, df: pd.DataFrame, prop: str, 
                                   threshold: PropertyThreshold, stats: dict, 
                                   sample_size: int) -> None:
        """Create Q-Q plot and histogram with threshold visualization."""
        import matplotlib.pyplot as plt
        from scipy import stats as scipy_stats
        
        values = df[prop].dropna()
        if len(values) == 0:
            return
        
        # Sample if too large
        if len(values) > sample_size:
            values = values.sample(sample_size, random_state=42)
        
        # Compute threshold bounds
        lower_bound, upper_bound = self._get_threshold_bounds(prop,
            threshold, stats.get('mean'), stats.get('std'), values
        )
        
        # Compute robust axis limits (1st to 99th percentile)
        x_min = values.quantile(0.01)
        x_max = values.quantile(0.99)
        
        # Extend limits slightly for visual padding
        x_range = x_max - x_min
        x_min -= 0.05 * x_range
        x_max += 0.05 * x_range
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        
        # Q-Q plot
        scipy_stats.probplot(values, dist="norm", plot=ax1)
        # access points through matplot
        points = ax1.get_lines()[0]
        points.set_marker('o')
        points.set_markerfacecolor(self.color_primary)
        points.set_markeredgecolor(self.color_primary)
        points.set_markersize(6.0)
        points.set_alpha(1.0)
        
        line = ax1.get_lines()[1]
        line.set_linewidth(1.0)
        line.set_color(self.color_threshold)

        ax1.set_title(f'Q-Q Plot: {prop}')
        ax1.grid(True, alpha=0.1)
        
        # Add threshold lines to Q-Q plot
        if lower_bound is not None or upper_bound is not None:
            ylim = ax1.get_ylim()
            xlim = ax1.get_xlim()  # Get xlim once and reuse
            
            if lower_bound is not None:
                ax1.axhline(y=lower_bound, color=self.color_threshold, linestyle='--', 
                        linewidth=1.0, label=f'Lower: {lower_bound:.2f}')
                # Shade filtered region below
                ax1.fill_between(xlim, ylim[0], lower_bound, 
                            color=self.color_filtered, alpha=0.2)
            
            if upper_bound is not None:
                ax1.axhline(y=upper_bound, color=self.color_threshold, linestyle='--', 
                        linewidth=1.0, label=f'Upper: {upper_bound:.2f}')
                # Shade filtered region above
                ax1.fill_between(xlim, upper_bound, ylim[1], 
                            color=self.color_filtered, alpha=0.2)
            
            # make sure limits were not moved by filling.
            ax1.set_ylim(ylim)
            ax1.set_xlim(xlim)
            ax1.legend(fontsize=8)
        
        # Histogram with normal overlay
        values_plot = values[(values >= x_min) & (values <= x_max)]
        ax2.hist(values_plot, bins=50, density=True, alpha=0.7, edgecolor='black', color=self.color_primary)
        mu, sigma = values.mean(), values.std()
        
        if sigma > 0:
            x = np.linspace(x_min, x_max, 100)
            ax2.plot(x, scipy_stats.norm.pdf(x, mu, sigma), color=self.color_threshold, ls='-', lw=1, label='Normal fit')
        
        # Add threshold lines to histogram
        if lower_bound is not None or upper_bound is not None:
            ylim = ax2.get_ylim()
            if lower_bound is not None and lower_bound >= x_min:
                ax2.axvline(x=lower_bound, color=self.color_threshold, linestyle='--', linewidth=1.0)
                # Shade filtered region
                ax2.axvspan(x_min, lower_bound, color=self.color_filtered, alpha=0.5)
            if upper_bound is not None and upper_bound <= x_max:
                ax2.axvline(x=upper_bound, color=self.color_threshold, linestyle='--', linewidth=1.0)
                # Shade filtered region
                ax2.axvspan(upper_bound, x_max, color=self.color_filtered, alpha=0.5)
                
            ax2.set_ylim(ylim)
        
        ax2.set_xlim(x_min, x_max)

        if sigma > 0:
            ax2.legend(loc='upper right')
        
        ax2.set_title(f'Distribution: {prop}')
        ax2.set_xlabel(prop)
        ax2.set_ylabel('Density')
        ax2.grid(True, alpha=0.1)
        
        # Add statistics text
        sk = scipy_stats.skew(values)
        kt = scipy_stats.kurtosis(values)
        n_outliers = len(values) - len(values_plot)
        stats_text = f'μ={mu:.2f}, σ={sigma:.2f}\nskew={sk:.2f}, kurt={kt:.2f}'
        
        if n_outliers > 0:
            stats_text += f'\n({n_outliers} outliers hidden)'
        
        ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes,
                verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        # Save figure
        output_path = os.path.join(self.distributions_dir, f'{prop}.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        self.log(f"  Saved plot: {output_path}", level='DEBUG')

    def _get_threshold_bounds(self, prop: str, threshold: Optional[PropertyThreshold],
                             mean: Optional[float], std: Optional[float],
                             values: pd.Series) -> Tuple[Optional[float], Optional[float]]:
        """Compute lower and upper threshold bounds for visualization."""
        if threshold is None or not self._has_active_threshold(threshold):
            return None, None
        
        lower_bound = self._compute_bound(
            threshold.min_value,
            threshold.statistical_lower,
            threshold.quantile_lower,
            mean, std, values
        )
        
        upper_bound = self._compute_bound(
            threshold.max_value,
            threshold.statistical_upper,
            threshold.quantile_upper,
            mean, std, values
        )

        lower_bound, upper_bound = self._enforce_non_negative_bounds(
            prop, lower_bound, upper_bound
        )
                
        return lower_bound, upper_bound

    @staticmethod
    def _enforce_non_negative_bounds(prop: str, lower_bound: Optional[float], 
                                    upper_bound: Optional[float]) -> Tuple[Optional[float], Optional[float]]:
        """Ensure properties that cannot be negative have non-negative bounds."""
        # Only logP can be negative
        non_negative_properties = {
            'num_atoms', 'num_rings', 'size_largest_ring', 'num_tokens',
            'tokens_atom_ratio', 'c_atom_ratio', 'longest_aliph_carbon',
            'molecular_weight', 'tpsa', 'num_rotatable_bonds',
            'num_h_donors', 'num_h_acceptors', 'fsp3', 'num_aromatic_rings',
            'num_stereocenters', 'num_heteroatoms', 'heteroatom_ratio'
        }
        
        if prop in non_negative_properties:
            if lower_bound is not None:
                lower_bound = max(0, lower_bound)
            if upper_bound is not None:
                upper_bound = max(0, upper_bound)
        
        return lower_bound, upper_bound

    def _plot_token_distribution(self, token_freq_df: pd.DataFrame, 
                                threshold_percent: float) -> None:
        """Create token frequency distribution plot with threshold visualization."""
        import matplotlib.pyplot as plt
        
        # Limit to top tokens for readability
        max_tokens = min(100, len(token_freq_df))
        plot_df = token_freq_df.head(max_tokens).copy()
        
        # Find cutoff token (last token above threshold)
        above_threshold = plot_df[plot_df['percent'] >= threshold_percent]
        if len(above_threshold) > 0:
            cutoff_idx = above_threshold.index.max()
        else:
            cutoff_idx = -1
        
        self.log(f"Creating token distribution plot with {len(plot_df)} tokens...")
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(8, len(plot_df) * 0.08)))
        
        # Left plot: Frequency bar chart (percentage)
        colors = [self.color_primary if i <= cutoff_idx else self.color_filtered 
                for i in range(len(plot_df))]
        ax1.barh(range(len(plot_df)), plot_df['percent'], color=colors, edgecolor='black', linewidth=0.5)
        ax1.set_yticks(range(len(plot_df)))
        ax1.set_yticklabels(plot_df['token'])
        ax1.set_xlabel('Frequency (%)', fontsize=12)
        ax1.set_ylabel('Token', fontsize=12)
        ax1.set_title('Token Distribution (Top 100)', fontsize=14)
        ax1.grid(axis='x', alpha=0.3)
        ax1.invert_yaxis()
        
        # Add threshold line (horizontal at cutoff token)
        if cutoff_idx >= 0:
            ax1.axhline(y=cutoff_idx + 0.5, color=self.color_threshold, linestyle='--', 
                    linewidth=1, label=f'Threshold: {threshold_percent}%')
            ax1.legend()
        
        # Right plot: Log-scale count
        ax2.plot(plot_df['count'], range(len(plot_df)), marker='o', markersize=4, 
                linewidth=1.5, color=self.color_primary)
        ax2.set_xscale('log')
        ax2.set_yticks(range(len(plot_df)))
        ax2.set_yticklabels(plot_df['token'])
        ax2.set_xlabel('Count (log scale)', fontsize=12)
        ax2.set_title('Token Counts', fontsize=14)
        ax2.grid(alpha=0.3)
        ax2.invert_yaxis()
        
        # Highlight cutoff (horizontal line)
        if cutoff_idx >= 0:
            ax2.axhline(y=cutoff_idx + 0.5, color=self.color_threshold, linestyle='--', linewidth=1)
        
        plt.tight_layout()
        
        # Save figure
        output_path = os.path.join(self.distributions_dir, 'token_distribution.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        self.log(f"Saved token distribution plot: {output_path}")

    def _compute_pca_coordinates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute PCA coordinates and add them to dataframe."""
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        
        # Select numeric properties
        properties = [p for p in self._properties_to_compute if p in df.columns]
        
        if len(properties) < 2:
            self.log("Need at least 2 properties for PCA. Skipping.", level='WARNING')
            return df
        
        # Get complete cases
        df_complete = df[properties].dropna()
        complete_idx = df_complete.index
        
        if len(df_complete) == 0:
            self.log("No complete cases for PCA. Skipping.", level='WARNING')
            return df
        
        self.log(f"Computing PCA coordinates for {len(df_complete)}/{len(df)} molecules...")
        
        # Standardize and fit PCA
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(df_complete)
        pca = PCA(n_components=min(2, len(properties)))
        X_pca = pca.fit_transform(X_scaled)
        
        # Add coordinates to dataframe
        df['distribution_pca_x'] = np.nan
        df['distribution_pca_y'] = np.nan
        df.loc[complete_idx, 'distribution_pca_x'] = X_pca[:, 0]
        if X_pca.shape[1] > 1:
            df.loc[complete_idx, 'distribution_pca_y'] = X_pca[:, 1]
        
        # Store PCA objects for plotting
        self._pca_data = {
            'pca': pca,
            'scaler': scaler,
            'properties': properties,
            'complete_idx': complete_idx
        }
        
        missing_pca = df['distribution_pca_x'].isna().sum()
        if missing_pca > 0:
            self.log(f"  {missing_pca} molecules missing PCA coordinates (incomplete data)")
        
        return df

    def _plot_pca(self, df: pd.DataFrame, sample_size: int = 5_000_000) -> None:
        """Create PCA visualization of property space."""
        import matplotlib.pyplot as plt
        
        if not hasattr(self, '_pca_data'):
            self.log("PCA data not available. Skipping plot.", level='WARNING')
            return
        
        pca_data = self._pca_data
        pca = pca_data['pca']
        complete_idx = pca_data['complete_idx']
        
        # Get PCA coordinates
        df_plot = df.loc[complete_idx, ['distribution_pca_x', 'distribution_pca_y']].dropna()
        
        if len(df_plot) == 0:
            self.log("No PCA coordinates available. Skipping plot.", level='WARNING')
            return
        
        # Sample if too large
        if len(df_plot) > sample_size:
            df_plot = df_plot.sample(sample_size, random_state=42)
        
        self.log(f"Creating PCA plot with {len(df_plot):,} molecules...")
        
        # Create plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Scatter plot
        ax1.scatter(df_plot['distribution_pca_x'], df_plot['distribution_pca_y'], 
                   alpha=0.3, s=10, c=self.color_primary)
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        ax1.set_title('PCA of Molecular Properties')
        ax1.grid(True, alpha=0.3)
        
        # Variance explained
        n_components = len(pca.explained_variance_ratio_)
        ax2.bar(range(1, n_components + 1), pca.explained_variance_ratio_ * 100, 
               alpha=0.7, color=self.color_primary, edgecolor='black')
        ax2.set_xlabel('Principal Component')
        ax2.set_ylabel('Variance Explained (%)')
        ax2.set_title('Variance Explained by Components')
        ax2.set_xticks(range(1, n_components + 1))
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Cumulative variance
        ax2_twin = ax2.twinx()
        cumsum = np.cumsum(pca.explained_variance_ratio_ * 100)
        ax2_twin.plot(range(1, n_components + 1), cumsum, 'r-o', linewidth=2, markersize=6)
        ax2_twin.set_ylabel('Cumulative Variance (%)', color='r')
        ax2_twin.tick_params(axis='y', labelcolor='r')
        
        plt.tight_layout()
        
        # Save
        output_path = os.path.join(self.distributions_dir, 'pca_properties.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        self.log(f"Saved PCA plot: {output_path}")
        self.log(f"  Total variance explained (PC1+PC2): {cumsum[1]:.1f}%")
    
    # ==================== Helper Methods ====================

    def _save_curation_results(self) -> None:
        """Save curation configuration and results to JSON file."""
        import json
        from datetime import datetime

        results = {
            'timestamp': datetime.now().isoformat(),
            'properties': {},
        }

        # Add threshold and statistics for each property
        for prop in self._properties_to_compute:
            threshold = self.thresholds.get(prop)
            if not threshold or not self._has_active_threshold(threshold):
                continue

            stats = self._statistics_cache.get(prop, {})
            prop_data = {
                'thresholds': {},
                'statistics': {k: float(v) if isinstance(v, (int, float, np.number)) else v
                              for k, v in stats.items()}
            }

            # Add threshold values
            if threshold.min_value is not None:
                prop_data['thresholds']['min_value'] = float(threshold.min_value)
            if threshold.max_value is not None:
                prop_data['thresholds']['max_value'] = float(threshold.max_value)
            if threshold.statistical_lower is not None:
                prop_data['thresholds']['statistical_lower'] = float(threshold.statistical_lower)
            if threshold.statistical_upper is not None:
                prop_data['thresholds']['statistical_upper'] = float(threshold.statistical_upper)
            if threshold.quantile_lower is not None:
                prop_data['thresholds']['quantile_lower'] = float(threshold.quantile_lower)
            if threshold.quantile_upper is not None:
                prop_data['thresholds']['quantile_upper'] = float(threshold.quantile_upper)

            results['properties'][prop] = prop_data

        # Save to file
        with open(self.curation_results_file, 'w') as f:
            json.dump(results, f, indent=2)

        self.log(f"Curation results saved to {self.curation_results_file}")

    def _log_configuration(self) -> None:
        """Log configuration details."""
        props = ' > '.join(self.properties) if isinstance(self.properties, list) else self.properties
        self.log(f"Initialized with properties: {props}", level='INFO')
        
        self.log("Thresholds:", level='INFO')
        for prop in self._properties_to_compute:
            threshold = self.thresholds[prop]
            if threshold and self._has_active_threshold(threshold):
                bounds = []
                if threshold.min_value is not None:
                    bounds.append(f"min={threshold.min_value}")
                if threshold.max_value is not None:
                    bounds.append(f"max={threshold.max_value}")
                if threshold.statistical_lower is not None:
                    bounds.append(f"stat_lower={threshold.statistical_lower}σ")
                if threshold.statistical_upper is not None:
                    bounds.append(f"stat_upper={threshold.statistical_upper}σ")
                if threshold.quantile_lower is not None:
                    bounds.append(f"quant_lower={threshold.quantile_lower}")
                if threshold.quantile_upper is not None:
                    bounds.append(f"quant_upper={threshold.quantile_upper}")
                
                self.log(f"  {prop}: {', '.join(bounds)}", level='INFO')
    
    def _report_progress(self, completed: int, total: int, start_time: float, 
                        chunk_size: int) -> None:
        """Report multiprocessing progress."""
        elapsed = time.time() - start_time
        progress_pct = 100 * completed / total
        molecules_processed = completed * chunk_size
        
        if completed > 0:
            eta = (elapsed * total / completed) - elapsed
            rate = molecules_processed / elapsed
            
            self.log(
                f"Progress: {completed}/{total} chunks ({progress_pct:.1f}%) | "
                f"Elapsed: {self._format_duration(elapsed)} | "
                f"ETA: {self._format_duration(eta)} | "
                f"Rate: {rate:.0f} mol/s"
            )
        else:
            self.log(f"Progress: {completed}/{total} chunks ({progress_pct:.1f}%)")