# =============================================================================
# Statistical Validator (statistical_validator.py)
# =============================================================================

"""
Statistical validation system for genomic analysis results.
Provides confidence intervals, significance testing, and result validation.
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, List, Tuple, Optional
from scipy import stats
from scipy.stats import bootstrap
import warnings

from ..logger import get_logger

logger = get_logger()

class StatisticalValidator:
    """
    Statistical validation system that provides robust statistical testing
    and confidence estimation for genomic analysis results.
    """
    
    def __init__(self, config):
        """
        Initialize statistical validator.
        
        Args:
            config: Configuration object with validation parameters
        """
        self.config = config
        self.enable_validation = config.get('enable_statistical_validation', False)
        self.confidence_level = config.get('statistical_confidence_level', 0.95)
        self.bootstrap_iterations = config.get('bootstrap_iterations', 1000)
        self.significance_threshold = config.get('significance_threshold', 0.05)
        
        logger.info("Statistical validator initialized")
        logger.info(f"  Statistical validation: {'enabled' if self.enable_validation else 'disabled'}")
        logger.info(f"  Confidence level: {self.confidence_level}")
    
    def validate_synteny_results(self, synteny_df: pd.DataFrame, 
                                ortholog_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Statistical validation of synteny analysis results.
        
        Args:
            synteny_df: DataFrame with synteny blocks
            ortholog_df: DataFrame with ortholog pairs
            
        Returns:
            Dictionary with validation results
        """
        if not self.enable_validation or len(synteny_df) == 0:
            return {'validated': False, 'reason': 'validation disabled or no data'}
        
        logger.info("Performing statistical validation of synteny results...")
        
        validation_results = {}
        
        # Validate correlation significance
        correlation_validation = self._validate_correlations(synteny_df)
        validation_results['correlation_validation'] = correlation_validation
        
        # Validate block size distribution
        block_size_validation = self._validate_block_sizes(synteny_df)
        validation_results['block_size_validation'] = block_size_validation
        
        # Validate strand consistency
        strand_validation = self._validate_strand_consistency(ortholog_df)
        validation_results['strand_validation'] = strand_validation
        
        # Overall validation score
        validation_score = self._calculate_validation_score(validation_results)
        validation_results['overall_validation'] = {
            'score': validation_score,
            'validated': validation_score > 0.5,
            'confidence_level': self.confidence_level
        }
        
        logger.info(f"  Synteny validation score: {validation_score:.3f}")
        
        return validation_results
    
    def validate_inversion_results(self, inversion_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Statistical validation of inversion detection results.
        
        Args:
            inversion_df: DataFrame with detected inversions
            
        Returns:
            Dictionary with validation results
        """
        if not self.enable_validation or len(inversion_df) == 0:
            return {'validated': False, 'reason': 'validation disabled or no data'}
        
        logger.info("Performing statistical validation of inversion results...")
        
        validation_results = {}
        
        # Validate inversion size distribution
        size_validation = self._validate_inversion_sizes(inversion_df)
        validation_results['size_validation'] = size_validation
        
        # Validate confidence scores
        confidence_validation = self._validate_confidence_scores(inversion_df)
        validation_results['confidence_validation'] = confidence_validation
        
        # Validate type distribution
        type_validation = self._validate_inversion_types(inversion_df)
        validation_results['type_validation'] = type_validation
        
        # Overall validation
        validation_score = self._calculate_validation_score(validation_results)
        validation_results['overall_validation'] = {
            'score': validation_score,
            'validated': validation_score > 0.5,
            'confidence_level': self.confidence_level
        }
        
        logger.info(f"  Inversion validation score: {validation_score:.3f}")
        
        return validation_results
    
    def calculate_confidence_intervals(self, data: List[float], 
                                     statistic_func: Optional[callable] = None) -> Dict[str, float]:
        """
        Calculate bootstrap confidence intervals for a statistic.
        
        Args:
            data: List of numerical data
            statistic_func: Function to calculate statistic (default: mean)
            
        Returns:
            Dictionary with confidence interval results
        """
        if not self.enable_validation or len(data) < 10:
            return {'ci_lower': np.nan, 'ci_upper': np.nan, 'statistic': np.nan}
        
        if statistic_func is None:
            statistic_func = np.mean
        
        try:
            # Perform bootstrap
            data_array = np.array(data)
            rng = np.random.default_rng(42)  # Fixed seed for reproducibility
            
            # Use scipy's bootstrap if available
            res = bootstrap((data_array,), statistic_func, 
                          n_resamples=self.bootstrap_iterations,
                          confidence_level=self.confidence_level,
                          random_state=rng)
            
            return {
                'ci_lower': res.confidence_interval.low,
                'ci_upper': res.confidence_interval.high,
                'statistic': statistic_func(data_array)
            }
            
        except Exception as e:
            logger.warning(f"Bootstrap confidence interval calculation failed: {e}")
            
            # Fallback to simple percentile method
            bootstrap_stats = []
            for _ in range(self.bootstrap_iterations):
                sample = np.random.choice(data, size=len(data), replace=True)
                bootstrap_stats.append(statistic_func(sample))
            
            alpha = 1 - self.confidence_level
            lower_percentile = (alpha / 2) * 100
            upper_percentile = (1 - alpha / 2) * 100
            
            return {
                'ci_lower': np.percentile(bootstrap_stats, lower_percentile),
                'ci_upper': np.percentile(bootstrap_stats, upper_percentile),
                'statistic': statistic_func(data)
            }
    
    def _validate_correlations(self, synteny_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate statistical significance of position correlations."""
        if 'position_correlation' not in synteny_df.columns:
            return {'validated': False, 'reason': 'no correlation data'}
        
        correlations = synteny_df['position_correlation'].dropna()
        if len(correlations) == 0:
            return {'validated': False, 'reason': 'no valid correlations'}
        
        # Test if correlations are significantly different from zero
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t_stat, p_value = stats.ttest_1samp(correlations, 0)
        
        # Calculate confidence interval for mean correlation
        ci_results = self.calculate_confidence_intervals(correlations.tolist())
        
        # Count significant correlations
        if 'p_value' in synteny_df.columns:
            significant_correlations = (synteny_df['p_value'] < self.significance_threshold).sum()
            significance_rate = significant_correlations / len(synteny_df)
        else:
            significance_rate = np.nan
        
        return {
            'validated': p_value < self.significance_threshold,
            'mean_correlation': correlations.mean(),
            'correlation_ci': ci_results,
            't_statistic': t_stat,
            'p_value': p_value,
            'significance_rate': significance_rate,
            'n_blocks': len(correlations)
        }
    
    def _validate_block_sizes(self, synteny_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate synteny block size distribution."""
        if 'block_size' not in synteny_df.columns:
            return {'validated': False, 'reason': 'no block size data'}
        
        block_sizes = synteny_df['block_size'].dropna()
        if len(block_sizes) == 0:
            return {'validated': False, 'reason': 'no valid block sizes'}
        
        # Test for reasonable distribution
        mean_size = block_sizes.mean()
        median_size = block_sizes.median()
        
        # Calculate confidence interval for mean block size
        ci_results = self.calculate_confidence_intervals(block_sizes.tolist())
        
        # Test for log-normal distribution (common for biological data)
        log_sizes = np.log(block_sizes[block_sizes > 0])
        if len(log_sizes) > 0:
            _, normality_p = stats.normaltest(log_sizes)
            log_normal_fit = normality_p > self.significance_threshold
        else:
            log_normal_fit = False
            normality_p = np.nan
        
        return {
            'validated': mean_size > 1 and median_size >= 1,
            'mean_size': mean_size,
            'median_size': median_size,
            'size_ci': ci_results,
            'log_normal_fit': log_normal_fit,
            'normality_p_value': normality_p,
            'n_blocks': len(block_sizes)
        }
    
    def _validate_strand_consistency(self, ortholog_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate strand consistency patterns."""
        if len(ortholog_df) == 0:
            return {'validated': False, 'reason': 'no ortholog data'}
        
        # Calculate strand consistency per chromosome pair
        strand_consistencies = []
        for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
            consistency = (group['first_strand'] == group['second_strand']).mean()
            strand_consistencies.append(consistency)
        
        if not strand_consistencies:
            return {'validated': False, 'reason': 'no strand data'}
        
        # Test if strand consistency is significantly > 0.5 (random)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t_stat, p_value = stats.ttest_1samp(strand_consistencies, 0.5)
        
        # Calculate confidence interval
        ci_results = self.calculate_confidence_intervals(strand_consistencies)
        
        mean_consistency = np.mean(strand_consistencies)
        
        return {
            'validated': p_value < self.significance_threshold and mean_consistency > 0.5,
            'mean_consistency': mean_consistency,
            'consistency_ci': ci_results,
            't_statistic': t_stat,
            'p_value': p_value,
            'n_chromosome_pairs': len(strand_consistencies)
        }
    
    def _validate_inversion_sizes(self, inversion_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate inversion size distribution."""
        if 'size_genes' not in inversion_df.columns:
            return {'validated': False, 'reason': 'no size data'}
        
        sizes = inversion_df['size_genes'].dropna()
        if len(sizes) == 0:
            return {'validated': False, 'reason': 'no valid sizes'}
        
        # Calculate basic statistics
        mean_size = sizes.mean()
        median_size = sizes.median()
        
        # Calculate confidence interval
        ci_results = self.calculate_confidence_intervals(sizes.tolist())
        
        # Test for reasonable size distribution
        reasonable_sizes = (sizes >= 1) & (sizes <= 100)  # Reasonable range
        reasonable_rate = reasonable_sizes.mean()
        
        return {
            'validated': reasonable_rate > 0.8 and mean_size >= 1,
            'mean_size': mean_size,
            'median_size': median_size,
            'size_ci': ci_results,
            'reasonable_rate': reasonable_rate,
            'n_inversions': len(sizes)
        }
    
    def _validate_confidence_scores(self, result_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate confidence score distribution."""
        if 'confidence' not in result_df.columns:
            return {'validated': False, 'reason': 'no confidence data'}
        
        confidences = result_df['confidence'].dropna()
        if len(confidences) == 0:
            return {'validated': False, 'reason': 'no valid confidence scores'}
        
        # Check if confidences are in valid range [0,1]
        valid_range = (confidences >= 0) & (confidences <= 1)
        valid_rate = valid_range.mean()
        
        # Calculate statistics
        mean_confidence = confidences.mean()
        
        # Calculate confidence interval
        ci_results = self.calculate_confidence_intervals(confidences.tolist())
        
        # Test if confidence scores are reasonable (not all the same)
        confidence_variance = confidences.var()
        
        return {
            'validated': valid_rate > 0.95 and confidence_variance > 0.01,
            'mean_confidence': mean_confidence,
            'confidence_ci': ci_results,
            'valid_range_rate': valid_rate,
            'variance': confidence_variance,
            'n_scores': len(confidences)
        }
    
    def _validate_inversion_types(self, inversion_df: pd.DataFrame) -> Dict[str, Any]:
        """Validate inversion type distribution."""
        if 'inversion_type' not in inversion_df.columns:
            return {'validated': False, 'reason': 'no type data'}
        
        type_counts = inversion_df['inversion_type'].value_counts()
        
        # Check for reasonable type diversity
        n_types = len(type_counts)
        entropy = stats.entropy(type_counts.values)
        max_entropy = np.log(n_types) if n_types > 1 else 0
        
        return {
            'validated': n_types > 0 and entropy >= 0,
            'n_types': n_types,
            'entropy': entropy,
            'max_entropy': max_entropy,
            'normalized_entropy': entropy / max_entropy if max_entropy > 0 else 0,
            'type_distribution': type_counts.to_dict()
        }
    
    def _calculate_validation_score(self, validation_results: Dict[str, Any]) -> float:
        """Calculate overall validation score from individual validation results."""
        scores = []
        
        for key, result in validation_results.items():
            if isinstance(result, dict) and 'validated' in result:
                scores.append(1.0 if result['validated'] else 0.0)
        
        return np.mean(scores) if scores else 0.0