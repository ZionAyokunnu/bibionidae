
# =============================================================================
# Quality Assessor (quality_assessor.py)
# =============================================================================

"""
Assembly quality assessment and adaptive parameter adjustment system.
Evaluates genome assembly quality and suggests optimal analysis parameters.
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, Tuple, List
from pathlib import Path
from Bio import SeqIO

from ..logger import get_logger

logger = get_logger()

class QualityAssessor:
    """
    Comprehensive assembly quality assessor that evaluates genome quality
    and provides adaptive parameter recommendations for analysis.
    """
    
    def __init__(self, config):
        """
        Initialize quality assessor.
        
        Args:
            config: Configuration object with quality assessment parameters
        """
        self.config = config
        self.enable_assessment = config.get('enable_assembly_quality_assessment', True)
        self.enable_adaptive_thresholds = config.get('enable_adaptive_thresholds', True)
        
        # Quality thresholds
        self.high_quality_busco_threshold = config.get('high_quality_busco_threshold', 0.95)
        self.medium_quality_busco_threshold = config.get('medium_quality_busco_threshold', 0.85)
        self.low_quality_busco_threshold = config.get('low_quality_busco_threshold', 0.70)
        self.high_quality_n50_threshold = config.get('high_quality_n50_threshold', 10000000)
        self.medium_quality_n50_threshold = config.get('medium_quality_n50_threshold', 1000000)
        
        logger.info("Quality assessor initialized")
        logger.info(f"  Quality assessment: {'enabled' if self.enable_assessment else 'disabled'}")
        logger.info(f"  Adaptive thresholds: {'enabled' if self.enable_adaptive_thresholds else 'disabled'}")
    
    def assess_assembly_quality(self, fasta_path: str, busco_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Comprehensive assembly quality assessment.
        
        Args:
            fasta_path: Path to genome FASTA file
            busco_df: DataFrame with BUSCO data
            
        Returns:
            Dictionary with quality assessment results
        """
        if not self.enable_assessment:
            return {
                'quality_score': 1.0, 
                'quality_class': 'medium', 
                'metrics': {},
                'adjustments': {}
            }
        
        logger.info(f"Assessing assembly quality for {fasta_path}")
        
        # Calculate assembly metrics
        assembly_metrics = self._calculate_assembly_metrics(fasta_path)
        
        # Calculate BUSCO metrics
        busco_metrics = self._calculate_busco_metrics(busco_df)
        
        # Combine all metrics
        all_metrics = {**assembly_metrics, **busco_metrics}
        
        # Calculate overall quality score
        quality_score = self._calculate_quality_score(all_metrics)
        
        # Classify quality
        quality_class = self._classify_assembly_quality(quality_score, all_metrics)
        
        # Suggest parameter adjustments
        adjustments = self._suggest_parameter_adjustments(quality_class, all_metrics)
        
        logger.info(f"  Assembly quality: {quality_class} (score: {quality_score:.3f})")
        logger.info(f"  Suggested {len(adjustments)} parameter adjustments")
        
        return {
            'quality_score': quality_score,
            'quality_class': quality_class,
            'metrics': all_metrics,
            'adjustments': adjustments
        }
    
    def _calculate_assembly_metrics(self, fasta_path: str) -> Dict[str, Any]:
        """Calculate basic assembly statistics."""
        try:
            sequences = list(SeqIO.parse(fasta_path, "fasta"))
            contig_lengths = [len(seq) for seq in sequences]
            total_length = sum(contig_lengths)
            n_contigs = len(contig_lengths)
            
            if not contig_lengths:
                return {}
            
            # Sort lengths for N-statistics
            sorted_lengths = sorted(contig_lengths, reverse=True)
            cumsum = np.cumsum(sorted_lengths)
            
            # Calculate N50
            n50_idx = np.where(cumsum >= total_length * 0.5)[0]
            n50 = sorted_lengths[n50_idx[0]] if len(n50_idx) > 0 else 0
            
            # Calculate N90
            n90_idx = np.where(cumsum >= total_length * 0.9)[0]
            n90 = sorted_lengths[n90_idx[0]] if len(n90_idx) > 0 else 0
            
            # Calculate L50 (number of contigs that make up 50% of assembly)
            l50 = len(n50_idx) + 1 if len(n50_idx) > 0 else n_contigs
            
            metrics = {
                'total_length': total_length,
                'n_contigs': n_contigs,
                'n50': n50,
                'n90': n90,
                'l50': l50,
                'max_contig': max(contig_lengths),
                'min_contig': min(contig_lengths),
                'mean_contig_length': total_length / n_contigs,
                'median_contig_length': np.median(contig_lengths),
                'contig_length_std': np.std(contig_lengths)
            }
            
            # Calculate additional quality indicators
            metrics['assembly_fragmentation'] = n_contigs / (total_length / 1000000)  # Contigs per Mb
            metrics['contig_length_cv'] = metrics['contig_length_std'] / metrics['mean_contig_length']
            
            return metrics
            
        except Exception as e:
            logger.error(f"Failed to calculate assembly metrics: {e}")
            return {}
    
    def _calculate_busco_metrics(self, busco_df: pd.DataFrame) -> Dict[str, Any]:
        """Calculate BUSCO completeness metrics."""
        if len(busco_df) == 0:
            return {}
        
        # Count BUSCO statuses
        status_counts = busco_df['status'].value_counts()
        total_buscos = len(busco_df)
        
        complete = status_counts.get('Complete', 0)
        fragmented = status_counts.get('Fragmented', 0)
        missing = status_counts.get('Missing', 0)
        duplicated = status_counts.get('Duplicated', 0)
        
        metrics = {
            'busco_total': total_buscos,
            'busco_complete': complete,
            'busco_fragmented': fragmented,
            'busco_missing': missing,
            'busco_duplicated': duplicated,
            'busco_completeness': complete / total_buscos if total_buscos > 0 else 0,
            'busco_fragmentation': fragmented / total_buscos if total_buscos > 0 else 0,
            'busco_duplication': duplicated / total_buscos if total_buscos > 0 else 0,
            'busco_missing_rate': missing / total_buscos if total_buscos > 0 else 0
        }
        
        # Calculate single-copy BUSCO rate
        single_copy = complete - duplicated
        metrics['busco_single_copy_rate'] = single_copy / total_buscos if total_buscos > 0 else 0
        
        return metrics
    
    def _calculate_quality_score(self, metrics: Dict[str, Any]) -> float:
        """Calculate comprehensive assembly quality score."""
        score_components = []
        
        # N50 component (0-1 scale)
        if 'n50' in metrics:
            if metrics['n50'] >= self.high_quality_n50_threshold:
                score_components.append(1.0)
            elif metrics['n50'] >= self.medium_quality_n50_threshold:
                score_components.append(0.7)
            else:
                score_components.append(max(0.0, metrics['n50'] / self.medium_quality_n50_# Phase 5: Analysis Modules (Synteny, Inversions, and Rearrangements)
# This phase handles the core genomic analysis for synteny, inversions, and rearrangements