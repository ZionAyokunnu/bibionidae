
# =============================================================================
# Inversion Detector (inversion_detector.py)
# =============================================================================

"""
Advanced inversion detection system with confidence scoring.
Identifies different types of inversions and assesses their reliability.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional
from collections import defaultdict

from ..logger import get_logger

logger = get_logger()

class InversionDetector:
    """
    Advanced inversion detection system that identifies various types of inversions
    and provides confidence scoring for detected events.
    """
    
    def __init__(self, config):
        """
        Initialize inversion detector.
        
        Args:
            config: Configuration object with inversion parameters
        """
        self.config = config
        self.min_inversion_size = config.get('base_min_inversion_size', 2)
        self.enable_single_gene = config.get('enable_single_gene_inversions', True)
        self.enable_micro_inversions = config.get('enable_micro_inversions', True)
        self.enable_confidence = config.get('enable_inversion_confidence', True)
        self.confidence_threshold = config.get('inversion_confidence_threshold', 0.7)
        
        if self.enable_single_gene:
            self.min_inversion_size = config.get('micro_inversion_size', 1)
        
        logger.info("Inversion detector initialized")
        logger.info(f"  Min inversion size: {self.min_inversion_size}")
        logger.info(f"  Single gene inversions: {'enabled' if self.enable_single_gene else 'disabled'}")
        logger.info(f"  Confidence scoring: {'enabled' if self.enable_confidence else 'disabled'}")
    
    def detect_inversions(self, synteny_df: pd.DataFrame, 
                         ortholog_df: pd.DataFrame) -> pd.DataFrame:
        """
        Detect inversions from synteny blocks and ortholog data.
        
        Args:
            synteny_df: DataFrame with synteny blocks
            ortholog_df: DataFrame with ortholog pairs
            
        Returns:
            DataFrame with detected inversions
        """
        if len(ortholog_df) == 0:
            logger.warning("No ortholog pairs found for inversion analysis")
            return pd.DataFrame()
        
        logger.info(f"Detecting inversions from {len(ortholog_df)} ortholog pairs...")
        
        inversions = []
        
        # Method 1: Detect inversions within synteny blocks
        if len(synteny_df) > 0 and self.enable_micro_inversions:
            block_inversions = self._detect_block_inversions(synteny_df, ortholog_df)
            inversions.extend(block_inversions)
        
        # Method 2: Detect single gene inversions
        if self.enable_single_gene:
            single_inversions = self._detect_single_gene_inversions(ortholog_df)
            inversions.extend(single_inversions)
        
        # Method 3: Detect strand inconsistency patterns
        pattern_inversions = self._detect_strand_patterns(ortholog_df)
        inversions.extend(pattern_inversions)
        
        inversion_df = pd.DataFrame(inversions)
        
        # Filter by confidence threshold if enabled
        if self.enable_confidence and len(inversion_df) > 0:
            inversion_df = inversion_df[inversion_df['confidence'] >= self.confidence_threshold]
        
        # Log detection results
        self._log_inversion_results(inversion_df, ortholog_df)
        
        return inversion_df
    
    def _detect_block_inversions(self, synteny_df: pd.DataFrame, 
                               ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect inversions within synteny blocks."""
        inversions = []
        
        for _, block in synteny_df.iterrows():
            # Get genes in this synteny block
            block_genes = ortholog_df[
                (ortholog_df['first_chr'] == block['first_chr']) & 
                (ortholog_df['second_chr'] == block['second_chr'])
            ].sort_values('first_start')
            
            if len(block_genes) >= self.min_inversion_size:
                # Detect inversion regions within this block
                block_inversions = self._identify_inversion_regions(block_genes, block)
                inversions.extend(block_inversions)
        
        return inversions
    
    def _identify_inversion_regions(self, genes: pd.DataFrame, 
                                  synteny_block: pd.Series) -> List[Dict]:
        """Identify inversion regions within a set of genes."""
        inversions = []
        genes_sorted = genes.sort_values('first_start')
        
        current_region = []
        current_inversion_type = None
        
        for _, gene in genes_sorted.iterrows():
            strand_flipped = gene['first_strand'] != gene['second_strand']
            
            if strand_flipped:
                # This gene is inverted
                current_region.append(gene)
                if current_inversion_type is None:
                    current_inversion_type = 'strand_inversion'
            else:
                # This gene is not inverted - end current region if it exists
                if len(current_region) >= self.min_inversion_size:
                    inversion = self._create_inversion_record(
                        current_region, current_inversion_type, synteny_block
                    )
                    inversions.append(inversion)
                
                current_region = []
                current_inversion_type = None
        
        # Handle final region
        if len(current_region) >= self.min_inversion_size:
            inversion = self._create_inversion_record(
                current_region, current_inversion_type, synteny_block
            )
            inversions.append(inversion)
        
        return inversions
    
    def _create_inversion_record(self, genes: List[pd.Series], inversion_type: str,
                               synteny_block: pd.Series) -> Dict:
        """Create an inversion record from a list of genes."""
        confidence = self._calculate_inversion_confidence(genes) if self.enable_confidence else 1.0
        
        first_gene = genes[0]
        last_gene = genes[-1]
        
        return {
            'first_chr': synteny_block['first_chr'],
            'second_chr': synteny_block['second_chr'],
            'start_gene': first_gene['busco_id'],
            'end_gene': last_gene['busco_id'],
            'first_start': first_gene['first_start'],
            'first_end': last_gene['first_start'],
            'second_start': min(gene['second_start'] for gene in genes),
            'second_end': max(gene['second_start'] for gene in genes),
            'size_genes': len(genes),
            'size_bp_first': last_gene['first_start'] - first_gene['first_start'],
            'size_bp_second': max(gene['second_start'] for gene in genes) - min(gene['second_start'] for gene in genes),
            'inversion_type': inversion_type,
            'strand_pattern': 'consistent_flip',
            'confidence': confidence,
            'synteny_context': synteny_block.get('synteny_type', 'unknown'),
            'average_similarity': np.mean([gene['similarity'] for gene in genes if 'similarity' in gene]),
            'detection_method': 'block_analysis'
        }
    
    def _detect_single_gene_inversions(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect single gene inversions."""
        inversions = []
        
        for _, gene in ortholog_df.iterrows():
            if gene['first_strand'] != gene['second_strand']:
                confidence = gene.get('confidence', 0.5)
                
                # Only include high-confidence single gene inversions
                if confidence >= self.confidence_threshold:
                    inversion = {
                        'first_chr': gene['first_chr'],
                        'second_chr': gene['second_chr'],
                        'start_gene': gene['busco_id'],
                        'end_gene': gene['busco_id'],
                        'first_start': gene['first_start'],
                        'first_end': gene['first_end'],
                        'second_start': gene['second_start'],
                        'second_end': gene['second_end'],
                        'size_genes': 1,
                        'size_bp_first': gene['first_end'] - gene['first_start'],
                        'size_bp_second': gene['second_end'] - gene['second_start'],
                        'inversion_type': 'single_gene_inversion',
                        'strand_pattern': 'single_flip',
                        'confidence': confidence,
                        'synteny_context': 'single_gene',
                        'average_similarity': gene.get('similarity', None),
                        'detection_method': 'single_gene'
                    }
                    inversions.append(inversion)
        
        return inversions
    
    def _detect_strand_patterns(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect inversion patterns based on strand inconsistencies."""
        inversions = []
        
        # Group by chromosome pairs
        for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
            if len(group) >= 3:  # Need enough genes to detect patterns
                pattern_inversions = self._analyze_strand_patterns(group, first_chr, second_chr)
                inversions.extend(pattern_inversions)
        
        return inversions
    
    def _analyze_strand_patterns(self, group: pd.DataFrame, first_chr: str, 
                               second_chr: str) -> List[Dict]:
        """Analyze strand patterns to detect complex inversions."""
        inversions = []
        
        # Sort by position
        group_sorted = group.sort_values('first_start')
        
        # Calculate sliding window of strand consistency
        window_size = 3
        for i in range(len(group_sorted) - window_size + 1):
            window = group_sorted.iloc[i:i + window_size]
            
            # Check if this window has consistent strand flipping
            strand_flips = (window['first_strand'] != window['second_strand']).sum()
            
            if strand_flips >= window_size - 1:  # Allow one gene to be non-flipped
                confidence = self._calculate_window_confidence(window)
                
                if confidence >= self.confidence_threshold:
                    inversion = {
                        'first_chr': first_chr,
                        'second_chr': second_chr,
                        'start_gene': window.iloc[0]['busco_id'],
                        'end_gene': window.iloc[-1]['busco_id'],
                        'first_start': window['first_start'].min(),
                        'first_end': window['first_start'].max(),
                        'second_start': window['second_start'].min(),
                        'second_end': window['second_start'].max(),
                        'size_genes': len(window),
                        'size_bp_first': window['first_start'].max() - window['first_start'].min(),
                        'size_bp_second': window['second_start'].max() - window['second_start'].min(),
                        'inversion_type': 'pattern_inversion',
                        'strand_pattern': f'window_{window_size}',
                        'confidence': confidence,
                        'synteny_context': 'pattern_based',
                        'average_similarity': window['similarity'].mean() if 'similarity' in window.columns else None,
                        'detection_method': 'pattern_analysis'
                    }
                    inversions.append(inversion)
        
        return inversions
    
    def _calculate_inversion_confidence(self, genes: List[pd.Series]) -> float:
        """Calculate confidence score for an inversion."""
        factors = []
        
        # Size factor (larger inversions are more confident)
        size_factor = min(1.0, len(genes) / 5.0)
        factors.append(size_factor)
        
        # Quality factor (average similarity of genes in inversion)
        similarities = [gene.get('similarity', 0.5) for gene in genes]
        if similarities:
            similarity_factor = np.mean(similarities)
            factors.append(similarity_factor)
        
        # Consistency factor (how consistently genes are inverted)
        strand_flips = sum(1 for gene in genes if gene.get('first_strand') != gene.get('second_strand'))
        consistency_factor = strand_flips / len(genes) if genes else 0
        factors.append(consistency_factor)
        
        # Confidence factor (average confidence of gene alignments)
        confidences = [gene.get('confidence', 0.5) for gene in genes]
        if confidences:
            conf_factor = np.mean(confidences)
            factors.append(conf_factor)
        
        return np.mean(factors) if factors else 0.5
    
    def _calculate_window_confidence(self, window: pd.DataFrame) -> float:
        """Calculate confidence for a sliding window inversion."""
        factors = []
        
        # Strand consistency factor
        strand_flips = (window['first_strand'] != window['second_strand']).sum()
        consistency_factor = strand_flips / len(window)
        factors.append(consistency_factor)
        
        # Average similarity factor
        if 'similarity' in window.columns:
            similarity_factor = window['similarity'].mean()
            factors.append(similarity_factor)
        
        # Average confidence factor
        if 'confidence' in window.columns:
            conf_factor = window['confidence'].mean()
            factors.append(conf_factor)
        
        return np.mean(factors) if factors else 0.5
    
    def _log_inversion_results(self, inversion_df: pd.DataFrame, ortholog_df: pd.DataFrame):
        """Log comprehensive inversion detection results."""
        logger.info("  Inversion detection results:")
        logger.info(f"    Total inversions detected: {len(inversion_df)}")
        
        if len(inversion_df) > 0:
            # Size distribution
            avg_size = inversion_df['size_genes'].mean()
            max_size = inversion_df['size_genes'].max()
            logger.info(f"    Average inversion size: {avg_size:.1f} genes")
            logger.info(f"    Largest inversion: {max_size} genes")
            
            # Type distribution
            if 'inversion_type' in inversion_df.columns:
                type_counts = inversion_df['inversion_type'].value_counts()
                logger.info(f"    Inversion types: {type_counts.to_dict()}")
            
            # Method distribution
            if 'detection_method' in inversion_df.columns:
                method_counts = inversion_df['detection_method'].value_counts()
                logger.info(f"    Detection methods: {method_counts.to_dict()}")
            
            # Quality metrics
            if 'confidence' in inversion_df.columns:
                avg_confidence = inversion_df['confidence'].mean()
                high_conf = (inversion_df['confidence'] >= 0.8).sum()
                logger.info(f"    Average confidence: {avg_confidence:.3f}")
                logger.info(f"    High confidence inversions: {high_conf}")