
# =============================================================================
# Ortholog Converter (converter.py)
# =============================================================================

"""
Converts alignment results to ortholog pair format for downstream analysis.
Handles coordinate mapping and metadata preservation.
"""

import pandas as pd
from typing import List, Dict, Any
from ..alignment.alignment_result import AlignmentResult
from ..utils import standardize_sequence_id
from ..logger import get_logger

logger = get_logger()

class OrthologConverter:
    """
    Converts alignment results to ortholog pairs format suitable for synteny analysis.
    Preserves all relevant metadata and coordinates.
    """
    
    def __init__(self, config):
        """
        Initialize ortholog converter.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.enable_debug = config.get('enable_debug_output', False)
        
        logger.info("Ortholog converter initialized")
        logger.info(f"  Debug output: {'enabled' if self.enable_debug else 'disabled'}")
    
    def convert_to_ortholog_pairs(self, alignment_results: List[AlignmentResult], 
                                sequence_pairs: List[Dict]) -> pd.DataFrame:
        """
        Convert alignment results to ortholog pairs DataFrame.
        
        Args:
            alignment_results: List of filtered AlignmentResult objects
            sequence_pairs: Original sequence pairs with coordinate information
            
        Returns:
            DataFrame with ortholog pairs
        """
        logger.info(f"Converting {len(alignment_results)} alignment results to ortholog pairs...")
        
        # Create mapping from pair index to sequence pair
        pair_map = {i: pair for i, pair in enumerate(sequence_pairs)}
        
        ortholog_pairs = []
        
        for result in alignment_results:
            pair_idx = result.pair_index
            
            if pair_idx in pair_map:
                ortholog_pair = self._create_ortholog_pair(result, pair_map[pair_idx])
                ortholog_pairs.append(ortholog_pair)
            else:
                if self.enable_debug:
                    logger.warning(f"Pair index {pair_idx} not found in sequence pairs")
        
        # Create DataFrame
        ortholog_df = pd.DataFrame(ortholog_pairs)
        
        # Log conversion statistics
        self._log_conversion_statistics(alignment_results, ortholog_df)
        
        return ortholog_df
    
    def _create_ortholog_pair(self, result: AlignmentResult, sequence_pair: Dict) -> Dict[str, Any]:
        """
        Create ortholog pair dictionary from alignment result and sequence pair.
        
        Args:
            result: AlignmentResult object
            sequence_pair: Original sequence pair data
            
        Returns:
            Dictionary with ortholog pair information
        """
        first_gene = sequence_pair['first_gene']
        second_gene = sequence_pair['second_gene']
        
        # Create base ortholog pair
        ortholog_pair = {
            # Identification
            'busco_id': result.busco_id,
            
            # First genome coordinates
            'first_chr': standardize_sequence_id(first_gene['sequence_id']),
            'first_start': first_gene['gene_start'],
            'first_end': first_gene['gene_end'],
            'first_strand': first_gene['strand'],
            'first_length': first_gene['gene_length'],
            
            # Second genome coordinates
            'second_chr': standardize_sequence_id(second_gene['sequence_id']),
            'second_start': second_gene['gene_start'],
            'second_end': second_gene['gene_end'],
            'second_strand': second_gene['strand'],
            'second_length': second_gene['gene_length'],
            
            # Alignment quality metrics
            'similarity': result.similarity,
            'confidence': result.confidence,
            'identity': result.identity,
            'coverage': result.min_coverage,
            'alignment_length': result.alignment_length,
            'matches': result.matches,
            
            # Method information
            'alignment_method': result.method,
            
            # Derived metrics
            'length_ratio': min(first_gene['gene_length'], second_gene['gene_length']) / max(first_gene['gene_length'], second_gene['gene_length']),
            'mapping_type': 'ortholog'
        }
        
        # Add method-specific metadata
        if result.method == 'minimap2' and result.mapq is not None:
            ortholog_pair['mapq'] = result.mapq
        
        if result.method == 'biopython' and result.score is not None:
            ortholog_pair['alignment_score'] = result.score
        
        # Add validation information if available
        if result.validation_method:
            ortholog_pair['validation_method'] = result.validation_method
        
        if result.alternative_confidence is not None:
            ortholog_pair['alternative_confidence'] = result.alternative_confidence
        
        # Add paralog information if available
        ortholog_pair['first_paralog_rank'] = first_gene.get('paralog_rank', 1)
        ortholog_pair['second_paralog_rank'] = second_gene.get('paralog_rank', 1)
        ortholog_pair['first_paralog_count'] = first_gene.get('paralog_count', 1)
        ortholog_pair['second_paralog_count'] = second_gene.get('paralog_count', 1)
        
        # Add debug information if enabled
        if self.enable_debug:
            ortholog_pair['pair_index'] = result.pair_index
            ortholog_pair['query_coverage'] = result.query_coverage
            ortholog_pair['target_coverage'] = result.target_coverage
        
        return ortholog_pair
    
    def _log_conversion_statistics(self, alignment_results: List[AlignmentResult], 
                                 ortholog_df: pd.DataFrame):
        """Log conversion statistics."""
        logger.info("  Conversion statistics:")
        logger.info(f"    Input alignment results: {len(alignment_results)}")
        logger.info(f"    Output ortholog pairs: {len(ortholog_df)}")
        
        if len(ortholog_df) > 0:
            # Quality statistics
            avg_similarity = ortholog_df['similarity'].mean()
            avg_confidence = ortholog_df['confidence'].mean()
            avg_identity = ortholog_df['identity'].mean()
            
            logger.info(f"    Average similarity: {avg_similarity:.3f}")
            logger.info(f"    Average confidence: {avg_confidence:.3f}")
            logger.info(f"    Average identity: {avg_identity:.3f}")
            
            # Method distribution
            method_counts = ortholog_df['alignment_method'].value_counts()
            logger.info(f"    Method distribution: {method_counts.to_dict()}")
            
            # Chromosome distribution
            unique_first_chrs = ortholog_df['first_chr'].nunique()
            unique_second_chrs = ortholog_df['second_chr'].nunique()
            logger.info(f"    Chromosomes covered: {unique_first_chrs} (first) × {unique_second_chrs} (second)")
    
    def get_conversion_statistics(self, ortholog_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Get comprehensive conversion statistics.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            
        Returns:
            Dictionary with conversion statistics
        """
        if len(ortholog_df) == 0:
            return {'total_pairs': 0}
        
        stats = {
            'total_pairs': len(ortholog_df),
            'unique_buscos': ortholog_df['busco_id'].nunique(),
            'first_chromosomes': ortholog_df['first_chr'].nunique(),
            'second_chromosomes': ortholog_df['second_chr'].nunique(),
            'quality_metrics': {
                'similarity': {
                    'mean': ortholog_df['similarity'].mean(),
                    'std': ortholog_df['similarity'].std(),
                    'min': ortholog_df['similarity'].min(),
                    'max': ortholog_df['similarity'].max()
                },
                'confidence': {
                    'mean': ortholog_df['confidence'].mean(),
                    'std': ortholog_df['confidence'].std(),
                    'min': ortholog_df['confidence'].min(),
                    'max': ortholog_df['confidence'].max()
                },
                'identity': {
                    'mean': ortholog_df['identity'].mean(),
                    'std': ortholog_df['identity'].std(),
                    'min': ortholog_df['identity'].min(),
                    'max': ortholog_df['identity'].max()
                }
            },
            'method_distribution': ortholog_df['alignment_method'].value_counts().to_dict(),
            'chromosome_pairs': len(ortholog_df.groupby(['first_chr', 'second_chr']))
        }
        
        return stats