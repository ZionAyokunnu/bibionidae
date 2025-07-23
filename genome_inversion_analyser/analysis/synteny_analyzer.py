# =============================================================================
# Synteny Analyzer (synteny_analyzer.py)
# =============================================================================

"""
Comprehensive synteny block detection and analysis system.
Identifies conserved gene order and calculates synteny confidence scores.
"""

import pandas as pd
import numpy as np
from typing import Dict, Any, Tuple, List
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict

from ..logger import get_logger

logger = get_logger()

class SyntenyAnalyzer:
    """
    Comprehensive synteny analyzer that detects conserved gene order blocks
    and provides detailed synteny analysis with confidence scoring.
    """
    
    def __init__(self, config):
        """
        Initialize synteny analyzer.
        
        Args:
            config: Configuration object with synteny parameters
        """
        self.config = config
        
        # Core synteny parameters
        self.min_genes_per_block = config.get('base_min_genes_per_chromosome', 3)
        self.min_synteny_block_size = config.get('base_min_synteny_block_size', 3)
        self.correlation_threshold = config.get('base_synteny_correlation_threshold', 0.5)
        self.strand_consistency_threshold = config.get('strand_consistency_threshold', 0.6)
        self.max_gap_in_synteny = config.get('base_max_gap_in_synteny', 1000000)
        
        # Enhanced features
        self.enable_small_synteny_blocks = config.get('enable_small_synteny_blocks', True)
        self.enable_synteny_confidence_scoring = config.get('enable_synteny_confidence_scoring', True)
        self.strict_correlation_threshold = config.get('strict_correlation_threshold', 0.8)
        self.relaxed_correlation_threshold = config.get('relaxed_correlation_threshold', 0.3)
        
        # Adaptive thresholds
        if self.enable_small_synteny_blocks:
            self.min_synteny_block_size = config.get('micro_synteny_block_size', 1)
        
        logger.info("Synteny analyzer initialized")
        logger.info(f"  Min genes per block: {self.min_genes_per_block}")
        logger.info(f"  Min synteny block size: {self.min_synteny_block_size}")
        logger.info(f"  Correlation threshold: {self.correlation_threshold}")
        logger.info(f"  Strand consistency threshold: {self.strand_consistency_threshold}")
    
    def analyze_synteny_blocks(self, ortholog_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Analyze synteny blocks from ortholog data.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            
        Returns:
            Tuple of (synteny_blocks_df, chromosome_mappings_df)
        """
        if len(ortholog_df) == 0:
            logger.warning("No ortholog pairs found for synteny analysis")
            return pd.DataFrame(), pd.DataFrame()
        
        logger.info(f"Analyzing synteny blocks from {len(ortholog_df)} ortholog pairs...")
        
        synteny_blocks = []
        chromosome_mappings = []
        
        # Group by chromosome pairs
        for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
            if len(group) >= self.min_genes_per_block:
                # Analyze this chromosome pair
                blocks, mapping = self._analyze_chromosome_pair(group, first_chr, second_chr)
                synteny_blocks.extend(blocks)
                if mapping:
                    chromosome_mappings.append(mapping)
        
        synteny_df = pd.DataFrame(synteny_blocks)
        mapping_df = pd.DataFrame(chromosome_mappings)
        
        # Log results
        self._log_synteny_results(synteny_df, mapping_df, ortholog_df)
        
        return synteny_df, mapping_df
    
    def _analyze_chromosome_pair(self, group: pd.DataFrame, first_chr: str, 
                                second_chr: str) -> Tuple[List[Dict], Dict]:
        """Analyze synteny for a specific chromosome pair."""
        # Sort by position in first genome
        group_sorted = group.sort_values('first_start').reset_index(drop=True)
        
        # Calculate overall correlation
        if len(group_sorted) > 1:
            try:
                correlation, p_value = pearsonr(group_sorted['first_start'], group_sorted['second_start'])
                spearman_corr, spearman_p = spearmanr(group_sorted['first_start'], group_sorted['second_start'])
            except:
                correlation, p_value = 0.0, 1.0
                spearman_corr, spearman_p = 0.0, 1.0
        else:
            correlation, p_value = 1.0, 0.0
            spearman_corr, spearman_p = 1.0, 0.0
        
        # Calculate strand consistency
        strand_consistency = self._calculate_strand_consistency(group_sorted)
        
        # Create chromosome mapping record
        mapping = {
            'first_chr': first_chr,
            'second_chr': second_chr,
            'gene_count': len(group),
            'position_correlation': correlation,
            'spearman_correlation': spearman_corr,
            'strand_consistency': strand_consistency,
            'p_value': p_value,
            'spearman_p_value': spearman_p,
            'total_length_first': group_sorted['first_end'].max() - group_sorted['first_start'].min(),
            'total_length_second': group_sorted['second_end'].max() - group_sorted['second_start'].min()
        }
        
        # Detect synteny blocks within this chromosome pair
        blocks = self._detect_synteny_blocks_in_pair(group_sorted, first_chr, second_chr, correlation, strand_consistency)
        
        return blocks, mapping
    
    def _detect_synteny_blocks_in_pair(self, group_sorted: pd.DataFrame, first_chr: str, 
                                     second_chr: str, overall_correlation: float, 
                                     overall_strand_consistency: float) -> List[Dict]:
        """Detect individual synteny blocks within a chromosome pair."""
        blocks = []
        
        # Method 1: Sliding window approach for local synteny detection
        window_blocks = self._sliding_window_synteny_detection(group_sorted, first_chr, second_chr)
        blocks.extend(window_blocks)
        
        # Method 2: Gap-based block detection
        gap_blocks = self._gap_based_block_detection(group_sorted, first_chr, second_chr)
        blocks.extend(gap_blocks)
        
        # Method 3: Overall block if chromosome pair shows good synteny
        if (abs(overall_correlation) >= self.correlation_threshold and 
            overall_strand_consistency >= self.strand_consistency_threshold and 
            len(group_sorted) >= self.min_synteny_block_size):
            
            confidence = self._calculate_synteny_confidence(group_sorted) if self.enable_synteny_confidence_scoring else 1.0
            synteny_type = self._classify_synteny_type(overall_correlation, overall_strand_consistency)
            
            overall_block = {
                'first_chr': first_chr,
                'second_chr': second_chr,
                'block_id': f"{first_chr}_{second_chr}_overall",
                'start_gene': group_sorted.iloc[0]['busco_id'],
                'end_gene': group_sorted.iloc[-1]['busco_id'],
                'first_start': group_sorted['first_start'].min(),
                'first_end': group_sorted['first_end'].max(),
                'second_start': group_sorted['second_start'].min(), 
                'second_end': group_sorted['second_end'].max(),
                'block_size': len(group_sorted),
                'position_correlation': overall_correlation,
                'strand_consistency': overall_strand_consistency,
                'synteny_type': synteny_type,
                'confidence': confidence,
                'detection_method': 'overall_chromosome',
                'block_length_first': group_sorted['first_end'].max() - group_sorted['first_start'].min(),
                'block_length_second': group_sorted['second_end'].max() - group_sorted['second_start'].min(),
                'gene_density': len(group_sorted) / max(1, (group_sorted['first_end'].max() - group_sorted['first_start'].min()) / 1000000)
            }
            blocks.append(overall_block)
        
        # Remove duplicate/overlapping blocks and return best ones
        return self._filter_best_blocks(blocks)
    
    def _sliding_window_synteny_detection(self, group_sorted: pd.DataFrame, 
                                        first_chr: str, second_chr: str) -> List[Dict]:
        """Detect synteny blocks using sliding window approach."""
        blocks = []
        window_size = max(self.min_synteny_block_size, 3)
        
        for i in range(len(group_sorted) - window_size + 1):
            window = group_sorted.iloc[i:i + window_size]
            
            # Calculate local correlation and strand consistency
            if len(window) >= 3:
                try:
                    local_corr, _ = pearsonr(window['first_start'], window['second_start'])
                except:
                    local_corr = 0.0
                
                local_strand_consistency = self._calculate_strand_consistency(window)
                
                # Check if this window shows synteny
                if (abs(local_corr) >= self.correlation_threshold and 
                    local_strand_consistency >= self.strand_consistency_threshold):
                    
                    confidence = self._calculate_synteny_confidence(window) if self.enable_synteny_confidence_scoring else 1.0
                    synteny_type = self._classify_synteny_type(local_corr, local_strand_consistency)
                    
                    block = {
                        'first_chr': first_chr,
                        'second_chr': second_chr,
                        'block_id': f"{first_chr}_{second_chr}_window_{i}",
                        'start_gene': window.iloc[0]['busco_id'],
                        'end_gene': window.iloc[-1]['busco_id'],
                        'first_start': window['first_start'].min(),
                        'first_end': window['first_end'].max(),
                        'second_start': window['second_start'].min(),
                        'second_end': window['second_start'].max(),
                        'block_size': len(window),
                        'position_correlation': local_corr,
                        'strand_consistency': local_strand_consistency,
                        'synteny_type': synteny_type,
                        'confidence': confidence,
                        'detection_method': 'sliding_window',
                        'window_start': i,
                        'block_length_first': window['first_end'].max() - window['first_start'].min(),
                        'block_length_second': window['second_end'].max() - window['second_start'].min(),
                        'gene_density': len(window) / max(1, (window['first_end'].max() - window['first_start'].min()) / 1000000)
                    }
                    blocks.append(block)
        
        return blocks
    
    def _gap_based_block_detection(self, group_sorted: pd.DataFrame, 
                                 first_chr: str, second_chr: str) -> List[Dict]:
        """Detect synteny blocks based on genomic gaps."""
        blocks = []
        
        # Calculate gaps between consecutive genes
        gaps_first = group_sorted['first_start'].diff().fillna(0)
        gaps_second = group_sorted['second_start'].diff().fillna(0)
        
        # Identify potential break points based on large gaps
        large_gaps = (gaps_first > self.max_gap_in_synteny) | (gaps_second > self.max_gap_in_synteny)
        
        # Split into blocks at break points
        block_starts = [0] + (large_gaps[large_gaps].index.tolist())
        block_ends = block_starts[1:] + [len(group_sorted)]
        
        for start_idx, end_idx in zip(block_starts, block_ends):
            block_genes = group_sorted.iloc[start_idx:end_idx]
            
            if len(block_genes) >= self.min_synteny_block_size:
                # Calculate block statistics
                try:
                    block_corr, _ = pearsonr(block_genes['first_start'], block_genes['second_start'])
                except:
                    block_corr = 0.0
                
                block_strand_consistency = self._calculate_strand_consistency(block_genes)
                
                if (abs(block_corr) >= self.relaxed_correlation_threshold and 
                    block_strand_consistency >= self.strand_consistency_threshold * 0.8):  # Slightly relaxed
                    
                    confidence = self._calculate_synteny_confidence(block_genes) if self.enable_synteny_confidence_scoring else 1.0
                    synteny_type = self._classify_synteny_type(block_corr, block_strand_consistency)
                    
                    block = {
                        'first_chr': first_chr,
                        'second_chr': second_chr,
                        'block_id': f"{first_chr}_{second_chr}_gap_{start_idx}_{end_idx}",
                        'start_gene': block_genes.iloc[0]['busco_id'],
                        'end_gene': block_genes.iloc[-1]['busco_id'],
                        'first_start': block_genes['first_start'].min(),
                        'first_end': block_genes['first_end'].max(),
                        'second_start': block_genes['second_start'].min(),
                        'second_end': block_genes['second_start'].max(),
                        'block_size': len(block_genes),
                        'position_correlation': block_corr,
                        'strand_consistency': block_strand_consistency,
                        'synteny_type': synteny_type,
                        'confidence': confidence,
                        'detection_method': 'gap_based',
                        'block_length_first': block_genes['first_end'].max() - block_genes['first_start'].min(),
                        'block_length_second': block_genes['second_end'].max() - block_genes['second_start'].min(),
                        'gene_density': len(block_genes) / max(1, (block_genes['first_end'].max() - block_genes['first_start'].min()) / 1000000)
                    }
                    blocks.append(block)
        
        return blocks
    
    def _filter_best_blocks(self, blocks: List[Dict]) -> List[Dict]:
        """Filter overlapping blocks and return the best ones."""
        if not blocks:
            return blocks
        
        # Sort blocks by confidence score (descending) and size (descending)
        blocks.sort(key=lambda x: (x.get('confidence', 0), x['block_size']), reverse=True)
        
        # Remove lower-quality duplicates and heavily overlapping blocks
        filtered_blocks = []
        for block in blocks:
            # Check for significant overlap with existing blocks
            overlap_found = False
            for existing_block in filtered_blocks:
                if self._blocks_significantly_overlap(block, existing_block):
                    overlap_found = True
                    break
            
            if not overlap_found:
                filtered_blocks.append(block)
        
        return filtered_blocks
    
    def _blocks_significantly_overlap(self, block1: Dict, block2: Dict, threshold: float = 0.7) -> bool:
        """Check if two blocks overlap significantly."""
        if block1['first_chr'] != block2['first_chr'] or block1['second_chr'] != block2['second_chr']:
            return False
        
        # Calculate overlap in first genome
        overlap_start_first = max(block1['first_start'], block2['first_start'])
        overlap_end_first = min(block1['first_end'], block2['first_end'])
        overlap_first = max(0, overlap_end_first - overlap_start_first)
        
        # Calculate overlap in second genome  
        overlap_start_second = max(block1['second_start'], block2['second_start'])
        overlap_end_second = min(block1['second_end'], block2['second_end'])
        overlap_second = max(0, overlap_end_second - overlap_start_second)
        
        # Calculate relative overlap
        length1_first = block1['first_end'] - block1['first_start']
        length1_second = block1['second_end'] - block1['second_start']
        length2_first = block2['first_end'] - block2['first_start']
        length2_second = block2['second_end'] - block2['second_start']
        
        if length1_first > 0 and length1_second > 0 and length2_first > 0 and length2_second > 0:
            overlap_ratio_first = min(overlap_first / length1_first, overlap_first / length2_first)
            overlap_ratio_second = min(overlap_second / length1_second, overlap_second / length2_second)
            
            return max(overlap_ratio_first, overlap_ratio_second) > threshold
        
        return False
    
    def _calculate_strand_consistency(self, genes: pd.DataFrame) -> float:
        """Calculate strand consistency for a group of genes."""
        if len(genes) == 0:
            return 0.0
        
        return (genes['first_strand'] == genes['second_strand']).mean()
    
    def _classify_synteny_type(self, correlation: float, strand_consistency: float) -> str:
        """Classify synteny type based on correlation and strand consistency."""
        if correlation > self.correlation_threshold:
            if strand_consistency > self.strand_consistency_threshold:
                return 'colinear'
            else:
                return 'colinear_inverted'
        elif correlation < -self.correlation_threshold:
            if strand_consistency < (1 - self.strand_consistency_threshold):
                return 'inverted'
            else:
                return 'inverted_mixed'
        else:
            if strand_consistency > self.strand_consistency_threshold:
                return 'rearranged_conserved_strand'
            else:
                return 'rearranged'
    
    def _calculate_synteny_confidence(self, genes: pd.DataFrame) -> float:
        """Calculate confidence score for synteny block."""
        factors = []
        
        # Size factor (larger blocks more confident)
        size_factor = min(1.0, len(genes) / 10.0)
        factors.append(size_factor)
        
        # Quality factor (average similarity of genes in block)
        if 'similarity' in genes.columns:
            similarity_factor = genes['similarity'].mean()
            factors.append(similarity_factor)
        
        # Alignment confidence factor
        if 'confidence' in genes.columns:
            conf_factor = genes['confidence'].mean()
            factors.append(conf_factor)
        
        # Position consistency factor
        if len(genes) > 2:
            try:
                # Check how well ordered the genes are
                first_ranks = genes['first_start'].rank()
                second_ranks = genes['second_start'].rank()
                rank_corr, _ = pearsonr(first_ranks, second_ranks)
                position_factor = abs(rank_corr)
                factors.append(position_factor)
            except:
                pass
        
        # Length consistency factor
        if 'first_length' in genes.columns and 'second_length' in genes.columns:
            length_ratios = genes['first_length'] / genes['second_length'].replace(0, 1)
            length_consistency = 1.0 - np.std(length_ratios) / max(np.mean(length_ratios), 1)
            length_factor = max(0.0, min(1.0, length_consistency))
            factors.append(length_factor)
        
        return np.mean(factors) if factors else 0.5
    
    def _log_synteny_results(self, synteny_df: pd.DataFrame, mapping_df: pd.DataFrame, 
                           ortholog_df: pd.DataFrame):
        """Log comprehensive synteny analysis results."""
        logger.info("  Synteny analysis results:")
        logger.info(f"    Chromosome pairs analyzed: {len(mapping_df)}")
        logger.info(f"    Synteny blocks detected: {len(synteny_df)}")
        
        if len(synteny_df) > 0:
            # Size distribution
            avg_size = synteny_df['block_size'].mean()
            max_size = synteny_df['block_size'].max()
            total_genes_in_blocks = synteny_df['block_size'].sum()
            
            logger.info(f"    Average block size: {avg_size:.1f} genes")
            logger.info(f"    Largest block: {max_size} genes")
            logger.info(f"    Genes in synteny blocks: {total_genes_in_blocks}/{len(ortholog_df)} ({total_genes_in_blocks/len(ortholog_df)*100:.1f}%)")
            
            # Type distribution
            if 'synteny_type' in synteny_df.columns:
                type_counts = synteny_df['synteny_type'].value_counts()
                logger.info(f"    Synteny types: {type_counts.to_dict()}")
            
            # Method distribution
            if 'detection_method' in synteny_df.columns:
                method_counts = synteny_df['detection_method'].value_counts()
                logger.info(f"    Detection methods: {method_counts.to_dict()}")
            
            # Quality metrics
            if 'confidence' in synteny_df.columns:
                avg_confidence = synteny_df['confidence'].mean()
                high_conf = (synteny_df['confidence'] >= 0.8).sum()
                logger.info(f"    Average confidence: {avg_confidence:.3f}")
                logger.info(f"    High confidence blocks: {high_conf}")
        
        if len(mapping_df) > 0:
            # Correlation statistics
            if 'position_correlation' in mapping_df.columns:
                avg_corr = mapping_df['position_correlation'].mean()
                strong_corr = (mapping_df['position_correlation'].abs() >= self.correlation_threshold).sum()
                logger.info(f"    Average position correlation: {avg_corr:.3f}")
                logger.info(f"    Strong correlations: {strong_corr}/{len(mapping_df)}")
            
            # Strand consistency
            if 'strand_consistency' in mapping_df.columns:
                avg_strand = mapping_df['strand_consistency'].mean()
                consistent_pairs = (mapping_df['strand_consistency'] >= self.strand_consistency_threshold).sum()
                logger.info(f"    Average strand consistency: {avg_strand:.3f}")
                logger.info(f"    Consistent strand pairs: {consistent_pairs}/{len(mapping_df)}")
    
    def get_synteny_statistics(self, synteny_df: pd.DataFrame, mapping_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Get comprehensive synteny analysis statistics.
        
        Args:
            synteny_df: DataFrame with synteny blocks
            mapping_df: DataFrame with chromosome mappings
            
        Returns:
            Dictionary with synteny statistics
        """
        stats = {
            'total_synteny_blocks': len(synteny_df),
            'total_chromosome_pairs': len(mapping_df)
        }
        
        if len(synteny_df) > 0:
            # Block size statistics
            stats['block_size_stats'] = {
                'mean': synteny_df['block_size'].mean(),
                'std': synteny_df['block_size'].std(),
                'min': synteny_df['block_size'].min(),
                'max': synteny_df['block_size'].max(),
                'median': synteny_df['block_size'].median(),
                'total_genes': synteny_df['block_size'].sum()
            }
            
            # Type distribution
            if 'synteny_type' in synteny_df.columns:
                stats['type_distribution'] = synteny_df['synteny_type'].value_counts().to_dict()
            
            # Method distribution
            if 'detection_method' in synteny_df.columns:
                stats['method_distribution'] = synteny_df['detection_method'].value_counts().to_dict()
            
            # Confidence statistics
            if 'confidence' in synteny_df.columns:
                stats['confidence_stats'] = {
                    'mean': synteny_df['confidence'].mean(),
                    'std': synteny_df['confidence'].std(),
                    'high_confidence_count': (synteny_df['confidence'] >= 0.8).sum(),
                    'medium_confidence_count': ((synteny_df['confidence'] >= 0.5) & (synteny_df['confidence'] < 0.8)).sum(),
                    'low_confidence_count': (synteny_df['confidence'] < 0.5).sum()
                }
            
            # Length statistics
            if 'block_length_first' in synteny_df.columns:
                stats['block_length_stats'] = {
                    'mean_length_first': synteny_df['block_length_first'].mean(),
                    'mean_length_second': synteny_df['block_length_second'].mean(),
                    'max_length_first': synteny_df['block_length_first'].max(),
                    'max_length_second': synteny_df['block_length_second'].max()
                }
        
        if len(mapping_df) > 0:
            # Correlation statistics
            if 'position_correlation' in mapping_df.columns:
                stats['correlation_stats'] = {
                    'mean_correlation': mapping_df['position_correlation'].mean(),
                    'std_correlation': mapping_df['position_correlation'].std(),
                    'strong_positive_correlations': (mapping_df['position_correlation'] >= self.correlation_threshold).sum(),
                    'strong_negative_correlations': (mapping_df['position_correlation'] <= -self.correlation_threshold).sum(),
                    'weak_correlations': (mapping_df['position_correlation'].abs() < self.correlation_threshold).sum()
                }
            
            # Strand consistency statistics
            if 'strand_consistency' in mapping_df.columns:
                stats['strand_stats'] = {
                    'mean_strand_consistency': mapping_df['strand_consistency'].mean(),
                    'std_strand_consistency': mapping_df['strand_consistency'].std(),
                    'consistent_pairs': (mapping_df['strand_consistency'] >= self.strand_consistency_threshold).sum(),
                    'inconsistent_pairs': (mapping_df['strand_consistency'] < self.strand_consistency_threshold).sum()
                }
            
            # Gene count statistics
            if 'gene_count' in mapping_df.columns:
                stats['gene_count_stats'] = {
                    'total_genes_mapped': mapping_df['gene_count'].sum(),
                    'mean_genes_per_pair': mapping_df['gene_count'].mean(),
                    'max_genes_per_pair': mapping_df['gene_count'].max(),
                    'pairs_with_many_genes': (mapping_df['gene_count'] >= 10).sum()
                }
        
        return stats