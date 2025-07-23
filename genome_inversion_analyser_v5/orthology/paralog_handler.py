
# =============================================================================
# Paralog Handler (paralog_handler.py)
# =============================================================================

"""
Advanced paralog detection and handling for complex ortholog relationships.
Identifies and manages paralogous relationships within and between genomes.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Set, Tuple, Any
from collections import defaultdict, Counter
from ..logger import get_logger

logger = get_logger()

class ParalogHandler:
    """
    Advanced paralog detection and handling system.
    Manages complex many-to-many orthologous relationships and paralog clusters.
    """
    
    def __init__(self, config):
        """
        Initialize paralog handler.
        
        Args:
            config: Configuration object with paralog settings
        """
        self.config = config
        self.enable_paralog_detection = config.get('enable_paralog_detection', True)
        self.enable_paralog_mapping = config.get('enable_paralog_ortholog_mapping', True)
        self.max_paralogs = config.get('max_paralogs_per_busco', 3)
        self.paralog_similarity_threshold = config.get('paralog_similarity_threshold', 0.8)
        
        logger.info("Paralog handler initialized")
        logger.info(f"  Paralog detection: {'enabled' if self.enable_paralog_detection else 'disabled'}")
        logger.info(f"  Paralog mapping: {'enabled' if self.enable_paralog_mapping else 'disabled'}")
        logger.info(f"  Max paralogs per BUSCO: {self.max_paralogs}")
    
    def detect_paralogs_in_busco_data(self, busco_df: pd.DataFrame) -> pd.DataFrame:
        """
        Detect and annotate paralogs in BUSCO data.
        
        Args:
            busco_df: DataFrame with BUSCO data
            
        Returns:
            DataFrame with paralog annotations
        """
        if not self.enable_paralog_detection:
            return busco_df
        
        logger.info(f"Detecting paralogs in {len(busco_df)} BUSCO entries...")
        
        # Count occurrences of each BUSCO ID
        busco_counts = busco_df['busco_id'].value_counts()
        paralogs = busco_counts[busco_counts > 1]
        
        # Add paralog annotations
        busco_df_annotated = busco_df.copy()
        busco_df_annotated['is_paralog'] = busco_df_annotated['busco_id'].isin(paralogs.index)
        busco_df_annotated['paralog_count'] = busco_df_annotated['busco_id'].map(busco_counts)
        
        # Initialize paralog metadata
        busco_df_annotated['paralog_rank'] = 1
        busco_df_annotated['paralog_cluster'] = busco_df_annotated['busco_id'] + '_singleton'
        
        # Process each BUSCO ID with paralogs
        for busco_id in paralogs.index:
            paralog_group = busco_df_annotated[busco_df_annotated['busco_id'] == busco_id].copy()
            
            # Rank paralogs using multiple criteria
            ranks, clusters = self._rank_and_cluster_paralogs(paralog_group)
            
            # Update main dataframe
            busco_df_annotated.loc[paralog_group.index, 'paralog_rank'] = ranks
            busco_df_annotated.loc[paralog_group.index, 'paralog_cluster'] = clusters
        
        # Log paralog statistics
        self._log_paralog_statistics(busco_df_annotated, paralogs)
        
        return busco_df_annotated
    
    def _rank_and_cluster_paralogs(self, paralog_group: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
        """
        Rank and cluster paralogs within a BUSCO group.
        
        Args:
            paralog_group: DataFrame with paralogs for single BUSCO ID
            
        Returns:
            Tuple of (ranks, cluster_ids)
        """
        busco_id = paralog_group['busco_id'].iloc[0]
        
        # Multi-criteria ranking
        ranking_criteria = []
        weights = []
        
        # Primary: BUSCO score (if available)
        if 'score' in paralog_group.columns and paralog_group['score'].notna().any():
            ranking_criteria.append(paralog_group['score'].fillna(0))
            weights.append(0.4)
        
        # Secondary: sequence length
        if 'gene_length' in paralog_group.columns and paralog_group['gene_length'].notna().any():
            ranking_criteria.append(paralog_group['gene_length'].fillna(0))
            weights.append(0.3)
        elif 'length' in paralog_group.columns and paralog_group['length'].notna().any():
            ranking_criteria.append(paralog_group['length'].fillna(0))
            weights.append(0.3)
        
        # Tertiary: genomic position (earlier positions get higher rank)
        if 'gene_start' in paralog_group.columns and paralog_group['gene_start'].notna().any():
            ranking_criteria.append(-paralog_group['gene_start'].fillna(float('inf')))
            weights.append(0.2)
        
        # Quaternary: chromosome preference (smaller chromosome names often more important)
        if 'sequence_id' in paralog_group.columns:
            chr_scores = paralog_group['sequence_id'].apply(self._calculate_chromosome_score)
            ranking_criteria.append(chr_scores)
            weights.append(0.1)
        
        # Calculate combined ranking score
        if ranking_criteria:
            # Normalize weights
            weights = [w / sum(weights) for w in weights]
            
            combined_score = sum(w * criteria for w, criteria in zip(weights, ranking_criteria))
            ranks = combined_score.rank(ascending=False, method='first')
        else:
            # Fallback: rank by dataframe order
            ranks = pd.Series(range(1, len(paralog_group) + 1), index=paralog_group.index)
        
        # Create cluster IDs
        clusters = pd.Series([f"{busco_id}_cluster"] * len(paralog_group), index=paralog_group.index)
        
        return ranks, clusters
    
    def _calculate_chromosome_score(self, chromosome_id: str) -> float:
        """Calculate chromosome preference score (smaller/main chromosomes get higher scores)."""
        chr_str = str(chromosome_id).lower()
        
        # Prefer main chromosomes (those with shorter names or numbers)
        if any(main_chr in chr_str for main_chr in ['chr1', 'chromosome1', 'chr01', 'scaffold1']):
            return 10.0
        elif any(chr_str.startswith(prefix) for prefix in ['chr', 'chromosome']):
            try:
                # Extract number and prefer smaller numbers
                import re
                numbers = re.findall(r'\d+', chr_str)
                if numbers:
                    return 10.0 / (int(numbers[0]) + 1)
            except:
                pass
        
        # Default score for other chromosomes
        return 1.0 / (len(chr_str) + 1)
    
    def _log_paralog_statistics(self, busco_df: pd.DataFrame, paralogs: pd.Series):
        """Log detailed paralog statistics."""
        total_paralogs = len(busco_df[busco_df['is_paralog']])
        
        logger.info("  Paralog detection results:")
        logger.info(f"    BUSCO IDs with paralogs: {len(paralogs)}")
        logger.info(f"    Total paralogous genes: {total_paralogs}")
        logger.info(f"    Max paralogs for single BUSCO: {paralogs.max()}")
        
        # Distribution of paralog counts
        paralog_distribution = paralogs.value_counts().sort_index()
        logger.info("    Paralog count distribution:")
        for count, frequency in paralog_distribution.items():
            logger.info(f"      {count} paralogs: {frequency} BUSCO IDs")
    
    def analyze_ortholog_paralog_relationships(self, ortholog_df: pd.DataFrame) -> pd.DataFrame:
        """
        Analyze complex ortholog-paralog relationships.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            
        Returns:
            DataFrame with paralog relationship analysis
        """
        if not self.enable_paralog_mapping or len(ortholog_df) == 0:
            return pd.DataFrame()
        
        logger.info(f"Analyzing ortholog-paralog relationships in {len(ortholog_df)} pairs...")
        
        paralog_relationships = []
        
        # Group by BUSCO ID to find complex relationships
        for busco_id, group in ortholog_df.groupby('busco_id'):
            if len(group) > 1:
                # Multiple ortholog pairs for same BUSCO - analyze relationships
                relationships = self._analyze_busco_paralog_group(busco_id, group)
                paralog_relationships.extend(relationships)
        
        # Analyze within-genome paralogs
        within_genome_paralogs = self._analyze_within_genome_paralogs(ortholog_df)
        paralog_relationships.extend(within_genome_paralogs)
        
        paralog_df = pd.DataFrame(paralog_relationships)
        
        if len(paralog_df) > 0:
            logger.info(f"  Found {len(paralog_df)} complex paralog relationships")
        
        return paralog_df
    
    def _analyze_busco_paralog_group(self, busco_id: str, group: pd.DataFrame) -> List[Dict]:
        """Analyze paralog relationships within a BUSCO group."""
        relationships = []
        
        # Check for many-to-many relationships
        first_chrs = group['first_chr'].unique()
        second_chrs = group['second_chr'].unique()
        
        if len(first_chrs) > 1 or len(second_chrs) > 1:
            relationship = {
                'busco_id': busco_id,
                'relationship_type': 'many_to_many',
                'first_paralogs': len(first_chrs),
                'second_paralogs': len(second_chrs),
                'total_pairs': len(group),
                'avg_similarity': group['similarity'].mean(),
                'avg_confidence': group['confidence'].mean(),
                'chromosomes_first': list(first_chrs),
                'chromosomes_second': list(second_chrs)
            }
            relationships.append(relationship)
        
        return relationships
    
    def _analyze_within_genome_paralogs(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Analyze within-genome paralog patterns."""
        relationships = []
        
        # Analyze first genome paralogs
        first_paralogs = self._find_within_genome_paralogs(
            ortholog_df, 'first_chr', 'first_paralog_count'
        )
        relationships.extend(first_paralogs)
        
        # Analyze second genome paralogs
        second_paralogs = self._find_within_genome_paralogs(
            ortholog_df, 'second_chr', 'second_paralog_count'
        )
        relationships.extend(second_paralogs)
        
        return relationships
    
    def _find_within_genome_paralogs(self, ortholog_df: pd.DataFrame, 
                                   chr_col: str, count_col: str) -> List[Dict]:
        """Find within-genome paralog patterns."""
        relationships = []
        
        if count_col not in ortholog_df.columns:
            return relationships
        
        # Find genes with paralogs
        paralog_genes = ortholog_df[ortholog_df[count_col] > 1]
        
        for busco_id, group in paralog_genes.groupby('busco_id'):
            if len(group) > 0:
                relationship = {
                    'busco_id': busco_id,
                    'relationship_type': 'within_genome_paralogs',
                    'genome': 'first' if 'first' in chr_col else 'second',
                    'paralog_count': group[count_col].iloc[0],
                    'chromosomes': list(group[chr_col].unique()),
                    'avg_similarity': group['similarity'].mean(),
                    'avg_confidence': group['confidence'].mean()
                }
                relationships.append(relationship)
        
        return relationships
    
    def get_paralog_statistics(self, busco_df: pd.DataFrame = None, 
                             ortholog_df: pd.DataFrame = None,
                             paralog_df: pd.DataFrame = None) -> Dict[str, Any]:
        """
        Get comprehensive paralog statistics.
        
        Args:
            busco_df: Optional BUSCO DataFrame with paralog annotations
            ortholog_df: Optional ortholog DataFrame
            paralog_df: Optional paralog relationship DataFrame
            
        Returns:
            Dictionary with paralog statistics
        """
        stats = {}
        
        # BUSCO-level paralog statistics
        if busco_df is not None and 'is_paralog' in busco_df.columns:
            total_genes = len(busco_df)
            paralog_genes = len(busco_df[busco_df['is_paralog']])
            unique_buscos = busco_df['busco_id'].nunique()
            paralog_buscos = busco_df[busco_df['is_paralog']]['busco_id'].nunique()
            
            stats['busco_level'] = {
                'total_genes': total_genes,
                'paralog_genes': paralog_genes,
                'paralog_rate': paralog_genes / total_genes if total_genes > 0 else 0,
                'unique_buscos': unique_buscos,
                'paralog_buscos': paralog_buscos,
                'busco_paralog_rate': paralog_buscos / unique_buscos if unique_buscos > 0 else 0
            }
            
            # Paralog count distribution
            if 'paralog_count' in busco_df.columns:
                paralog_counts = busco_df[busco_df['is_paralog']]['paralog_count'].value_counts()
                stats['busco_level']['paralog_distribution'] = paralog_counts.to_dict()
        
        # Ortholog-level paralog statistics
        if ortholog_df is not None:
            stats['ortholog_level'] = {
                'total_pairs': len(ortholog_df),
                'unique_buscos': ortholog_df['busco_id'].nunique()
            }
            
            # Check for many-to-many relationships
            busco_pair_counts = ortholog_df['busco_id'].value_counts()
            many_to_many = busco_pair_counts[busco_pair_counts > 1]
            
            stats['ortholog_level']['many_to_many_buscos'] = len(many_to_many)
            stats['ortholog_level']['many_to_many_rate'] = len(many_to_many) / len(busco_pair_counts) if len(busco_pair_counts) > 0 else 0
        
        # Complex relationship statistics
        if paralog_df is not None and len(paralog_df) > 0:
            relationship_types = paralog_df['relationship_type'].value_counts()
            stats['complex_relationships'] = {
                'total_relationships': len(paralog_df),
                'relationship_types': relationship_types.to_dict()
            }
        
        return stats