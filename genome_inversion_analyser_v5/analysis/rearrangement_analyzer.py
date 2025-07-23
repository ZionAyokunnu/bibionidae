# =============================================================================
# Rearrangement Analyzer (rearrangement_analyzer.py)
# =============================================================================

"""
Comprehensive chromosome rearrangement detection system.
Identifies splits, fusions, translocations, and complex rearrangements.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Any, Set, Tuple, Optional
from collections import defaultdict, Counter

from ..logger import get_logger

logger = get_logger()

class RearrangementAnalyzer:
    """
    Comprehensive chromosome rearrangement analyzer that detects various types
    of chromosomal structural variations including splits, fusions, and translocations.
    """
    
    def __init__(self, config):
        """
        Initialize rearrangement analyzer.
        
        Args:
            config: Configuration object with rearrangement parameters
        """
        self.config = config
        self.enable_translocation_detection = config.get('enable_translocation_detection', True)
        self.min_genes_for_rearrangement = config.get('min_genes_for_rearrangement', 3)
        self.rearrangement_confidence_threshold = config.get('rearrangement_confidence_threshold', 0.5)
        
        logger.info("Rearrangement analyzer initialized")
        logger.info(f"  Translocation detection: {'enabled' if self.enable_translocation_detection else 'disabled'}")
        logger.info(f"  Min genes for rearrangement: {self.min_genes_for_rearrangement}")
    
    def analyze_chromosome_rearrangements(self, ortholog_df: pd.DataFrame) -> pd.DataFrame:
        """
        Analyze chromosome rearrangements from ortholog data.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            
        Returns:
            DataFrame with detected rearrangements
        """
        if len(ortholog_df) == 0:
            logger.warning("No ortholog pairs found for rearrangement analysis")
            return pd.DataFrame()
        
        logger.info(f"Analyzing chromosome rearrangements from {len(ortholog_df)} ortholog pairs...")
        
        rearrangements = []
        
        # Detect chromosome splits
        splits = self._detect_chromosome_splits(ortholog_df)
        rearrangements.extend(splits)
        
        # Detect chromosome fusions
        fusions = self._detect_chromosome_fusions(ortholog_df)
        rearrangements.extend(fusions)
        
        # Detect translocations if enabled
        if self.enable_translocation_detection:
            translocations = self._detect_translocations(ortholog_df)
            rearrangements.extend(translocations)
        
        # Detect complex rearrangements
        complex_rearrangements = self._detect_complex_rearrangements(ortholog_df)
        rearrangements.extend(complex_rearrangements)
        
        rearrangement_df = pd.DataFrame(rearrangements)
        
        # Log detection results
        self._log_rearrangement_results(rearrangement_df, ortholog_df)
        
        return rearrangement_df
    
    def _detect_chromosome_splits(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect chromosome split events (1 → many)."""
        splits = []
        
        # Group by first chromosome to find splits
        for first_chr in ortholog_df['first_chr'].unique():
            first_genes = ortholog_df[ortholog_df['first_chr'] == first_chr]
            target_chromosomes = first_genes['second_chr'].value_counts()
            
            if len(target_chromosomes) > 1:  # This chromosome maps to multiple targets
                # Check if this represents a genuine split
                total_genes = len(first_genes)
                
                if total_genes >= self.min_genes_for_rearrangement:
                    confidence = self._calculate_rearrangement_confidence(first_genes, 'split')
                    
                    if confidence >= self.rearrangement_confidence_threshold:
                        split = {
                            'type': 'chromosome_split',
                            'source_chr': first_chr,
                            'target_chrs': target_chromosomes.index.tolist(),
                            'source_genome': 'first',
                            'target_genome': 'second',
                            'gene_counts': target_chromosomes.values.tolist(),
                            'total_genes': total_genes,
                            'confidence': confidence,
                            'coverage': sum(target_chromosomes.values) / total_genes,
                            'split_ratio': len(target_chromosomes),
                            'largest_fragment': target_chromosomes.max(),
                            'smallest_fragment': target_chromosomes.min(),
                            'fragment_balance': target_chromosomes.min() / target_chromosomes.max()
                        }
                        splits.append(split)
        
        return splits
    
    def _detect_chromosome_fusions(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect chromosome fusion events (many → 1)."""
        fusions = []
        
        # Group by second chromosome to find fusions
        for second_chr in ortholog_df['second_chr'].unique():
            second_genes = ortholog_df[ortholog_df['second_chr'] == second_chr]
            source_chromosomes = second_genes['first_chr'].value_counts()
            
            if len(source_chromosomes) > 1:  # Multiple sources map to this chromosome
                total_genes = len(second_genes)
                
                if total_genes >= self.min_genes_for_rearrangement:
                    confidence = self._calculate_rearrangement_confidence(second_genes, 'fusion')
                    
                    if confidence >= self.rearrangement_confidence_threshold:
                        fusion = {
                            'type': 'chromosome_fusion',
                            'source_chrs': source_chromosomes.index.tolist(),
                            'target_chr': second_chr,
                            'source_genome': 'first',
                            'target_genome': 'second',
                            'gene_counts': source_chromosomes.values.tolist(),
                            'total_genes': total_genes,
                            'confidence': confidence,
                            'coverage': sum(source_chromosomes.values) / total_genes,
                            'fusion_ratio': len(source_chromosomes),
                            'largest_contributor': source_chromosomes.max(),
                            'smallest_contributor': source_chromosomes.min(),
                            'contribution_balance': source_chromosomes.min() / source_chromosomes.max()
                        }
                        fusions.append(fusion)
        
        return fusions
    
    def _detect_translocations(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect reciprocal translocation events."""
        translocations = []
        
        # Look for reciprocal patterns (A→B and B→A)
        chromosome_pairs = set()
        
        for _, row in ortholog_df.iterrows():
            pair = tuple(sorted([row['first_chr'], row['second_chr']]))
            chromosome_pairs.add(pair)
        
        # Analyze each chromosome pair for translocation patterns
        for chr1, chr2 in chromosome_pairs:
            translocation = self._analyze_translocation_pattern(ortholog_df, chr1, chr2)
            if translocation:
                translocations.append(translocation)
        
        return translocations
    
    def _analyze_translocation_pattern(self, ortholog_df: pd.DataFrame, 
                                     chr1: str, chr2: str) -> Optional[Dict]:
        """Analyze a specific chromosome pair for translocation patterns."""
        
        # Get all genes involving these two chromosomes
        involved_genes = ortholog_df[
            ((ortholog_df['first_chr'] == chr1) | (ortholog_df['first_chr'] == chr2)) &
            ((ortholog_df['second_chr'] == chr1) | (ortholog_df['second_chr'] == chr2))
        ]
        
        if len(involved_genes) < self.min_genes_for_rearrangement:
            return None
        
        # Count mapping patterns
        pattern_counts = defaultdict(int)
        for _, gene in involved_genes.iterrows():
            pattern = f"{gene['first_chr']}→{gene['second_chr']}"
            pattern_counts[pattern] += 1
        
        # Check for reciprocal pattern
        pattern1 = f"{chr1}→{chr2}"
        pattern2 = f"{chr2}→{chr1}"
        
        if pattern1 in pattern_counts and pattern2 in pattern_counts:
            # Potential reciprocal translocation
            confidence = self._calculate_translocation_confidence(involved_genes, pattern_counts)
            
            if confidence >= self.rearrangement_confidence_threshold:
                return {
                    'type': 'reciprocal_translocation',
                    'chromosome_1': chr1,
                    'chromosome_2': chr2,
                    'pattern_1_count': pattern_counts[pattern1],
                    'pattern_2_count': pattern_counts[pattern2],
                    'total_genes': len(involved_genes),
                    'confidence': confidence,
                    'reciprocity_balance': min(pattern_counts[pattern1], pattern_counts[pattern2]) / max(pattern_counts[pattern1], pattern_counts[pattern2]),
                    'pattern_counts': dict(pattern_counts)
                }
        
        return None
    
    def _detect_complex_rearrangements(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Detect complex multi-chromosome rearrangements."""
        complex_rearrangements = []
        
        # Analyze chromosome mapping complexity
        chromosome_complexity = self._analyze_chromosome_complexity(ortholog_df)
        
        for complexity_info in chromosome_complexity:
            if complexity_info['complexity_score'] > 2.0:  # Threshold for complex rearrangement
                confidence = self._calculate_complexity_confidence(complexity_info)
                
                if confidence >= self.rearrangement_confidence_threshold:
                    complex_rearr = {
                        'type': 'complex_rearrangement',
                        'involved_chromosomes': complexity_info['chromosomes'],
                        'complexity_score': complexity_info['complexity_score'],
                        'total_genes': complexity_info['total_genes'],
                        'mapping_patterns': complexity_info['patterns'],
                        'confidence': confidence,
                        'rearrangement_index': len(complexity_info['chromosomes']) / 2.0  # Normalized complexity
                    }
                    complex_rearrangements.append(complex_rearr)
        
        return complex_rearrangements
    
    def _analyze_chromosome_complexity(self, ortholog_df: pd.DataFrame) -> List[Dict]:
        """Analyze the complexity of chromosome mappings."""
        complexity_results = []
        
        # Group chromosomes that are involved in multiple mappings
        chromosome_networks = defaultdict(set)
        
        for _, row in ortholog_df.iterrows():
            first_chr = row['first_chr']
            second_chr = row['second_chr']
            chromosome_networks[first_chr].add(second_chr)
            chromosome_networks[second_chr].add(first_chr)
        
        # Find highly connected chromosome groups
        processed = set()
        
        for chr_name, connections in chromosome_networks.items():
            if chr_name in processed or len(connections) <= 1:
                continue
            
            # Find connected component
            component = {chr_name}
            to_process = list(connections)
            
            while to_process:
                current = to_process.pop()
                if current not in component:
                    component.add(current)
                    to_process.extend(chromosome_networks[current] - component)
            
            if len(component) > 2:  # Complex if more than 2 chromosomes involved
                # Calculate complexity metrics
                total_genes = len(ortholog_df[
                    ortholog_df['first_chr'].isin(component) | 
                    ortholog_df['second_chr'].isin(component)
                ])
                
                # Count unique mapping patterns
                patterns = set()
                for _, row in ortholog_df.iterrows():
                    if row['first_chr'] in component or row['second_chr'] in component:
                        patterns.add(f"{row['first_chr']}→{row['second_chr']}")
                
                complexity_score = len(patterns) / len(component)  # Average patterns per chromosome
                
                complexity_results.append({
                    'chromosomes': list(component),
                    'complexity_score': complexity_score,
                    'total_genes': total_genes,
                    'patterns': list(patterns),
                    'pattern_count': len(patterns)
                })
                
                processed.update(component)
        
        return complexity_results
    
    def _calculate_rearrangement_confidence(self, genes: pd.DataFrame, rearr_type: str) -> float:
        """Calculate confidence for a rearrangement event."""
        factors = []
        
        # Gene count factor (more genes = higher confidence)
        count_factor = min(1.0, len(genes) / 10.0)
        factors.append(count_factor)
        
        # Quality factor (average similarity of involved genes)
        if 'similarity' in genes.columns:
            similarity_factor = genes['similarity'].mean()
            factors.append(similarity_factor)
        
        # Confidence factor (average confidence of gene alignments)
        if 'confidence' in genes.columns:
            conf_factor = genes['confidence'].mean()
            factors.append(conf_factor)
        
        # Distribution factor (how evenly genes are distributed across targets)
        if rearr_type in ['split', 'fusion']:
            target_col = 'second_chr' if rearr_type == 'split' else 'first_chr'
            target_counts = genes[target_col].value_counts()
            # Higher entropy (more even distribution) = higher confidence
            entropy = -sum((count / len(genes)) * np.log2(count / len(genes)) for count in target_counts)
            max_entropy = np.log2(len(target_counts))
            distribution_factor = entropy / max_entropy if max_entropy > 0 else 0
            factors.append(distribution_factor)
        
        return np.mean(factors) if factors else 0.5
    
    def _calculate_translocation_confidence(self, genes: pd.DataFrame, 
                                          pattern_counts: Dict[str, int]) -> float:
        """Calculate confidence for a translocation event."""
        factors = []
        
        # Gene count factor
        count_factor = min(1.0, len(genes) / 10.0)
        factors.append(count_factor)
        
        # Reciprocity balance factor
        counts = list(pattern_counts.values())
        if len(counts) >= 2:
            balance_factor = min(counts) / max(counts)
            factors.append(balance_factor)
        
        # Quality factor
        if 'similarity' in genes.columns:
            similarity_factor = genes['similarity'].mean()
            factors.append(similarity_factor)
        
        return np.mean(factors) if factors else 0.5
    
    def _calculate_complexity_confidence(self, complexity_info: Dict) -> float:
        """Calculate confidence for a complex rearrangement."""
        factors = []
        
        # Complexity factor (moderate complexity is more confident than extreme)
        complexity_score = complexity_info['complexity_score']
        # Optimal complexity around 2-3 patterns per chromosome
        if 2.0 <= complexity_score <= 4.0:
            complexity_factor = 1.0
        else:
            complexity_factor = 1.0 / (1.0 + abs(complexity_score - 3.0) * 0.2)
        factors.append(complexity_factor)
        
        # Gene count factor
        count_factor = min(1.0, complexity_info['total_genes'] / 20.0)
        factors.append(count_factor)
        
        # Pattern diversity factor
        pattern_diversity = complexity_info['pattern_count'] / len(complexity_info['chromosomes'])
        diversity_factor = min(1.0, pattern_diversity / 2.0)
        factors.append(diversity_factor)
        
        return np.mean(factors) if factors else 0.5
    
    def _log_rearrangement_results(self, rearrangement_df: pd.DataFrame, 
                                 ortholog_df: pd.DataFrame):
        """Log comprehensive rearrangement detection results."""
        logger.info("  Chromosome rearrangement results:")
        logger.info(f"    Total rearrangements detected: {len(rearrangement_df)}")
        
        if len(rearrangement_df) > 0:
            # Type distribution
            type_counts = rearrangement_df['type'].value_counts()
            logger.info(f"    Rearrangement types: {type_counts.to_dict()}")
            
            # Quality metrics
            if 'confidence' in rearrangement_df.columns:
                avg_confidence = rearrangement_df['confidence'].mean()
                high_conf = (rearrangement_df['confidence'] >= 0.8).sum()
                logger.info(f"    Average confidence: {avg_confidence:.3f}")
                logger.info(f"    High confidence events: {high_conf}")
            
            # Gene involvement
            total_genes_involved = 0
            for _, rearr in rearrangement_df.iterrows():
                if 'total_genes' in rearr:
                    total_genes_involved += rearr['total_genes']
            
            involvement_rate = total_genes_involved / len(ortholog_df) if len(ortholog_df) > 0 else 0
            logger.info(f"    Gene involvement: {involvement_rate:.1%} ({total_genes_involved}/{len(ortholog_df)} genes)")
    
    def get_rearrangement_statistics(self, rearrangement_df: pd.DataFrame) -> Dict[str, Any]:
        """
        Get comprehensive rearrangement statistics.
        
        Args:
            rearrangement_df: DataFrame with detected rearrangements
            
        Returns:
            Dictionary with rearrangement statistics
        """
        if len(rearrangement_df) == 0:
            return {'total_rearrangements': 0}
        
        stats = {
            'total_rearrangements': len(rearrangement_df),
            'type_distribution': rearrangement_df['type'].value_counts().to_dict()
        }
        
        # Confidence statistics
        if 'confidence' in rearrangement_df.columns:
            stats['confidence_stats'] = {
                'mean': rearrangement_df['confidence'].mean(),
                'std': rearrangement_df['confidence'].std(),
                'min': rearrangement_df['confidence'].min(),
                'max': rearrangement_df['confidence'].max(),
                'high_confidence_count': (rearrangement_df['confidence'] >= 0.8).sum()
            }
        
        # Type-specific statistics
        for rearr_type in rearrangement_df['type'].unique():
            type_data = rearrangement_df[rearrangement_df['type'] == rearr_type]
            stats[f'{rearr_type}_stats'] = {
                'count': len(type_data),
                'avg_confidence': type_data['confidence'].mean() if 'confidence' in type_data.columns else None
            }
            
            # Add type-specific metrics
            if rearr_type == 'chromosome_split' and 'split_ratio' in type_data.columns:
                stats[f'{rearr_type}_stats']['avg_split_ratio'] = type_data['split_ratio'].mean()
            elif rearr_type == 'chromosome_fusion' and 'fusion_ratio' in type_data.columns:
                stats[f'{rearr_type}_stats']['avg_fusion_ratio'] = type_data['fusion_ratio'].mean()
            elif rearr_type == 'complex_rearrangement' and 'complexity_score' in type_data.columns:
                stats[f'{rearr_type}_stats']['avg_complexity'] = type_data['complexity_score'].mean()
        
        return stats