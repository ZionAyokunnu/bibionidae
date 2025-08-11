"""
Synteny analysis module for the Genome Inversion Analyzer
Handles synteny block detection, chromosome rearrangements, and inversion analysis
"""

import pandas as pd
import numpy as np
import logging
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def analyze_enhanced_synteny_blocks(ortholog_df, config):
    """Enhanced synteny analysis with improved block detection"""
    synteny_blocks = []
    chromosome_mappings = []
    
    if len(ortholog_df) == 0:
        logger.warning("No ortholog pairs found for synteny analysis")
        return pd.DataFrame(synteny_blocks)ƒ, pd.DataFrame(chromosome_mappings)
    
    for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
        min_genes = config.get('base_min_genes_per_chromosome', 3)
        if config.get('enable_small_synteny_blocks', False):
            min_genes = config.get('micro_synteny_block_size', 1)
            
        if len(group) >= min_genes:
            group_sorted = group.sort_values('first_start')
            
            if len(group_sorted) > 1:
                try:
                    correlation, p_value = pearsonr(group_sorted['first_start'], group_sorted['second_start'])
                except:
                    correlation, p_value = 0.0, 1.0
            else:
                correlation, p_value = 1.0, 0.0
            
            strand_consistency = (group_sorted['first_strand'] == group_sorted['second_strand']).mean()
            
            synteny_blocks.append({
                'first_chr': first_chr,
                'second_chr': second_chr,
                'block_size': len(group),
                'start_gene': group_sorted.iloc[0]['busco_id'],
                'end_gene': group_sorted.iloc[-1]['busco_id'],
                'position_correlation': correlation,
                'strand_consistency': strand_consistency,
                'synteny_type': classify_enhanced_synteny_type(correlation, strand_consistency, config),
                'confidence': calculate_synteny_confidence(group_sorted, config)
            })
            
            chromosome_mappings.append({
                'first_chr': first_chr,
                'second_chr': second_chr,
                'gene_count': len(group),
                'position_correlation': correlation,
                'strand_consistency': strand_consistency,
                'p_value': p_value
            })
    
    return pd.DataFrame(synteny_blocks), pd.DataFrame(chromosome_mappings)


def classify_enhanced_synteny_type(correlation, strand_consistency, config):
    """Enhanced synteny type classification"""
    correlation_threshold = config.get('base_synteny_correlation_threshold', 0.5)
    strand_threshold = config.get('strand_consistency_threshold', 0.7)
    
    if correlation > correlation_threshold:
        if strand_consistency > strand_threshold:
            return 'colinear'
        else:
            return 'colinear_inverted'
    elif correlation < -correlation_threshold:
        if strand_consistency < (1 - strand_threshold):
            return 'inverted'
        else:
            return 'inverted_mixed'
    else:
        return 'rearranged'


def calculate_synteny_confidence(group, config):
    """Calculate confidence score for synteny block"""
    if config.get('enable_synteny_confidence_scoring', False):
        factors = []
        
        # Size factor
        size_factor = min(1.0, len(group) / 10)
        factors.append(size_factor)
        
        # Similarity factor
        if 'similarity' in group.columns:
            sim_factor = group['similarity'].mean()
            factors.append(sim_factor)
        
        # Confidence factor
        if 'confidence' in group.columns:
            conf_factor = group['confidence'].mean()
            factors.append(conf_factor)
        
        return np.mean(factors) if factors else 0.5
    else:
        return 1.0


def analyze_enhanced_chromosome_rearrangements(ortholog_df, config):
    """Enhanced chromosome rearrangement analysis"""
    rearrangements = []
    
    if len(ortholog_df) == 0:
        logger.warning("No ortholog pairs found for rearrangement analysis")
        return pd.DataFrame(rearrangements)
    
    # Enhanced split detection
    for first_chr in ortholog_df['first_chr'].unique():
        first_genes = ortholog_df[ortholog_df['first_chr'] == first_chr]
        target_chromosomes = first_genes['second_chr'].value_counts()
        
        if len(target_chromosomes) > 1:
            total_genes = len(first_genes)
            confidence = calculate_rearrangement_confidence(first_genes, 'split', config)
            
            rearrangements.append({
                'type': 'chromosome_split',
                'first_chr': first_chr,
                'second_chrs': target_chromosomes.index.tolist(),
                'gene_counts': target_chromosomes.values.tolist(),
                'total_genes': total_genes,
                'confidence': confidence,
                'coverage': sum(target_chromosomes.values) / total_genes
            })
    
    # Enhanced fusion detection
    for second_chr in ortholog_df['second_chr'].unique():
        second_genes = ortholog_df[ortholog_df['second_chr'] == second_chr]
        source_chromosomes = second_genes['first_chr'].value_counts()
        
        if len(source_chromosomes) > 1:
            total_genes = len(second_genes)
            confidence = calculate_rearrangement_confidence(second_genes, 'fusion', config)
            
            rearrangements.append({
                'type': 'chromosome_fusion',
                'first_chrs': source_chromosomes.index.tolist(),
                'second_chr': second_chr,
                'gene_counts': source_chromosomes.values.tolist(),
                'total_genes': total_genes,
                'confidence': confidence,
                'coverage': sum(source_chromosomes.values) / total_genes
            })
    
    return pd.DataFrame(rearrangements)


def calculate_rearrangement_confidence(genes, rearr_type, config):
    """Calculate confidence for rearrangement detection"""
    if config.get('enable_rearrangement_scoring', False):
        factors = []
        
        # Gene count factor
        count_factor = min(1.0, len(genes) / 10)
        factors.append(count_factor)
        
        # Similarity factor
        if 'similarity' in genes.columns:
            sim_factor = genes['similarity'].mean()
            factors.append(sim_factor)
        
        # Mapping confidence factor
        if 'confidence' in genes.columns:
            conf_factor = genes['confidence'].mean()
            factors.append(conf_factor)
        
        return np.mean(factors) if factors else 0.5
    else:
        return 1.0


def analyze_enhanced_inversions(synteny_df, ortholog_df, config):
    """Enhanced inversion analysis"""
    inversions = []
    
    if len(synteny_df) == 0 or len(ortholog_df) == 0:
        logger.warning("No synteny blocks or ortholog pairs found for inversion analysis")
        return pd.DataFrame(inversions)
    
    if config.get('enable_micro_inversions', False):
        # Analyze inversions in all synteny types
        for _, block in synteny_df.iterrows():
            block_genes = ortholog_df[
                (ortholog_df['first_chr'] == block['first_chr']) & 
                (ortholog_df['second_chr'] == block['second_chr'])
            ].sort_values('first_start')
            
            if len(block_genes) >= config.get('base_min_inversion_size', 2):
                inversion_regions = identify_enhanced_inversions(block_genes, config)
                
                for region in inversion_regions:
                    inversions.append({
                        'first_chr': block['first_chr'],
                        'second_chr': block['second_chr'],
                        'start_gene': region['start_gene'],
                        'end_gene': region['end_gene'],
                        'size_genes': region['size'],
                        'inversion_type': region['type'],
                        'strand_pattern': region['strand_pattern'],
                        'confidence': region.get('confidence', 0.5),
                        'synteny_context': block['synteny_type']
                    })
    
    # Single gene inversions if enabled
    if config.get('enable_single_gene_inversions', False):
        single_inversions = detect_single_gene_inversions(ortholog_df, config)
        inversions.extend(single_inversions)
    
    return pd.DataFrame(inversions)


def identify_enhanced_inversions(genes, config):
    """Enhanced inversion identification with confidence scoring"""
    regions = []
    genes_sorted = genes.sort_values('first_start')
    
    current_region = []
    for _, gene in genes_sorted.iterrows():
        strand_flipped = gene['first_strand'] != gene['second_strand']
        
        if strand_flipped:
            current_region.append(gene)
        else:
            if len(current_region) >= config.get('base_min_inversion_size', 2):
                confidence = calculate_inversion_confidence(current_region, config)
                
                regions.append({
                    'start_gene': current_region[0]['busco_id'],
                    'end_gene': current_region[-1]['busco_id'],
                    'size': len(current_region),
                    'type': 'strand_inversion',
                    'strand_pattern': 'consistent_flip',
                    'confidence': confidence
                })
            current_region = []
    
    # Handle final region
    if len(current_region) >= config.get('base_min_inversion_size', 2):
        confidence = calculate_inversion_confidence(current_region, config)
        regions.append({
            'start_gene': current_region[0]['busco_id'],
            'end_gene': current_region[-1]['busco_id'],
            'size': len(current_region),
            'type': 'strand_inversion',
            'strand_pattern': 'consistent_flip',
            'confidence': confidence
        })
    
    return regions


def calculate_inversion_confidence(genes, config):
    """Calculate confidence for inversion detection"""
    if config.get('enable_inversion_confidence', False):
        factors = []
        
        # Size factor (larger inversions more confident)
        size_factor = min(1.0, len(genes) / 5)
        factors.append(size_factor)
        
        # Similarity consistency
        if hasattr(genes[0], 'similarity'):
            similarities = [g['similarity'] for g in genes]
            sim_std = np.std(similarities)
            sim_factor = max(0.0, 1.0 - sim_std)
            factors.append(sim_factor)
        
        return np.mean(factors) if factors else 0.5
    else:
        return 1.0


def detect_single_gene_inversions(ortholog_df, config):
    """Detect single gene inversions"""
    single_inversions = []
    
    for _, gene in ortholog_df.iterrows():
        if gene['first_strand'] != gene['second_strand']:
            confidence = gene.get('confidence', 0.5)
            
            if confidence >= config.get('inversion_confidence_threshold', 0.7):
                single_inversions.append({
                    'first_chr': gene['first_chr'],
                    'second_chr': gene['second_chr'],
                    'start_gene': gene['busco_id'],
                    'end_gene': gene['busco_id'],
                    'size_genes': 1,
                    'inversion_type': 'single_gene_inversion',
                    'strand_pattern': 'single_flip',
                    'confidence': confidence,
                    'synteny_context': 'single_gene'
                })
    
    return single_inversions