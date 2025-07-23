"""
Quality assessment module for the Genome Inversion Analyzer
Evaluates assembly quality and suggests parameter adjustments
"""

import numpy as np
import logging
from Bio import SeqIO

logger = logging.getLogger(__name__)


def assess_assembly_quality(fasta_path, busco_df, config):
    """Comprehensive assembly quality assessment"""
    if not config.get('enable_assembly_quality_assessment', False):
        return {'quality_score': 1.0, 'quality_class': 'medium', 'adjustments': {}}
    
    logger.info(f"Assessing assembly quality for {fasta_path}")
    
    quality_metrics = {}
    
    # Basic assembly statistics
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    total_length = sum(len(seq) for seq in sequences)
    n_contigs = len(sequences)
    contig_lengths = [len(seq) for seq in sequences]
    
    if contig_lengths:
        # Calculate assembly metrics
        sorted_lengths = sorted(contig_lengths, reverse=True)
        cumsum = np.cumsum(sorted_lengths)
        
        # N50
        n50_idx = np.where(cumsum >= total_length * 0.5)[0]
        n50 = sorted_lengths[n50_idx[0]] if len(n50_idx) > 0 else 0
        
        # N90
        n90_idx = np.where(cumsum >= total_length * 0.9)[0]
        n90 = sorted_lengths[n90_idx[0]] if len(n90_idx) > 0 else 0
        
        quality_metrics.update({
            'total_length': total_length,
            'n_contigs': n_contigs,
            'n50': n50,
            'n90': n90,
            'max_contig': max(contig_lengths),
            'mean_contig_length': np.mean(contig_lengths),
            'median_contig_length': np.median(contig_lengths),
            'contig_length_std': np.std(contig_lengths)
        })
    
    # BUSCO completeness analysis
    if len(busco_df) > 0:
        complete_buscos = len(busco_df[busco_df['status'] == 'Complete'])
        fragmented_buscos = len(busco_df[busco_df['status'] == 'Fragmented'])
        missing_buscos = len(busco_df[busco_df['status'] == 'Missing'])
        duplicated_buscos = len(busco_df[busco_df['status'] == 'Duplicated'])
        total_buscos = len(busco_df)
        
        completeness = complete_buscos / total_buscos if total_buscos > 0 else 0
        fragmentation = fragmented_buscos / total_buscos if total_buscos > 0 else 0
        duplication = duplicated_buscos / total_buscos if total_buscos > 0 else 0
        
        quality_metrics.update({
            'busco_completeness': completeness,
            'busco_fragmentation': fragmentation,
            'busco_duplication': duplication,
            'busco_missing': missing_buscos / total_buscos if total_buscos > 0 else 0
        })
    
    # Calculate overall quality score and classification
    quality_score = calculate_quality_score(quality_metrics, config)
    quality_class = classify_assembly_quality(quality_score, quality_metrics, config)
    
    # Suggest parameter adjustments
    adjustments = suggest_parameter_adjustments(quality_class, quality_metrics, config)
    
    logger.info(f"  Assembly quality: {quality_class} (score: {quality_score:.3f})")
    logger.info(f"  Suggested {len(adjustments)} parameter adjustments")
    
    return {
        'quality_score': quality_score,
        'quality_class': quality_class,
        'metrics': quality_metrics,
        'adjustments': adjustments
    }


def calculate_quality_score(metrics, config):
    """Calculate comprehensive assembly quality score"""
    score_components = []
    
    # Define default thresholds if not in config
    high_quality_n50_threshold = config.get('high_quality_n50_threshold', 10000000)
    medium_quality_n50_threshold = config.get('medium_quality_n50_threshold', 1000000)
    high_quality_busco_threshold = config.get('high_quality_busco_threshold', 0.95)
    medium_quality_busco_threshold = config.get('medium_quality_busco_threshold', 0.85)
    low_quality_busco_threshold = config.get('low_quality_busco_threshold', 0.70)
    
    # N50 component
    if 'n50' in metrics:
        if metrics['n50'] >= high_quality_n50_threshold:
            score_components.append(1.0)
        elif metrics['n50'] >= medium_quality_n50_threshold:
            score_components.append(0.7)
        else:
            score_components.append(0.3)
    
    # BUSCO completeness component
    if 'busco_completeness' in metrics:
        completeness = metrics['busco_completeness']
        if completeness >= high_quality_busco_threshold:
            score_components.append(1.0)
        elif completeness >= medium_quality_busco_threshold:
            score_components.append(0.7)
        elif completeness >= low_quality_busco_threshold:
            score_components.append(0.4)
        else:
            score_components.append(0.1)
    
    # Contig count penalty (fewer contigs is better for chromosomal assemblies)
    if 'n_contigs' in metrics:
        if metrics['n_contigs'] <= 50:  # Good chromosomal assembly
            score_components.append(1.0)
        elif metrics['n_contigs'] <= 1000:  # Reasonable scaffold assembly
            score_components.append(0.6)
        else:  # Fragmented assembly
            score_components.append(0.2)
    
    # Fragmentation penalty
    if 'busco_fragmentation' in metrics:
        fragmentation_score = max(0.0, 1.0 - metrics['busco_fragmentation'] * 2)
        score_components.append(fragmentation_score)
    
    return np.mean(score_components) if score_components else 0.5


def classify_assembly_quality(quality_score, metrics, config):
    """Classify assembly quality into categories"""
    if quality_score >= 0.8:
        return 'high'
    elif quality_score >= 0.6:
        return 'medium'
    elif quality_score >= 0.3:
        return 'low'
    else:
        return 'fragmented'


def suggest_parameter_adjustments(quality_class, metrics, config):
    """Suggest parameter adjustments based on assembly quality"""
    adjustments = {}
    
    # Define default thresholds if not in config
    high_quality_similarity_threshold = config.get('high_quality_similarity_threshold', 0.8)
    medium_quality_similarity_threshold = config.get('medium_quality_similarity_threshold', 0.6)
    low_quality_similarity_threshold = config.get('low_quality_similarity_threshold', 0.3)
    fragmented_assembly_similarity_threshold = config.get('fragmented_assembly_similarity_threshold', 0.2)
    
    high_quality_min_length = config.get('high_quality_min_length', 200)
    medium_quality_min_length = config.get('medium_quality_min_length', 150)
    low_quality_min_length = config.get('low_quality_min_length', 100)
    fragmented_assembly_min_length = config.get('fragmented_assembly_min_length', 50)
    
    base_min_synteny_block_size = config.get('base_min_synteny_block_size', 3)
    strict_correlation_threshold = config.get('strict_correlation_threshold', 0.9)
    base_synteny_correlation_threshold = config.get('base_synteny_correlation_threshold', 0.8)
    relaxed_correlation_threshold = config.get('relaxed_correlation_threshold', 0.6)
    
    base_max_gap_in_synteny = config.get('base_max_gap_in_synteny', 1000000)
    adaptive_max_gap_multiplier = config.get('adaptive_max_gap_multiplier', 2.0)
    sparse_genome_gap_factor = config.get('sparse_genome_gap_factor', 3.0)
    micro_synteny_block_size = config.get('micro_synteny_block_size', 1)
    
    if quality_class == 'high':
        adjustments.update({
            'similarity_threshold': high_quality_similarity_threshold,
            'min_busco_length': high_quality_min_length,
            'min_synteny_block_size': base_min_synteny_block_size,
            'correlation_threshold': strict_correlation_threshold,
            'max_gap_in_synteny': base_max_gap_in_synteny
        })
    elif quality_class == 'medium':
        adjustments.update({
            'similarity_threshold': medium_quality_similarity_threshold,
            'min_busco_length': medium_quality_min_length,
            'min_synteny_block_size': base_min_synteny_block_size,
            'correlation_threshold': base_synteny_correlation_threshold,
            'max_gap_in_synteny': base_max_gap_in_synteny
        })
    elif quality_class == 'low':
        adjustments.update({
            'similarity_threshold': low_quality_similarity_threshold,
            'min_busco_length': low_quality_min_length,
            'min_synteny_block_size': max(1, base_min_synteny_block_size - 1),
            'correlation_threshold': relaxed_correlation_threshold,
            'max_gap_in_synteny': int(base_max_gap_in_synteny * adaptive_max_gap_multiplier)
        })
    else:  # fragmented
        adjustments.update({
            'similarity_threshold': fragmented_assembly_similarity_threshold,
            'min_busco_length': fragmented_assembly_min_length,
            'min_synteny_block_size': micro_synteny_block_size,
            'correlation_threshold': relaxed_correlation_threshold,
            'max_gap_in_synteny': int(base_max_gap_in_synteny * sparse_genome_gap_factor)
        })
    
    return adjustments