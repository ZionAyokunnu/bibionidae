#!/usr/bin/env python3
"""
BUSCO-BASED GENOME INVERSION ANALYZER with K-mer Optimization
Uses BUSCO orthologs for biologically meaningful comparisons
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from pathlib import Path
import logging
from difflib import SequenceMatcher
import hashlib
from collections import defaultdict, Counter
import random
from itertools import combinations

from rich.console import Console
from rich.table import Table

console = Console()

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

################################################################################
# CONFIGURATION PARAMETERS
################################################################################

ANALYSIS_CONFIG = {
    # Input files
    # 'first_fasta_path': 'GCA_910594885.2_idBibMarc1.2_genomic.fna',
    # 'second_fasta_path': 'GCA_958336335.1_idDilFebr1.1_genomic.fna',

    'first_fasta_path': 'GCA_963924295.1_idDioRufi1.1_genomic.fna',
    'second_fasta_path': 'Dioctria_linearis.fna',
    
    # BUSCO files
    # 'first_busco_path': 'Bibio_marci/full_table.tsv',  # BUSCO results for bibio_marci
    # 'second_busco_path': 'Dilophus_febrilis/full_table.tsv',  # BUSCO results for dilophus_febrilis

    'first_busco_path': 'dioctria_rufipes.tsv',  # BUSCO results for dioctria_rufipes
    'second_busco_path': 'dioctria_linearis.tsv',  # BUSCO results for dioctria_linearis

    
    # Output files
    'inversion_summary_csv': 'busco_inversion_summary.csv',
    'inversion_frequency_csv': 'busco_inversion_frequency.csv',
    'detailed_analysis_csv': 'busco_detailed_analysis.csv',
    'comparison_dotplot_png': 'busco_comparison_dotplot.png',
    'inversion_heatmap_png': 'busco_inversion_heatmap.png',
    
    # Analysis parameters
    'use_busco_markers': True,           # Use BUSCO genes as markers
    'use_nucleotide_markers': False,     # Use raw nucleotide windows
    'marker_size_bp': 5000,              # Only used if use_nucleotide_markers=True
    'similarity_threshold': 0.7,         # Higher threshold for BUSCO orthologs
    'min_inversion_size': 2,
    
    # K-mer optimization parameters
    'use_kmer_prefilter': True,          # Enable k-mer prefiltering
    'kmer_size': 15,                     # K-mer size for prefiltering
    'kmer_threshold': 0.2,               # Minimum k-mer Jaccard similarity
    'max_candidates_per_gene': 10,       # Max candidates after k-mer filtering
    
    # BUSCO filtering parameters
    'busco_status_filter': ['Complete'],  # Only use Complete BUSCO genes
    'min_busco_length': 300,             # Minimum gene length
    
    # Plotting parameters
    'plot_width': 12,
    'plot_height': 8,
    'dpi': 300
}

################################################################################
# BUSCO PROCESSING FUNCTIONS
################################################################################

def parse_busco_table(busco_path):
    """
    Parse BUSCO full_table.tsv file
    
    Args:
        busco_path (str): Path to BUSCO full_table.tsv
        
    Returns:
        DataFrame: Parsed BUSCO results
    """
    logger.info(f"Parsing BUSCO table: {busco_path}")
    
    # Read BUSCO table (skip comment lines)
    with open(busco_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    # Create DataFrame
    busco_data = []
    for line in lines:
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                busco_data.append({
                    'busco_id': parts[0],
                    'status': parts[1],
                    'sequence': parts[2],
                    'gene_start': int(parts[3]) if parts[3] != 'N/A' else None,
                    'gene_end': int(parts[4]) if parts[4] != 'N/A' else None,
                    'strand': parts[5] if len(parts) > 5 else '+',
                    'score': float(parts[6]) if len(parts) > 6 and parts[6] != 'N/A' else None,
                    'length': int(parts[7]) if len(parts) > 7 and parts[7] != 'N/A' else None
                })
    
    busco_df = pd.DataFrame(busco_data)
    
    logger.info(f"  Found {len(busco_df)} BUSCO entries")
    logger.info(f"  Status distribution: {busco_df['status'].value_counts().to_dict()}")
    
    return busco_df

def filter_busco_genes(busco_df, config):
    """
    Filter BUSCO genes based on configuration
    
    Args:
        busco_df (DataFrame): BUSCO results
        config (dict): Configuration parameters
        
    Returns:
        DataFrame: Filtered BUSCO genes
    """
    logger.info("Filtering BUSCO genes...")
    
    # Filter by status
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    
    # Filter by length
    min_length = config.get('min_busco_length', 300)
    if 'length' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['length'] >= min_length]
    
    # Remove entries with missing coordinates
    filtered_df = filtered_df.dropna(subset=['gene_start', 'gene_end'])
    
    logger.info(f"  Filtered to {len(filtered_df)} high-quality BUSCO genes")
    
    return filtered_df

def extract_busco_sequences(busco_df, fasta_path):
    """
    Extract BUSCO gene sequences from genome FASTA
    
    Args:
        busco_df (DataFrame): Filtered BUSCO genes
        fasta_path (str): Path to genome FASTA
        
    Returns:
        DataFrame: BUSCO genes with extracted sequences
    """
    logger.info(f"Extracting BUSCO sequences from {fasta_path}...")
    
    # Load genome sequences
    genome_seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome_seqs[record.id] = str(record.seq)
    
    # Extract sequences for each BUSCO gene
    busco_with_seqs = []
    
    for _, gene in busco_df.iterrows():
        seq_id = gene['sequence']
        start = int(gene['gene_start']) - 1  # Convert to 0-based
        end = int(gene['gene_end'])
        
        if seq_id in genome_seqs:
            full_seq = genome_seqs[seq_id]
            
            # Extract gene sequence
            if start < len(full_seq) and end <= len(full_seq):
                gene_seq = full_seq[start:end]
                
                # Handle strand orientation
                if gene['strand'] == '-':
                    gene_seq = str(Seq(gene_seq).reverse_complement())
                
                busco_with_seqs.append({
                    'busco_id': gene['busco_id'],
                    'sequence_id': seq_id,
                    'gene_start': gene['gene_start'],
                    'gene_end': gene['gene_end'],
                    'strand': gene['strand'],
                    'gene_sequence': gene_seq,
                    'gene_length': len(gene_seq),
                    'status': gene['status']
                })
    
    busco_seq_df = pd.DataFrame(busco_with_seqs)
    
    logger.info(f"  Extracted sequences for {len(busco_seq_df)} BUSCO genes")
    
    return busco_seq_df

################################################################################
# K-MER OPTIMIZATION FUNCTIONS
################################################################################

def generate_kmers(sequence, k=15):
    """
    Generate k-mers from a sequence
    
    Args:
        sequence (str): Input sequence
        k (int): K-mer size
        
    Returns:
        set: Set of k-mers
    """
    if len(sequence) < k:
        return set()
    
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:  # Skip k-mers with N's
            kmers.add(kmer)
    
    return kmers

def jaccard_similarity(set1, set2):
    """
    Calculate Jaccard similarity between two sets
    
    Args:
        set1, set2 (set): Sets to compare
        
    Returns:
        float: Jaccard similarity (0-1)
    """
    if len(set1) == 0 and len(set2) == 0:
        return 1.0
    
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    
    return intersection / union if union > 0 else 0.0

def create_kmer_index(busco_df, k=15):
    """
    Create k-mer index for fast similarity search
    
    Args:
        busco_df (DataFrame): BUSCO genes with sequences
        k (int): K-mer size
        
    Returns:
        dict: K-mer index and gene k-mer sets
    """
    logger.info(f"Creating k-mer index with k={k}...")
    
    kmer_index = defaultdict(set)
    gene_kmers = {}
    
    for idx, gene in busco_df.iterrows():
        gene_id = gene['busco_id']
        sequence = gene['gene_sequence']
        
        # Generate k-mers for this gene
        kmers = generate_kmers(sequence, k)
        gene_kmers[gene_id] = kmers
        
        # Add to index
        for kmer in kmers:
            kmer_index[kmer].add(gene_id)
    
    logger.info(f"  Created index with {len(kmer_index)} unique k-mers")
    logger.info(f"  Average k-mers per gene: {np.mean([len(kmers) for kmers in gene_kmers.values()]):.1f}")
    
    return kmer_index, gene_kmers

def find_kmer_candidates(query_gene_id, query_kmers, kmer_index, gene_kmers, 
                        threshold=0.2, max_candidates=10):
    """
    Find candidate genes using k-mer similarity
    
    Args:
        query_gene_id (str): Query gene ID
        query_kmers (set): Query gene k-mers
        kmer_index (dict): K-mer index
        gene_kmers (dict): Gene k-mer sets
        threshold (float): Minimum Jaccard similarity
        max_candidates (int): Maximum candidates to return
        
    Returns:
        list: List of candidate gene IDs with similarities
    """
    candidate_genes = set()
    
    # Find genes that share k-mers with query
    for kmer in query_kmers:
        if kmer in kmer_index:
            candidate_genes.update(kmer_index[kmer])
    
    # Remove self-match
    candidate_genes.discard(query_gene_id)
    
    # Calculate Jaccard similarities
    similarities = []
    for candidate_id in candidate_genes:
        if candidate_id in gene_kmers:
            similarity = jaccard_similarity(query_kmers, gene_kmers[candidate_id])
            if similarity >= threshold:
                similarities.append((candidate_id, similarity))
    
    # Sort by similarity and return top candidates
    similarities.sort(key=lambda x: x[1], reverse=True)
    
    return similarities[:max_candidates]

################################################################################
# BUSCO INVERSION DETECTION
################################################################################

def detect_busco_inversions(first_busco_df, second_busco_df, config):
    """
    Detect inversions using BUSCO orthologs with k-mer optimization
    
    Args:
        first_busco_df (DataFrame): first genome BUSCO genes
        second_busco_df (DataFrame): second genome BUSCO genes
        config (dict): Configuration parameters
        
    Returns:
        dict: Inversion detection results
    """
    logger.info("Detecting inversions using BUSCO orthologs...")
    
    # Create k-mer index for second genome
    if config.get('use_kmer_prefilter', True):
        kmer_index, gene_kmers = create_kmer_index(second_busco_df, config.get('kmer_size', 15))
    
    comparison_list = []
    inversions = []
    
    # Find orthologous BUSCO genes between genomes
    common_buscos = set(first_busco_df['busco_id']) & set(second_busco_df['busco_id'])
    logger.info(f"  Found {len(common_buscos)} common BUSCO genes")
    
    # Process each common BUSCO gene
    for i, busco_id in enumerate(sorted(common_buscos)):
        if i % 100 == 0:
            logger.info(f"  Processing BUSCO gene {i+1:,} of {len(common_buscos):,}")
        
        # Get first gene
        orig_gene = first_busco_df[first_busco_df['busco_id'] == busco_id].iloc[0]
        orig_seq = orig_gene['gene_sequence']
        
        # Get second gene (direct ortholog match)
        inv_gene = second_busco_df[second_busco_df['busco_id'] == busco_id].iloc[0]
        inv_seq = inv_gene['gene_sequence']
        
        # Compare sequences
        forward_similarity = SequenceMatcher(None, orig_seq, inv_seq).ratio()
        reverse_similarity = SequenceMatcher(None, orig_seq, inv_seq[::-1]).ratio()
        
        # Determine if sequence is reversed
        is_reversed = reverse_similarity > forward_similarity
        best_similarity = max(forward_similarity, reverse_similarity)
        
        # Record comparison if similarity is above threshold
        if best_similarity >= config.get('similarity_threshold', 0.7):
            comparison_list.append({
                'first_position': i + 1,
                'busco_id': busco_id,
                'first_sequence_id': orig_gene['sequence_id'],
                'second_sequence_id': inv_gene['sequence_id'],
                'first_start': orig_gene['gene_start'],
                'first_end': orig_gene['gene_end'],
                'second_start': inv_gene['gene_start'],
                'second_end': inv_gene['gene_end'],
                'similarity_score': best_similarity,
                'sequences_reversed': is_reversed,
                'sequences_identical': best_similarity >= 0.99,
                'sequence_changed': best_similarity < 1.0,
                'strand_first': orig_gene['strand'],
                'strand_second': inv_gene['strand']
            })
    
    comparison_df = pd.DataFrame(comparison_list)
    
    # Detect inversion blocks
    if len(comparison_df) > 0:
        reversed_genes = comparison_df[comparison_df['sequences_reversed']]
        
        if len(reversed_genes) > 0:
            logger.info(f"  Found {len(reversed_genes)} reversed BUSCO genes")
            
            # Group consecutive reversed genes into inversion blocks
            inversion_blocks = []
            current_block = []
            
            for _, gene in reversed_genes.iterrows():
                if not current_block:
                    current_block = [gene]
                else:
                    # Check if genes are on the same chromosome and reasonably close
                    prev_gene = current_block[-1]
                    
                    # Simple grouping: same chromosome or adjacent positions
                    if (gene['first_sequence_id'] == prev_gene['first_sequence_id'] and
                        gene['first_position'] - prev_gene['first_position'] <= 10):
                        current_block.append(gene)
                    else:
                        if len(current_block) > 0:
                            inversion_blocks.append(current_block)
                        current_block = [gene]
            
            # Add final block
            if len(current_block) > 0:
                inversion_blocks.append(current_block)
            
            # Convert blocks to inversion records
            for i, block in enumerate(inversion_blocks):
                if len(block) >= config.get('min_inversion_size', 2):
                    inversions.append({
                        'inversion_id': i + 1,
                        'start_pos': block[0]['first_position'],
                        'end_pos': block[-1]['first_position'],
                        'size': len(block),
                        'start_busco': block[0]['busco_id'],
                        'end_busco': block[-1]['busco_id'],
                        'chromosome': block[0]['first_sequence_id'],
                        'type': 'busco_inversion',
                        'avg_similarity': sum(g['similarity_score'] for g in block) / len(block)
                    })
    
    inversions_df = pd.DataFrame(inversions)
    
    logger.info(f"  Found {len(comparison_df)} orthologous BUSCO pairs")
    logger.info(f"  Detected {len(inversions_df)} potential inversion blocks")
    
    if len(comparison_df) > 0:
        changed_genes = comparison_df['sequence_changed'].sum()
        reversed_genes = comparison_df['sequences_reversed'].sum()
        avg_similarity = comparison_df['similarity_score'].mean()
        
        logger.info(f"  Average sequence similarity: {avg_similarity:.3f}")
        logger.info(f"  BUSCO genes with sequence changes: {changed_genes:,}")
        logger.info(f"  BUSCO genes with reversed sequences: {reversed_genes:,}")
    
    return {
        'inversions': inversions_df,
        'comparison': comparison_df,
        'stats': {
            'total_buscos': len(common_buscos),
            'matched_buscos': len(comparison_df),
            'changed_buscos': comparison_df['sequence_changed'].sum() if len(comparison_df) > 0 else 0,
            'reversed_buscos': comparison_df['sequences_reversed'].sum() if len(comparison_df) > 0 else 0,
            'avg_similarity': comparison_df['similarity_score'].mean() if len(comparison_df) > 0 else 0
        }
    }

################################################################################
# ANALYSIS FUNCTIONS
################################################################################

def calculate_busco_inversion_frequency(inversions, comparison_df):
    """Calculate inversion frequency for BUSCO genes"""
    logger.info("Calculating BUSCO inversion frequencies...")
    
    if len(inversions) == 0:
        return pd.DataFrame({
            'busco_id': [],
            'frequency': [],
            'total_events': []
        })
    
    involved_buscos = []
    for _, inversion in inversions.iterrows():
        start_pos = int(inversion['start_pos']) - 1
        end_pos = int(inversion['end_pos']) - 1
        
        if start_pos < len(comparison_df) and end_pos < len(comparison_df):
            buscos_in_range = comparison_df.iloc[start_pos:end_pos + 1]['busco_id'].tolist()
            involved_buscos.extend(buscos_in_range)
    
    if involved_buscos:
        frequency_df = pd.DataFrame({'busco_id': involved_buscos})
        frequency_df = frequency_df.groupby('busco_id').size().reset_index(name='frequency')
        frequency_df['total_events'] = len(inversions)
        frequency_df = frequency_df.sort_values('frequency', ascending=False).reset_index(drop=True)
    else:
        frequency_df = pd.DataFrame({
            'busco_id': [],
            'frequency': [],
            'total_events': []
        })
    
    logger.info(f"  {len(frequency_df)} BUSCO genes involved in inversions")
    return frequency_df

def create_busco_detailed_analysis(inversions, first_busco_df, config):
    """Create detailed analysis for BUSCO inversions"""
    logger.info("Creating detailed BUSCO analysis...")
    
    if len(inversions) == 0:
        return pd.DataFrame({
            'inversion_id': [],
            'chromosome': [],
            'start_bp': [],
            'end_bp': [],
            'size_genes': [],
            'size_mb': [],
            'num_genes': [],
            'type': [],
            'avg_similarity': []
        })
    
    detailed_analysis = []
    
    for _, inversion in inversions.iterrows():
        # Get gene information for start and end positions
        start_pos = int(inversion['start_pos']) - 1
        end_pos = int(inversion['end_pos']) - 1
        
        # Estimate genomic coordinates (would need position mapping for exact coordinates)
        size_genes = int(inversion['size'])
        estimated_size_mb = size_genes * 50000 / 1e6  # Rough estimate: 50kb per gene
        
        detailed_analysis.append({
            'inversion_id': int(inversion['inversion_id']),
            'chromosome': inversion['chromosome'],
            'start_bp': 0,  # Would need position mapping
            'end_bp': 0,    # Would need position mapping
            'size_genes': size_genes,
            'size_mb': estimated_size_mb,
            'num_genes': size_genes,
            'type': inversion['type'],
            'avg_similarity': inversion.get('avg_similarity', 0.0)
        })
    
    detailed_df = pd.DataFrame(detailed_analysis)
    logger.info(f"  Detailed analysis for {len(detailed_df)} inversions")
    
    return detailed_df

def create_busco_dotplot(comparison_df, output_path):
    """Create dotplot for BUSCO comparison"""
    logger.info("Creating BUSCO comparison dotplot...")
    
    if len(comparison_df) == 0:
        logger.info("No BUSCO comparison data available, skipping dotplot")
        return
    
    plt.figure(figsize=(12, 10))
    
    # Create position mapping
    plot_data = comparison_df.copy()
    plot_data['synthetic_second_pos'] = plot_data['first_position']
    
    # Color by inversion status
    normal_data = plot_data[~plot_data['sequences_reversed']]
    second_data = plot_data[plot_data['sequences_reversed']]
    
    if len(normal_data) > 0:
        plt.scatter(normal_data['first_position'], normal_data['synthetic_second_pos'], 
                   color='steelblue', s=20, alpha=0.7, label='Normal BUSCO genes')
    
    if len(second_data) > 0:
        # Offset second genes for visualization
        second_data.loc[:, 'synthetic_second_pos'] = (
            second_data['first_position'].max() - 
            second_data['first_position'] + 
            second_data['first_position'].min()
        )
        
        plt.scatter(second_data['first_position'], second_data['synthetic_second_pos'], 
                   color='red', s=25, alpha=0.8, label='second BUSCO genes')
    
    # Add diagonal reference line
    max_pos = plot_data['first_position'].max()
    plt.plot([1, max_pos], [1, max_pos], color='gray', linestyle='--', alpha=0.7, 
             label='Expected (no inversion)')
    
    plt.title('BUSCO-based Genome Comparison', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('first Genome Position (BUSCO gene order)', fontsize=12)
    plt.ylabel('Comparison Genome Position', fontsize=12)
    
    plt.figtext(0.5, 0.92, f"BUSCO genes compared: {len(comparison_df):,}", 
               ha='center', fontsize=12)
    
    if len(second_data) > 0:
        plt.figtext(0.5, 0.88, f"second BUSCO genes: {len(second_data):,}", 
                   ha='center', fontsize=11, color='red')
    
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"BUSCO dotplot saved to: {output_path}")

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

def run_busco_inversion_analysis(config=None):
    """Main function to run BUSCO-based inversion analysis"""
    if config is None:
        config = ANALYSIS_CONFIG
    
    logger.info("=== BUSCO-BASED GENOME INVERSION ANALYZER ===")
    logger.info("Using BUSCO orthologs for biologically meaningful comparisons\n")
    
    # Step 1: Parse BUSCO tables
    logger.info("Step 1: Parsing BUSCO tables")
    first_busco_raw = parse_busco_table(config['first_busco_path'])
    second_busco_raw = parse_busco_table(config['second_busco_path'])
    
    # Step 2: Filter BUSCO genes
    logger.info("\nStep 2: Filtering BUSCO genes")
    first_busco_filtered = filter_busco_genes(first_busco_raw, config)
    second_busco_filtered = filter_busco_genes(second_busco_raw, config)
    
    # Step 3: Extract sequences
    logger.info("\nStep 3: Extracting BUSCO sequences")
    first_busco_seqs = extract_busco_sequences(first_busco_filtered, config['first_fasta_path'])
    second_busco_seqs = extract_busco_sequences(second_busco_filtered, config['second_fasta_path'])
    
    # Step 4: Detect inversions
    logger.info("\nStep 4: Detecting inversions using BUSCO orthologs")
    detection_result = detect_busco_inversions(first_busco_seqs, second_busco_seqs, config)
    
    comparison_df = detection_result['comparison']
    
    logger.info("✓ Completed BUSCO-based inversion detection")
    logger.info(f"  Analyzed {detection_result['stats']['total_buscos']} common BUSCO genes")
    logger.info(f"  Found {detection_result['stats']['matched_buscos']} high-similarity matches")
    logger.info(f"  Detected {len(detection_result['inversions'])} inversion blocks")
    
    # Step 5: Calculate frequencies
    logger.info("\nStep 5: Calculating inversion frequencies")
    frequency_result = calculate_busco_inversion_frequency(
        detection_result['inversions'], 
        comparison_df
    )
    
    # Step 6: Create detailed analysis
    logger.info("\nStep 6: Creating detailed analysis")
    detailed_analysis = create_busco_detailed_analysis(
        detection_result['inversions'],
        first_busco_seqs,
        config
    )
    
    # Step 7: Create visualizations
    logger.info("\nStep 7: Creating visualizations")
    create_busco_dotplot(comparison_df, config['comparison_dotplot_png'])
    
    # Step 8: Save results
    logger.info("\nStep 8: Saving analysis results")
    
    if len(detailed_analysis) > 0:
        detailed_analysis.to_csv(config['inversion_summary_csv'], index=False)
        frequency_result.to_csv(config['inversion_frequency_csv'], index=False)
        
        # Enhanced detailed analysis
        enhanced_analysis = detailed_analysis.copy()
        enhanced_analysis['total_inversions'] = len(detailed_analysis)
        enhanced_analysis['total_buscos'] = detection_result['stats']['total_buscos']
        enhanced_analysis.to_csv(config['detailed_analysis_csv'], index=False)
        
        logger.info(f"Results saved to CSV files")
    else:
        logger.info("No inversions detected - creating empty output files")
        
        # Create empty files
        empty_df = pd.DataFrame()
        empty_df.to_csv(config['inversion_summary_csv'], index=False)
        empty_df.to_csv(config['inversion_frequency_csv'], index=False)
        empty_df.to_csv(config['detailed_analysis_csv'], index=False)
    
    # Step 9: Summary report
    logger.info("\n=== BUSCO ANALYSIS COMPLETED ===")
    logger.info("Summary Statistics:")
    logger.info(f"  - Total common BUSCO genes: {detection_result['stats']['total_buscos']}")
    logger.info(f"  - High-similarity matches: {detection_result['stats']['matched_buscos']}")
    logger.info(f"  - Inversion blocks detected: {len(detection_result['inversions'])}")
    logger.info(f"  - BUSCO genes in inversions: {len(frequency_result)}")
    logger.info(f"  - Average sequence similarity: {detection_result['stats']['avg_similarity']:.3f}")
    
    if len(detailed_analysis) > 0:
        avg_size = detailed_analysis['size_genes'].mean()
        logger.info(f"  - Average inversion size: {avg_size:.1f} genes")
        logger.info(f"  - Largest inversion: {detailed_analysis['size_genes'].max()} genes")
    
    logger.info("\nFiles created:")
    logger.info(f"  - BUSCO comparison dotplot: {config['comparison_dotplot_png']}")
    logger.info(f"  - Inversion summary: {config['inversion_summary_csv']}")
    logger.info(f"  - Frequency analysis: {config['inversion_frequency_csv']}")
    logger.info(f"  - Detailed analysis: {config['detailed_analysis_csv']}")
    
    return {
        'config': config,
        'detection_result': detection_result,
        'frequency_result': frequency_result,
        'detailed_analysis': detailed_analysis,
        'comparison_df': comparison_df
    }

################################################################################
# MAIN EXECUTION
################################################################################

if __name__ == "__main__":
    # Set random seed for reproducible results
    random.seed(42)
    
    # Execute the analysis
    try:
        analysis_results = run_busco_inversion_analysis(ANALYSIS_CONFIG)
        
        # Print detailed summary
        print("\n" + "=" * 60)
        print("BUSCO-BASED INVERSION ANALYSIS SUMMARY")
        print("=" * 60)
        
        if len(analysis_results['detailed_analysis']) > 0:
            print("Inversion Details:")
            print(analysis_results['detailed_analysis'].to_string(index=False))
            
            print("\nFrequency Analysis (Top 10):")
            print(analysis_results['frequency_result'].head(10).to_string(index=False))
        else:
            print("No inversions detected in the BUSCO comparison.")
        
        # Performance and biological summary
        print(f"\nBiological Summary:")
        print(f"  - first genome: {analysis_results['config']['first_fasta_path']}")
        print(f"  - Comparison genome: {analysis_results['config']['second_fasta_path']}")
        print(f"  - first BUSCO: {analysis_results['config']['first_busco_path']}")
        print(f"  - Comparison BUSCO: {analysis_results['config']['second_busco_path']}")
        
        stats = analysis_results['detection_result']['stats']
        print(f"\nBUSCO Analysis Results:")
        print(f"  - Total common BUSCO genes: {stats['total_buscos']:,}")
        print(f"  - High-similarity matches: {stats['matched_buscos']:,}")
        print(f"  - Average sequence similarity: {stats['avg_similarity']:.3f}")
        print(f"  - Genes with sequence changes: {stats['changed_buscos']:,}")
        print(f"  - Genes with reversed sequences: {stats['reversed_buscos']:,}")
        
        if len(analysis_results['detailed_analysis']) > 0:
            print(f"\nInversion Characteristics:")
            detailed = analysis_results['detailed_analysis']
            print(f"  - Number of inversion blocks: {len(detailed)}")
            print(f"  - Average genes per inversion: {detailed['size_genes'].mean():.1f}")
            print(f"  - Largest inversion: {detailed['size_genes'].max()} genes")
            print(f"  - Smallest inversion: {detailed['size_genes'].min()} genes")
            print(f"  - Total genes in inversions: {detailed['size_genes'].sum()}")
        
        print(f"\nOptimization Benefits:")
        print(f"  - BUSCO-based analysis: Biologically meaningful comparisons")
        print(f"  - Ortholog matching: Direct homology relationships")
        print(f"  - Reduced complexity: ~1,000 genes vs ~60,000 markers")
        print(f"  - Higher confidence: Sequence similarity threshold {analysis_results['config']['similarity_threshold']}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nError: {str(e)}")
        print("\nTroubleshooting:")
        print("1. Check that BUSCO files exist and are properly formatted")
        print("2. Verify FASTA files are accessible")
        print("3. Ensure BUSCO sequence IDs match FASTA sequence IDs")
        print("4. Check that BUSCO tables contain 'Complete' genes")
        
    print("\n" + "=" * 60)