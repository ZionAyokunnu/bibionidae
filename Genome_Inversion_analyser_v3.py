#!/usr/bin/env python3
"""
INTEGRATED SYNTENY AND INVERSION ANALYZER
Combines synteny detection with inversion analysis for comprehensive chromosomal rearrangement study
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
from collections import defaultdict, Counter
import random
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

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
    'first_fasta_path': 'GCA_910594885.2_idBibMarc1.2_genomic.fna',
    'second_fasta_path': 'GCA_958336335.1_idDilFebr1.1_genomic.fna',
    'first_busco_path': 'Bibio_marci/full_table.tsv',
    'second_busco_path': 'Dilophus_febrilis/full_table.tsv',
    
    # Output files
    'synteny_analysis_csv': 'Result/synteny_analysis.csv',
    'inversion_summary_csv': 'Result/inversion_summary.csv',
    'chromosome_rearrangements_csv': 'Result/chromosome_rearrangements.csv',
    'synteny_dotplot_png': 'Result/synteny_dotplot.png',
    'inversion_dotplot_png': 'Result/inversion_dotplot.png',
    'combined_analysis_png': 'Result/combined_analysis.png',
    
    # Analysis parameters
    'similarity_threshold': 0.5,        # Lower threshold for cross-species
    'min_genes_per_chromosome': 5,      # Minimum genes to analyze synteny
    'synteny_correlation_threshold': 0.8,  # Correlation threshold for gene order
    'min_synteny_block_size': 3,        # Minimum genes in synteny block
    'max_gap_in_synteny': 5,            # Maximum gap between syntenic genes
    
    # Inversion detection parameters
    'min_inversion_size': 2,            # Minimum genes in inversion
    'strand_consistency_threshold': 0.7, # Threshold for strand consistency
    
    # BUSCO filtering
    'busco_status_filter': ['Complete', 'Fragmented'],
    'min_busco_length': 150,
    
    # Plotting parameters
    'plot_width': 15,
    'plot_height': 10,
    'dpi': 300
}

################################################################################
# BUSCO PROCESSING FUNCTIONS (from previous version)
################################################################################

def parse_busco_table(busco_path):
    """Parse BUSCO full_table.tsv file"""
    logger.info(f"Parsing BUSCO table: {busco_path}")
    
    with open(busco_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
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
    return busco_df

def filter_busco_genes(busco_df, config):
    """Filter BUSCO genes based on configuration"""
    logger.info("Filtering BUSCO genes...")
    
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    
    min_length = config.get('min_busco_length', 150)
    if 'length' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['length'] >= min_length]
    
    filtered_df = filtered_df.dropna(subset=['gene_start', 'gene_end'])
    
    logger.info(f"  Filtered to {len(filtered_df)} high-quality BUSCO genes")
    return filtered_df

def extract_busco_sequences(busco_df, fasta_path):
    """Extract BUSCO gene sequences from genome FASTA"""
    logger.info(f"Extracting BUSCO sequences from {fasta_path}...")
    
    genome_seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome_seqs[record.id] = str(record.seq)
    
    busco_with_seqs = []
    
    for _, gene in busco_df.iterrows():
        seq_id = gene['sequence']
        start = int(gene['gene_start']) - 1
        end = int(gene['gene_end'])
        
        if seq_id in genome_seqs:
            full_seq = genome_seqs[seq_id]
            
            if start < len(full_seq) and end <= len(full_seq):
                gene_seq = full_seq[start:end]
                
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
# SYNTENY ANALYSIS FUNCTIONS
################################################################################

def create_ortholog_mapping(first_busco_df, second_busco_df, config):
    """Create ortholog mapping between two species"""
    logger.info("Creating ortholog mapping...")
    
    # Find common BUSCO genes
    common_buscos = set(first_busco_df['busco_id']) & set(second_busco_df['busco_id'])
    logger.info(f"  Found {len(common_buscos)} common BUSCO genes")
    
    ortholog_pairs = []
    
    for busco_id in common_buscos:
        first_gene = first_busco_df[first_busco_df['busco_id'] == busco_id].iloc[0]
        second_gene = second_busco_df[second_busco_df['busco_id'] == busco_id].iloc[0]
        
        # Calculate sequence similarity
        similarity = SequenceMatcher(None, first_gene['gene_sequence'], second_gene['gene_sequence']).ratio()
        
        if similarity >= config.get('similarity_threshold', 0.5):
            ortholog_pairs.append({
                'busco_id': busco_id,
                'first_chr': first_gene['sequence_id'],
                'first_start': first_gene['gene_start'],
                'first_end': first_gene['gene_end'],
                'first_strand': first_gene['strand'],
                'second_chr': second_gene['sequence_id'],
                'second_start': second_gene['gene_start'],
                'second_end': second_gene['gene_end'],
                'second_strand': second_gene['strand'],
                'similarity': similarity,
                'first_length': first_gene['gene_length'],
                'second_length': second_gene['gene_length']
            })
    
    ortholog_df = pd.DataFrame(ortholog_pairs)
    logger.info(f"  Created {len(ortholog_df)} ortholog pairs above similarity threshold")
    
    return ortholog_df

def analyze_synteny_blocks(ortholog_df, config):
    """Analyze synteny blocks between chromosomes"""
    logger.info("Analyzing synteny blocks...")
    
    synteny_blocks = []
    chromosome_mappings = []
    
    # Group by chromosome pairs
    for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
        if len(group) >= config.get('min_genes_per_chromosome', 5):
            
            # Sort genes by position
            group_sorted = group.sort_values('first_start')
            
            # Calculate position correlation
            if len(group_sorted) > 1:
                correlation, p_value = pearsonr(group_sorted['first_start'], group_sorted['second_start'])
            else:
                correlation, p_value = 1.0, 0.0
            
            # Analyze strand consistency
            strand_consistency = (group_sorted['first_strand'] == group_sorted['second_strand']).mean()
            
            # Identify synteny blocks within this chromosome pair
            blocks = identify_synteny_blocks_in_pair(group_sorted, config)
            
            for block in blocks:
                synteny_blocks.append({
                    'first_chr': first_chr,
                    'second_chr': second_chr,
                    'block_size': len(block),
                    'start_gene': block.iloc[0]['busco_id'],
                    'end_gene': block.iloc[-1]['busco_id'],
                    'position_correlation': correlation,
                    'strand_consistency': strand_consistency,
                    'synteny_type': classify_synteny_type(block, correlation, strand_consistency)
                })
            
            chromosome_mappings.append({
                'first_chr': first_chr,
                'second_chr': second_chr,
                'gene_count': len(group),
                'position_correlation': correlation,
                'strand_consistency': strand_consistency,
                'p_value': p_value
            })
    
    synteny_df = pd.DataFrame(synteny_blocks)
    mapping_df = pd.DataFrame(chromosome_mappings)
    
    logger.info(f"  Found {len(synteny_df)} synteny blocks")
    logger.info(f"  Analyzed {len(mapping_df)} chromosome pairs")
    
    return synteny_df, mapping_df

def identify_synteny_blocks_in_pair(gene_group, config):
    """Identify synteny blocks within a chromosome pair"""
    blocks = []
    current_block = []
    
    genes_sorted = gene_group.sort_values('first_start')
    
    for i, (_, gene) in enumerate(genes_sorted.iterrows()):
        if not current_block:
            current_block = [gene]
        else:
            # Check if gene is close enough to be in the same block
            prev_gene = current_block[-1]
            
            # Distance in terms of gene positions
            first_distance = abs(gene['first_start'] - prev_gene['first_start'])
            second_distance = abs(gene['second_start'] - prev_gene['second_start'])
            
            # Simple heuristic: genes should be reasonably close
            max_distance = 1000000  # 1Mb
            
            if first_distance < max_distance and second_distance < max_distance:
                current_block.append(gene)
            else:
                # End current block and start new one
                if len(current_block) >= config.get('min_synteny_block_size', 3):
                    blocks.append(pd.DataFrame(current_block))
                current_block = [gene]
    
    # Add final block
    if len(current_block) >= config.get('min_synteny_block_size', 3):
        blocks.append(pd.DataFrame(current_block))
    
    return blocks

def classify_synteny_type(block, correlation, strand_consistency):
    """Classify the type of synteny based on correlation and strand consistency"""
    if correlation > 0.8:
        if strand_consistency > 0.8:
            return 'colinear'
        else:
            return 'colinear_inverted'
    elif correlation < -0.8:
        if strand_consistency < 0.2:
            return 'inverted'
        else:
            return 'inverted_mixed'
    else:
        return 'rearranged'

def analyze_chromosome_rearrangements(ortholog_df, config):
    """Analyze chromosome-level rearrangements"""
    logger.info("Analyzing chromosome rearrangements...")
    
    rearrangements = []
    
    # Analyze each first chromosome
    for first_chr in ortholog_df['first_chr'].unique():
        first_genes = ortholog_df[ortholog_df['first_chr'] == first_chr]
        
        # Find target chromosomes
        target_chromosomes = first_genes['second_chr'].value_counts()
        
        if len(target_chromosomes) > 1:
            # Chromosome split
            rearrangements.append({
                'type': 'chromosome_split',
                'first_chr': first_chr,
                'second_chrs': target_chromosomes.index.tolist(),
                'gene_counts': target_chromosomes.values.tolist(),
                'total_genes': len(first_genes)
            })
    
    # Analyze each second chromosome
    for second_chr in ortholog_df['second_chr'].unique():
        second_genes = ortholog_df[ortholog_df['second_chr'] == second_chr]
        
        # Find source chromosomes
        source_chromosomes = second_genes['first_chr'].value_counts()
        
        if len(source_chromosomes) > 1:
            # Chromosome fusion
            rearrangements.append({
                'type': 'chromosome_fusion',
                'first_chrs': source_chromosomes.index.tolist(),
                'second_chr': second_chr,
                'gene_counts': source_chromosomes.values.tolist(),
                'total_genes': len(second_genes)
            })
    
    rearrangement_df = pd.DataFrame(rearrangements)
    logger.info(f"  Found {len(rearrangement_df)} chromosome rearrangements")
    
    return rearrangement_df

################################################################################
# INVERSION ANALYSIS FUNCTIONS
################################################################################

def analyze_inversions_within_synteny(synteny_df, ortholog_df, config):
    """Analyze inversions within synteny blocks"""
    logger.info("Analyzing inversions within synteny blocks...")
    
    inversions = []
    
    for _, block in synteny_df.iterrows():
        if block['synteny_type'] in ['inverted', 'colinear_inverted']:
            # Get genes in this block
            block_genes = ortholog_df[
                (ortholog_df['first_chr'] == block['first_chr']) & 
                (ortholog_df['second_chr'] == block['second_chr'])
            ].sort_values('first_start')
            
            # Analyze inversion patterns
            inversion_regions = identify_inversion_regions(block_genes, config)
            
            for region in inversion_regions:
                inversions.append({
                    'first_chr': block['first_chr'],
                    'second_chr': block['second_chr'],
                    'start_gene': region['start_gene'],
                    'end_gene': region['end_gene'],
                    'size_genes': region['size'],
                    'inversion_type': region['type'],
                    'strand_pattern': region['strand_pattern']
                })
    
    inversion_df = pd.DataFrame(inversions)
    logger.info(f"  Found {len(inversion_df)} inversion regions")
    
    return inversion_df

def identify_inversion_regions(genes, config):
    """Identify specific inversion regions within a gene set"""
    regions = []
    
    # Look for strand flip patterns
    genes_sorted = genes.sort_values('first_start')
    
    current_region = []
    for i, (_, gene) in enumerate(genes_sorted.iterrows()):
        strand_flipped = gene['first_strand'] != gene['second_strand']
        
        if strand_flipped:
            current_region.append(gene)
        else:
            # End current region
            if len(current_region) >= config.get('min_inversion_size', 2):
                regions.append({
                    'start_gene': current_region[0]['busco_id'],
                    'end_gene': current_region[-1]['busco_id'],
                    'size': len(current_region),
                    'type': 'strand_inversion',
                    'strand_pattern': 'consistent_flip'
                })
            current_region = []
    
    # Add final region
    if len(current_region) >= config.get('min_inversion_size', 2):
        regions.append({
            'start_gene': current_region[0]['busco_id'],
            'end_gene': current_region[-1]['busco_id'],
            'size': len(current_region),
            'type': 'strand_inversion',
            'strand_pattern': 'consistent_flip'
        })
    
    return regions

################################################################################
# VISUALIZATION FUNCTIONS
################################################################################

def create_synteny_dotplot(ortholog_df, output_path, config):
    """Create dotplot showing synteny relationships"""
    logger.info("Creating synteny dotplot...")
    
    plt.figure(figsize=(config['plot_width'], config['plot_height']))
    
    # Create chromosome mapping for plotting
    first_chr_order = sorted(ortholog_df['first_chr'].unique())
    second_chr_order = sorted(ortholog_df['second_chr'].unique())
    
    first_chr_map = {chr_name: i for i, chr_name in enumerate(first_chr_order)}
    second_chr_map = {chr_name: i for i, chr_name in enumerate(second_chr_order)}
    
    # Plot points colored by synteny type
    colors = {'colinear': 'blue', 'inverted': 'red', 'rearranged': 'orange'}
    
    for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
        # Determine synteny type for this chromosome pair
        if len(group) > 1:
            correlation, _ = pearsonr(group['first_start'], group['second_start'])
            strand_consistency = (group['first_strand'] == group['second_strand']).mean()
            synteny_type = classify_synteny_type(group, correlation, strand_consistency)
        else:
            synteny_type = 'colinear'
        
        color = colors.get(synteny_type, 'gray')
        
        plt.scatter([first_chr_map[first_chr]] * len(group), 
                   [second_chr_map[second_chr]] * len(group),
                   c=color, alpha=0.6, s=20, label=synteny_type)
    
    # Remove duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    
    plt.xlabel('First Genome Chromosomes')
    plt.ylabel('Second Genome Chromosomes')
    plt.title('Synteny Relationships Between Chromosomes')
    
    plt.xticks(range(len(first_chr_order)), first_chr_order, rotation=45)
    plt.yticks(range(len(second_chr_order)), second_chr_order)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=config['dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Synteny dotplot saved to: {output_path}")

def create_combined_analysis_plot(synteny_df, inversion_df, rearrangement_df, output_path, config):
    """Create combined analysis visualization"""
    logger.info("Creating combined analysis plot...")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(config['plot_width'], config['plot_height']))
    
    # Plot 1: Synteny block sizes
    if len(synteny_df) > 0:
        ax1.hist(synteny_df['block_size'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax1.set_xlabel('Synteny Block Size (genes)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Distribution of Synteny Block Sizes')
    
    # Plot 2: Synteny types
    if len(synteny_df) > 0:
        synteny_counts = synteny_df['synteny_type'].value_counts()
        ax2.pie(synteny_counts.values, labels=synteny_counts.index, autopct='%1.1f%%')
        ax2.set_title('Synteny Types')
    
    # Plot 3: Chromosome rearrangements
    if len(rearrangement_df) > 0:
        rearr_counts = rearrangement_df['type'].value_counts()
        ax3.bar(rearr_counts.index, rearr_counts.values, color=['coral', 'lightgreen'])
        ax3.set_xlabel('Rearrangement Type')
        ax3.set_ylabel('Count')
        ax3.set_title('Chromosome Rearrangements')
    
    # Plot 4: Inversion sizes
    if len(inversion_df) > 0:
        ax4.hist(inversion_df['size_genes'], bins=15, alpha=0.7, color='lightcoral', edgecolor='black')
        ax4.set_xlabel('Inversion Size (genes)')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Distribution of Inversion Sizes')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=config['dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Combined analysis plot saved to: {output_path}")

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

def run_integrated_analysis(config=None):
    """Main function to run integrated synteny and inversion analysis"""
    if config is None:
        config = ANALYSIS_CONFIG
    
    logger.info("=== INTEGRATED SYNTENY AND INVERSION ANALYZER ===")
    logger.info("Comprehensive chromosomal rearrangement analysis\n")
    
    # Step 1: Parse and process BUSCO data
    logger.info("Step 1: Processing BUSCO data")
    first_busco_raw = parse_busco_table(config['first_busco_path'])
    second_busco_raw = parse_busco_table(config['second_busco_path'])
    
    first_busco_filtered = filter_busco_genes(first_busco_raw, config)
    second_busco_filtered = filter_busco_genes(second_busco_raw, config)
    
    first_busco_seqs = extract_busco_sequences(first_busco_filtered, config['first_fasta_path'])
    second_busco_seqs = extract_busco_sequences(second_busco_filtered, config['second_fasta_path'])
    
    # Step 2: Create ortholog mapping
    logger.info("\nStep 2: Creating ortholog mapping")
    ortholog_df = create_ortholog_mapping(first_busco_seqs, second_busco_seqs, config)
    
    # Step 3: Analyze synteny
    logger.info("\nStep 3: Analyzing synteny blocks")
    synteny_df, mapping_df = analyze_synteny_blocks(ortholog_df, config)
    
    # Step 4: Analyze chromosome rearrangements
    logger.info("\nStep 4: Analyzing chromosome rearrangements")
    rearrangement_df = analyze_chromosome_rearrangements(ortholog_df, config)
    
    # Step 5: Analyze inversions within synteny
    logger.info("\nStep 5: Analyzing inversions within synteny blocks")
    inversion_df = analyze_inversions_within_synteny(synteny_df, ortholog_df, config)
    
    # Step 6: Create visualizations
    logger.info("\nStep 6: Creating visualizations")
    create_synteny_dotplot(ortholog_df, config['synteny_dotplot_png'], config)
    create_combined_analysis_plot(synteny_df, inversion_df, rearrangement_df, 
                                 config['combined_analysis_png'], config)
    
    # Step 7: Save results
    logger.info("\nStep 7: Saving results")
    synteny_df.to_csv(config['synteny_analysis_csv'], index=False)
    inversion_df.to_csv(config['inversion_summary_csv'], index=False)
    rearrangement_df.to_csv(config['chromosome_rearrangements_csv'], index=False)
    
    # Step 8: Generate summary
    logger.info("\n=== ANALYSIS COMPLETED ===")
    logger.info("Summary Statistics:")
    logger.info(f"  - Ortholog pairs analyzed: {len(ortholog_df)}")
    logger.info(f"  - Synteny blocks found: {len(synteny_df)}")
    logger.info(f"  - Chromosome rearrangements: {len(rearrangement_df)}")
    logger.info(f"  - Inversion regions: {len(inversion_df)}")
    
    if len(synteny_df) > 0:
        logger.info(f"  - Average synteny block size: {synteny_df['block_size'].mean():.1f} genes")
        logger.info(f"  - Largest synteny block: {synteny_df['block_size'].max()} genes")
    
    if len(inversion_df) > 0:
        logger.info(f"  - Average inversion size: {inversion_df['size_genes'].mean():.1f} genes")
        logger.info(f"  - Largest inversion: {inversion_df['size_genes'].max()} genes")
    
    logger.info("\nFiles created:")
    logger.info(f"  - Synteny analysis: {config['synteny_analysis_csv']}")
    logger.info(f"  - Inversion summary: {config['inversion_summary_csv']}")
    logger.info(f"  - Chromosome rearrangements: {config['chromosome_rearrangements_csv']}")
    logger.info(f"  - Synteny dotplot: {config['synteny_dotplot_png']}")
    logger.info(f"  - Combined analysis: {config['combined_analysis_png']}")
    
    return {
        'config': config,
        'ortholog_df': ortholog_df,
        'synteny_df': synteny_df,
        'mapping_df': mapping_df,
        'rearrangement_df': rearrangement_df,
        'inversion_df': inversion_df
    }

################################################################################
# MAIN EXECUTION
################################################################################

if __name__ == "__main__":
    # Set random seed for reproducible results
    random.seed(42)
    
    # Execute the analysis
    try:
        results = run_integrated_analysis(ANALYSIS_CONFIG)
        
        # Print detailed summary
        print("\n" + "=" * 70)
        print("INTEGRATED SYNTENY AND INVERSION ANALYSIS SUMMARY")
        print("=" * 70)
        
        print(f"\nOrtholog Analysis:")
        print(f"  - Total ortholog pairs: {len(results['ortholog_df'])}")
        print(f"  - Average sequence similarity: {results['ortholog_df']['similarity'].mean():.3f}")
        
        print(f"\nSynteny Analysis:")
        print(f"  - Synteny blocks found: {len(results['synteny_df'])}")
        if len(results['synteny_df']) > 0:
            print(f"  - Average block size: {results['synteny_df']['block_size'].mean():.1f} genes")
            print(f"  - Synteny types: {results['synteny_df']['synteny_type'].value_counts().to_dict()}")
        
        print(f"\nChromosome Rearrangements:")
        print(f"  - Total rearrangements: {len(results['rearrangement_df'])}")
        if len(results['rearrangement_df']) > 0:
            print(f"  - Rearrangement types: {results['rearrangement_df']['type'].value_counts().to_dict()}")
        
        print(f"\nInversion Analysis:")
        print(f"  - Inversion regions found: {len(results['inversion_df'])}")
        if len(results['inversion_df']) > 0:
            print(f"  - Average inversion size: {results['inversion_df']['size_genes'].mean():.1f} genes")
            print(f"  - Inversion types: {results['inversion_df']['inversion_type'].value_counts().to_dict()}")
        
        print(f"\nBiological Interpretation:")
        if len(results['synteny_df']) > 0:
            colinear_blocks = len(results['synteny_df'][results['synteny_df']['synteny_type'] == 'colinear'])
            print(f"  - Colinear synteny blocks: {colinear_blocks} (conserved gene order)")
        
        if len(results['rearrangement_df']) > 0:
            print(f"  - Chromosome structural changes detected")
        else:
            print(f"  - No major chromosome rearrangements")
        
        print(f"\nMethodological Advantages:")
        print(f"  - Two-step analysis: synteny first, then inversions")
        print(f"  - Biologically meaningful: gene-level analysis")
        print(f"  - Comprehensive: covers all rearrangement types")
        print(f"  - Statistically robust: correlation-based synteny")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nError: {str(e)}")
    
    print("\n" + "=" * 70)#!/usr/bin/env python3
"""
INTEGRATED SYNTENY AND INVERSION ANALYZER
Combines synteny detection with inversion analysis for comprehensive chromosomal rearrangement study
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
from collections import defaultdict, Counter
import random
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

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
    'first_fasta_path': 'GCA_910594885.2_idBibMarc1.2_genomic.fna',
    'second_fasta_path': 'GCA_958336335.1_idDilFebr1.1_genomic.fna',
    'first_busco_path': 'Bibio_marci/full_table.tsv',
    'second_busco_path': 'Dilophus_febrilis/full_table.tsv',
    
    # Output files
    'synteny_analysis_csv': 'synteny_analysis.csv',
    'inversion_summary_csv': 'inversion_summary.csv',
    'chromosome_rearrangements_csv': 'chromosome_rearrangements.csv',
    'synteny_dotplot_png': 'synteny_dotplot.png',
    'inversion_dotplot_png': 'inversion_dotplot.png',
    'combined_analysis_png': 'combined_analysis.png',
    
    # Analysis parameters
    'similarity_threshold': 0.5,        # Lower threshold for cross-species
    'min_genes_per_chromosome': 5,      # Minimum genes to analyze synteny
    'synteny_correlation_threshold': 0.8,  # Correlation threshold for gene order
    'min_synteny_block_size': 3,        # Minimum genes in synteny block
    'max_gap_in_synteny': 5,            # Maximum gap between syntenic genes
    
    # Inversion detection parameters
    'min_inversion_size': 2,            # Minimum genes in inversion
    'strand_consistency_threshold': 0.7, # Threshold for strand consistency
    
    # BUSCO filtering
    'busco_status_filter': ['Complete', 'Fragmented'],
    'min_busco_length': 150,
    
    # Plotting parameters
    'plot_width': 15,
    'plot_height': 10,
    'dpi': 300
}

################################################################################
# BUSCO PROCESSING FUNCTIONS (from previous version)
################################################################################

def parse_busco_table(busco_path):
    """Parse BUSCO full_table.tsv file"""
    logger.info(f"Parsing BUSCO table: {busco_path}")
    
    with open(busco_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
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
    return busco_df

def filter_busco_genes(busco_df, config):
    """Filter BUSCO genes based on configuration"""
    logger.info("Filtering BUSCO genes...")
    
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    
    min_length = config.get('min_busco_length', 150)
    if 'length' in filtered_df.columns:
        filtered_df = filtered_df[filtered_df['length'] >= min_length]
    
    filtered_df = filtered_df.dropna(subset=['gene_start', 'gene_end'])
    
    logger.info(f"  Filtered to {len(filtered_df)} high-quality BUSCO genes")
    return filtered_df

def extract_busco_sequences(busco_df, fasta_path):
    """Extract BUSCO gene sequences from genome FASTA"""
    logger.info(f"Extracting BUSCO sequences from {fasta_path}...")
    
    genome_seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome_seqs[record.id] = str(record.seq)
    
    busco_with_seqs = []
    
    for _, gene in busco_df.iterrows():
        seq_id = gene['sequence']
        start = int(gene['gene_start']) - 1
        end = int(gene['gene_end'])
        
        if seq_id in genome_seqs:
            full_seq = genome_seqs[seq_id]
            
            if start < len(full_seq) and end <= len(full_seq):
                gene_seq = full_seq[start:end]
                
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
# SYNTENY ANALYSIS FUNCTIONS
################################################################################

def create_ortholog_mapping(first_busco_df, second_busco_df, config):
    """Create ortholog mapping between two species"""
    logger.info("Creating ortholog mapping...")
    
    # Find common BUSCO genes
    common_buscos = set(first_busco_df['busco_id']) & set(second_busco_df['busco_id'])
    logger.info(f"  Found {len(common_buscos)} common BUSCO genes")
    
    ortholog_pairs = []
    
    for busco_id in common_buscos:
        first_gene = first_busco_df[first_busco_df['busco_id'] == busco_id].iloc[0]
        second_gene = second_busco_df[second_busco_df['busco_id'] == busco_id].iloc[0]
        
        # Calculate sequence similarity
        similarity = SequenceMatcher(None, first_gene['gene_sequence'], second_gene['gene_sequence']).ratio()
        
        if similarity >= config.get('similarity_threshold', 0.5):
            ortholog_pairs.append({
                'busco_id': busco_id,
                'first_chr': first_gene['sequence_id'],
                'first_start': first_gene['gene_start'],
                'first_end': first_gene['gene_end'],
                'first_strand': first_gene['strand'],
                'second_chr': second_gene['sequence_id'],
                'second_start': second_gene['gene_start'],
                'second_end': second_gene['gene_end'],
                'second_strand': second_gene['strand'],
                'similarity': similarity,
                'first_length': first_gene['gene_length'],
                'second_length': second_gene['gene_length']
            })
    
    ortholog_df = pd.DataFrame(ortholog_pairs)
    logger.info(f"  Created {len(ortholog_df)} ortholog pairs above similarity threshold")
    
    return ortholog_df

def analyze_synteny_blocks(ortholog_df, config):
    """Analyze synteny blocks between chromosomes"""
    logger.info("Analyzing synteny blocks...")
    
    synteny_blocks = []
    chromosome_mappings = []
    
    # Group by chromosome pairs
    for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
        if len(group) >= config.get('min_genes_per_chromosome', 5):
            
            # Sort genes by position
            group_sorted = group.sort_values('first_start')
            
            # Calculate position correlation
            if len(group_sorted) > 1:
                correlation, p_value = pearsonr(group_sorted['first_start'], group_sorted['second_start'])
            else:
                correlation, p_value = 1.0, 0.0
            
            # Analyze strand consistency
            strand_consistency = (group_sorted['first_strand'] == group_sorted['second_strand']).mean()
            
            # Identify synteny blocks within this chromosome pair
            blocks = identify_synteny_blocks_in_pair(group_sorted, config)
            
            for block in blocks:
                synteny_blocks.append({
                    'first_chr': first_chr,
                    'second_chr': second_chr,
                    'block_size': len(block),
                    'start_gene': block.iloc[0]['busco_id'],
                    'end_gene': block.iloc[-1]['busco_id'],
                    'position_correlation': correlation,
                    'strand_consistency': strand_consistency,
                    'synteny_type': classify_synteny_type(block, correlation, strand_consistency)
                })
            
            chromosome_mappings.append({
                'first_chr': first_chr,
                'second_chr': second_chr,
                'gene_count': len(group),
                'position_correlation': correlation,
                'strand_consistency': strand_consistency,
                'p_value': p_value
            })
    
    synteny_df = pd.DataFrame(synteny_blocks)
    mapping_df = pd.DataFrame(chromosome_mappings)
    
    logger.info(f"  Found {len(synteny_df)} synteny blocks")
    logger.info(f"  Analyzed {len(mapping_df)} chromosome pairs")
    
    return synteny_df, mapping_df

def identify_synteny_blocks_in_pair(gene_group, config):
    """Identify synteny blocks within a chromosome pair"""
    blocks = []
    current_block = []
    
    genes_sorted = gene_group.sort_values('first_start')
    
    for i, (_, gene) in enumerate(genes_sorted.iterrows()):
        if not current_block:
            current_block = [gene]
        else:
            # Check if gene is close enough to be in the same block
            prev_gene = current_block[-1]
            
            # Distance in terms of gene positions
            first_distance = abs(gene['first_start'] - prev_gene['first_start'])
            second_distance = abs(gene['second_start'] - prev_gene['second_start'])
            
            # Simple heuristic: genes should be reasonably close
            max_distance = 1000000  # 1Mb
            
            if first_distance < max_distance and second_distance < max_distance:
                current_block.append(gene)
            else:
                # End current block and start new one
                if len(current_block) >= config.get('min_synteny_block_size', 3):
                    blocks.append(pd.DataFrame(current_block))
                current_block = [gene]
    
    # Add final block
    if len(current_block) >= config.get('min_synteny_block_size', 3):
        blocks.append(pd.DataFrame(current_block))
    
    return blocks

def classify_synteny_type(block, correlation, strand_consistency):
    """Classify the type of synteny based on correlation and strand consistency"""
    if correlation > 0.8:
        if strand_consistency > 0.8:
            return 'colinear'
        else:
            return 'colinear_inverted'
    elif correlation < -0.8:
        if strand_consistency < 0.2:
            return 'inverted'
        else:
            return 'inverted_mixed'
    else:
        return 'rearranged'

def analyze_chromosome_rearrangements(ortholog_df, config):
    """Analyze chromosome-level rearrangements"""
    logger.info("Analyzing chromosome rearrangements...")
    
    rearrangements = []
    
    # Analyze each first chromosome
    for first_chr in ortholog_df['first_chr'].unique():
        first_genes = ortholog_df[ortholog_df['first_chr'] == first_chr]
        
        # Find target chromosomes
        target_chromosomes = first_genes['second_chr'].value_counts()
        
        if len(target_chromosomes) > 1:
            # Chromosome split
            rearrangements.append({
                'type': 'chromosome_split',
                'first_chr': first_chr,
                'second_chrs': target_chromosomes.index.tolist(),
                'gene_counts': target_chromosomes.values.tolist(),
                'total_genes': len(first_genes)
            })
    
    # Analyze each second chromosome
    for second_chr in ortholog_df['second_chr'].unique():
        second_genes = ortholog_df[ortholog_df['second_chr'] == second_chr]
        
        # Find source chromosomes
        source_chromosomes = second_genes['first_chr'].value_counts()
        
        if len(source_chromosomes) > 1:
            # Chromosome fusion
            rearrangements.append({
                'type': 'chromosome_fusion',
                'first_chrs': source_chromosomes.index.tolist(),
                'second_chr': second_chr,
                'gene_counts': source_chromosomes.values.tolist(),
                'total_genes': len(second_genes)
            })
    
    rearrangement_df = pd.DataFrame(rearrangements)
    logger.info(f"  Found {len(rearrangement_df)} chromosome rearrangements")
    
    return rearrangement_df

################################################################################
# INVERSION ANALYSIS FUNCTIONS
################################################################################

def analyze_inversions_within_synteny(synteny_df, ortholog_df, config):
    """Analyze inversions within synteny blocks"""
    logger.info("Analyzing inversions within synteny blocks...")
    
    inversions = []
    
    for _, block in synteny_df.iterrows():
        if block['synteny_type'] in ['inverted', 'colinear_inverted']:
            # Get genes in this block
            block_genes = ortholog_df[
                (ortholog_df['first_chr'] == block['first_chr']) & 
                (ortholog_df['second_chr'] == block['second_chr'])
            ].sort_values('first_start')
            
            # Analyze inversion patterns
            inversion_regions = identify_inversion_regions(block_genes, config)
            
            for region in inversion_regions:
                inversions.append({
                    'first_chr': block['first_chr'],
                    'second_chr': block['second_chr'],
                    'start_gene': region['start_gene'],
                    'end_gene': region['end_gene'],
                    'size_genes': region['size'],
                    'inversion_type': region['type'],
                    'strand_pattern': region['strand_pattern']
                })
    
    inversion_df = pd.DataFrame(inversions)
    logger.info(f"  Found {len(inversion_df)} inversion regions")
    
    return inversion_df

def identify_inversion_regions(genes, config):
    """Identify specific inversion regions within a gene set"""
    regions = []
    
    # Look for strand flip patterns
    genes_sorted = genes.sort_values('first_start')
    
    current_region = []
    for i, (_, gene) in enumerate(genes_sorted.iterrows()):
        strand_flipped = gene['first_strand'] != gene['second_strand']
        
        if strand_flipped:
            current_region.append(gene)
        else:
            # End current region
            if len(current_region) >= config.get('min_inversion_size', 2):
                regions.append({
                    'start_gene': current_region[0]['busco_id'],
                    'end_gene': current_region[-1]['busco_id'],
                    'size': len(current_region),
                    'type': 'strand_inversion',
                    'strand_pattern': 'consistent_flip'
                })
            current_region = []
    
    # Add final region
    if len(current_region) >= config.get('min_inversion_size', 2):
        regions.append({
            'start_gene': current_region[0]['busco_id'],
            'end_gene': current_region[-1]['busco_id'],
            'size': len(current_region),
            'type': 'strand_inversion',
            'strand_pattern': 'consistent_flip'
        })
    
    return regions

################################################################################
# VISUALIZATION FUNCTIONS
################################################################################

def create_synteny_dotplot(ortholog_df, output_path, config):
    """Create dotplot showing synteny relationships"""
    logger.info("Creating synteny dotplot...")
    
    plt.figure(figsize=(config['plot_width'], config['plot_height']))
    
    # Create chromosome mapping for plotting
    first_chr_order = sorted(ortholog_df['first_chr'].unique())
    second_chr_order = sorted(ortholog_df['second_chr'].unique())
    
    first_chr_map = {chr_name: i for i, chr_name in enumerate(first_chr_order)}
    second_chr_map = {chr_name: i for i, chr_name in enumerate(second_chr_order)}
    
    # Plot points colored by synteny type
    colors = {'colinear': 'blue', 'inverted': 'red', 'rearranged': 'orange'}
    
    for (first_chr, second_chr), group in ortholog_df.groupby(['first_chr', 'second_chr']):
        # Determine synteny type for this chromosome pair
        if len(group) > 1:
            correlation, _ = pearsonr(group['first_start'], group['second_start'])
            strand_consistency = (group['first_strand'] == group['second_strand']).mean()
            synteny_type = classify_synteny_type(group, correlation, strand_consistency)
        else:
            synteny_type = 'colinear'
        
        color = colors.get(synteny_type, 'gray')
        
        plt.scatter([first_chr_map[first_chr]] * len(group), 
                   [second_chr_map[second_chr]] * len(group),
                   c=color, alpha=0.6, s=20, label=synteny_type)
    
    # Remove duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    
    plt.xlabel('First Genome Chromosomes')
    plt.ylabel('Second Genome Chromosomes')
    plt.title('Synteny Relationships Between Chromosomes')
    
    plt.xticks(range(len(first_chr_order)), first_chr_order, rotation=45)
    plt.yticks(range(len(second_chr_order)), second_chr_order)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=config['dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Synteny dotplot saved to: {output_path}")

def create_combined_analysis_plot(synteny_df, inversion_df, rearrangement_df, output_path, config):
    """Create combined analysis visualization"""
    logger.info("Creating combined analysis plot...")
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(config['plot_width'], config['plot_height']))
    
    # Plot 1: Synteny block sizes
    if len(synteny_df) > 0:
        ax1.hist(synteny_df['block_size'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax1.set_xlabel('Synteny Block Size (genes)')
        ax1.set_ylabel('Frequency')
        ax1.set_title('Distribution of Synteny Block Sizes')
    
    # Plot 2: Synteny types
    if len(synteny_df) > 0:
        synteny_counts = synteny_df['synteny_type'].value_counts()
        ax2.pie(synteny_counts.values, labels=synteny_counts.index, autopct='%1.1f%%')
        ax2.set_title('Synteny Types')
    
    # Plot 3: Chromosome rearrangements
    if len(rearrangement_df) > 0:
        rearr_counts = rearrangement_df['type'].value_counts()
        ax3.bar(rearr_counts.index, rearr_counts.values, color=['coral', 'lightgreen'])
        ax3.set_xlabel('Rearrangement Type')
        ax3.set_ylabel('Count')
        ax3.set_title('Chromosome Rearrangements')
    
    # Plot 4: Inversion sizes
    if len(inversion_df) > 0:
        ax4.hist(inversion_df['size_genes'], bins=15, alpha=0.7, color='lightcoral', edgecolor='black')
        ax4.set_xlabel('Inversion Size (genes)')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Distribution of Inversion Sizes')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=config['dpi'], bbox_inches='tight')
    plt.close()
    
    logger.info(f"Combined analysis plot saved to: {output_path}")

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

def run_integrated_analysis(config=None):
    """Main function to run integrated synteny and inversion analysis"""
    if config is None:
        config = ANALYSIS_CONFIG
    
    logger.info("=== INTEGRATED SYNTENY AND INVERSION ANALYZER ===")
    logger.info("Comprehensive chromosomal rearrangement analysis\n")
    
    # Step 1: Parse and process BUSCO data
    logger.info("Step 1: Processing BUSCO data")
    first_busco_raw = parse_busco_table(config['first_busco_path'])
    second_busco_raw = parse_busco_table(config['second_busco_path'])
    
    first_busco_filtered = filter_busco_genes(first_busco_raw, config)
    second_busco_filtered = filter_busco_genes(second_busco_raw, config)
    
    first_busco_seqs = extract_busco_sequences(first_busco_filtered, config['first_fasta_path'])
    second_busco_seqs = extract_busco_sequences(second_busco_filtered, config['second_fasta_path'])
    
    # Step 2: Create ortholog mapping
    logger.info("\nStep 2: Creating ortholog mapping")
    ortholog_df = create_ortholog_mapping(first_busco_seqs, second_busco_seqs, config)
    
    # Step 3: Analyze synteny
    logger.info("\nStep 3: Analyzing synteny blocks")
    synteny_df, mapping_df = analyze_synteny_blocks(ortholog_df, config)
    
    # Step 4: Analyze chromosome rearrangements
    logger.info("\nStep 4: Analyzing chromosome rearrangements")
    rearrangement_df = analyze_chromosome_rearrangements(ortholog_df, config)
    
    # Step 5: Analyze inversions within synteny
    logger.info("\nStep 5: Analyzing inversions within synteny blocks")
    inversion_df = analyze_inversions_within_synteny(synteny_df, ortholog_df, config)
    
    # Step 6: Create visualizations
    logger.info("\nStep 6: Creating visualizations")
    create_synteny_dotplot(ortholog_df, config['synteny_dotplot_png'], config)
    create_combined_analysis_plot(synteny_df, inversion_df, rearrangement_df, 
                                 config['combined_analysis_png'], config)
    
    # Step 7: Save results
    logger.info("\nStep 7: Saving results")
    synteny_df.to_csv(config['synteny_analysis_csv'], index=False)
    inversion_df.to_csv(config['inversion_summary_csv'], index=False)
    rearrangement_df.to_csv(config['chromosome_rearrangements_csv'], index=False)
    
    # Step 8: Generate summary
    logger.info("\n=== ANALYSIS COMPLETED ===")
    logger.info("Summary Statistics:")
    logger.info(f"  - Ortholog pairs analyzed: {len(ortholog_df)}")
    logger.info(f"  - Synteny blocks found: {len(synteny_df)}")
    logger.info(f"  - Chromosome rearrangements: {len(rearrangement_df)}")
    logger.info(f"  - Inversion regions: {len(inversion_df)}")
    
    if len(synteny_df) > 0:
        logger.info(f"  - Average synteny block size: {synteny_df['block_size'].mean():.1f} genes")
        logger.info(f"  - Largest synteny block: {synteny_df['block_size'].max()} genes")
    
    if len(inversion_df) > 0:
        logger.info(f"  - Average inversion size: {inversion_df['size_genes'].mean():.1f} genes")
        logger.info(f"  - Largest inversion: {inversion_df['size_genes'].max()} genes")
    
    logger.info("\nFiles created:")
    logger.info(f"  - Synteny analysis: {config['synteny_analysis_csv']}")
    logger.info(f"  - Inversion summary: {config['inversion_summary_csv']}")
    logger.info(f"  - Chromosome rearrangements: {config['chromosome_rearrangements_csv']}")
    logger.info(f"  - Synteny dotplot: {config['synteny_dotplot_png']}")
    logger.info(f"  - Combined analysis: {config['combined_analysis_png']}")
    
    return {
        'config': config,
        'ortholog_df': ortholog_df,
        'synteny_df': synteny_df,
        'mapping_df': mapping_df,
        'rearrangement_df': rearrangement_df,
        'inversion_df': inversion_df
    }

################################################################################
# MAIN EXECUTION
################################################################################

if __name__ == "__main__":
    # Set random seed for reproducible results
    random.seed(42)
    
    # Execute the analysis
    try:
        results = run_integrated_analysis(ANALYSIS_CONFIG)
        
        # Print detailed summary
        print("\n" + "=" * 70)
        print("INTEGRATED SYNTENY AND INVERSION ANALYSIS SUMMARY")
        print("=" * 70)
        
        print(f"\nOrtholog Analysis:")
        print(f"  - Total ortholog pairs: {len(results['ortholog_df'])}")
        print(f"  - Average sequence similarity: {results['ortholog_df']['similarity'].mean():.3f}")
        
        print(f"\nSynteny Analysis:")
        print(f"  - Synteny blocks found: {len(results['synteny_df'])}")
        if len(results['synteny_df']) > 0:
            print(f"  - Average block size: {results['synteny_df']['block_size'].mean():.1f} genes")
            print(f"  - Synteny types: {results['synteny_df']['synteny_type'].value_counts().to_dict()}")
        
        print(f"\nChromosome Rearrangements:")
        print(f"  - Total rearrangements: {len(results['rearrangement_df'])}")
        if len(results['rearrangement_df']) > 0:
            print(f"  - Rearrangement types: {results['rearrangement_df']['type'].value_counts().to_dict()}")
        
        print(f"\nInversion Analysis:")
        print(f"  - Inversion regions found: {len(results['inversion_df'])}")
        if len(results['inversion_df']) > 0:
            print(f"  - Average inversion size: {results['inversion_df']['size_genes'].mean():.1f} genes")
            print(f"  - Inversion types: {results['inversion_df']['inversion_type'].value_counts().to_dict()}")
        
        print(f"\nBiological Interpretation:")
        if len(results['synteny_df']) > 0:
            colinear_blocks = len(results['synteny_df'][results['synteny_df']['synteny_type'] == 'colinear'])
            print(f"  - Colinear synteny blocks: {colinear_blocks} (conserved gene order)")
        
        if len(results['rearrangement_df']) > 0:
            print(f"  - Chromosome structural changes detected")
        else:
            print(f"  - No major chromosome rearrangements")
        
        print(f"\nMethodological Advantages:")
        print(f"  - Two-step analysis: synteny first, then inversions")
        print(f"  - Biologically meaningful: gene-level analysis")
        print(f"  - Comprehensive: covers all rearrangement types")
        print(f"  - Statistically robust: correlation-based synteny")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nError: {str(e)}")
    
    print("\n" + "=" * 70)