#!/usr/bin/env python3
"""
MODULE 2: GENOME INVERSION ANALYZER (Python Version)
Purpose: Compare original and inverted FASTA files, detect inversions, 
         create analysis plots, and generate detailed CSV reports
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
    'original_fasta_path': 'GCA_963924295.1_idDioRufi1.1_genomic.fna',
    # 'inverted_fasta_path': 'Bibio_marci_inverted.fasta',
    # 'marker_order_file': 'Bibio_marci_inverted_marker_order.csv',  # Optional
    'inverted_fasta_path': 'Dioctria_linearis.fna',
    
    # Output files
    'inversion_summary_csv': 'inversion_summary.csv',
    'inversion_frequency_csv': 'inversion_frequency.csv',
    'detailed_analysis_csv': 'detailed_inversion_analysis.csv',
    'comparison_dotplot_png': 'comparison_dotplot.png',
    'inversion_heatmap_png': 'inversion_heatmap.png',
    
    # Analysis parameters
    'marker_size_bp': 5000,
    'detect_inversion_events': True,  # Boolean toggle for event detection
    'min_inversion_size': 2,          # Minimum markers for inversion detection
    'use_marker_order_file': False,   # Set to True if using marker order file
    
    # Plotting parameters
    'plot_width': 12,
    'plot_height': 8,
    'dpi': 300
}

################################################################################
# CORE ANALYSIS FUNCTIONS
################################################################################

def process_fasta_to_markers_with_sequences(fasta_path, marker_size=1000):
    """
    Process FASTA file into comparable marker format WITH SEQUENCE CONTENT
    
    Args:
        fasta_path (str): Path to FASTA file
        marker_size (int): Size of each marker in base pairs
        
    Returns:
        tuple: (genome_data, markers_df)
    """
    logger.info(f"Processing FASTA file with sequence content: {fasta_path}")
    
    # Read sequences
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    
    # Create genome data
    genome_data = []
    for seq_record in sequences:
        genome_data.append({
            'chromosome': seq_record.id,
            'chromosome_size_bp': len(seq_record.seq),
            'sequence': str(seq_record.seq)
        })
    
    genome_df = pd.DataFrame(genome_data)
    
    # Generate markers WITH SEQUENCE CONTENT
    markers_list = []
    
    for _, row in genome_df.iterrows():
        chromosome = row['chromosome']
        chr_size = row['chromosome_size_bp']
        full_sequence = row['sequence']
        
        marker_count = chr_size // marker_size
        
        if marker_count > 0:  # Skip very small scaffolds
            for marker_idx in range(1, marker_count + 1):
                start_bp = (marker_idx - 1) * marker_size + 1
                end_bp = marker_idx * marker_size
                marker_id = f"{chromosome}_{marker_idx:06d}"
                
                # Extract the actual sequence for this marker
                seq_start = start_bp - 1  # Convert to 0-based
                seq_end = min(end_bp, len(full_sequence))
                marker_sequence = full_sequence[seq_start:seq_end]
                
                # Create a sequence hash for comparison
                sequence_hash = hash(marker_sequence)
                
                markers_list.append({
                    'chromosome': chromosome,
                    'marker_index': marker_idx,
                    'start_bp': start_bp,
                    'end_bp': end_bp,
                    'marker_id': marker_id,
                    'size_mb': marker_size / 1e6,
                    'sequence': marker_sequence,
                    'sequence_hash': sequence_hash,
                    'sequence_length': len(marker_sequence)
                })
    
    markers_df = pd.DataFrame(markers_list)
    markers_df = markers_df.sort_values(['chromosome', 'marker_index']).reset_index(drop=True)
    markers_df['linear_position'] = range(1, len(markers_df) + 1)
    
    logger.info(f"  Generated {len(markers_df):,} markers with sequences")
    
    return {
        'genome_data': genome_df,
        'markers_df': markers_df
    }

def detect_sequence_based_inversions(original_markers_df, inverted_markers_df):
    """
    Detect inversions by comparing actual sequences, not just positions
    
    Args:
        original_markers_df (DataFrame): Original markers with sequences
        inverted_markers_df (DataFrame): Inverted markers with sequences
        
    Returns:
        dict: Dictionary with inversion detection results
    """
    logger.info("Detecting inversions using sequence-based comparison...")
    
    # Create comparison based on sequence content changes
    comparison_list = []
    inversions = []
    inversion_count = 0
    
    # Compare each marker's sequence
    for i, (_, orig_marker) in enumerate(original_markers_df.iterrows()):
        marker_id = orig_marker['marker_id']
        
        # Find corresponding marker in inverted genome
        inv_marker = inverted_markers_df[inverted_markers_df['marker_id'] == marker_id]
        
        if len(inv_marker) > 0:
            inv_marker = inv_marker.iloc[0]
            
            # Check if sequences are different
            orig_seq = orig_marker['sequence']
            inv_seq = inv_marker['sequence']
            
            sequences_identical = (orig_seq == inv_seq)
            sequences_reversed = (orig_seq == inv_seq[::-1])
            
            comparison_list.append({
                'original_position': i + 1,
                'marker_id': marker_id,
                'chromosome': orig_marker['chromosome'],
                'start_bp': orig_marker['start_bp'],
                'end_bp': orig_marker['end_bp'],
                'sequences_identical': sequences_identical,
                'sequences_reversed': sequences_reversed,
                'sequence_changed': not sequences_identical
            })
            
            # If sequences are reversed, this indicates a nucleotide-level inversion
            if sequences_reversed and not sequences_identical:
                # Check if this is part of a larger inversion block
                if (inversion_count == 0 or 
                    i - inversions[-1]['end_pos'] > 10):  # New inversion block
                    inversion_count += 1
                    inversions.append({
                        'inversion_id': inversion_count,
                        'start_pos': i + 1,
                        'end_pos': i + 1,
                        'size': 1,
                        'start_marker': marker_id,
                        'end_marker': marker_id,
                        'type': 'nucleotide_inversion'
                    })
                else:
                    # Extend existing inversion block
                    inversions[-1]['end_pos'] = i + 1
                    inversions[-1]['size'] = inversions[-1]['end_pos'] - inversions[-1]['start_pos'] + 1
                    inversions[-1]['end_marker'] = marker_id
    
    comparison_df = pd.DataFrame(comparison_list)
    inversions_df = pd.DataFrame(inversions)
    
    # Merge adjacent inversions into larger blocks
    if len(inversions_df) > 0:
        merged_inversions = []
        current_inversion = inversions_df.iloc[0].to_dict()
        
        for i in range(1, len(inversions_df)):
            next_inv = inversions_df.iloc[i]
            
            # If inversions are adjacent (within 5 markers), merge them
            if next_inv['start_pos'] - current_inversion['end_pos'] <= 5:
                current_inversion['end_pos'] = next_inv['end_pos']
                current_inversion['size'] = current_inversion['end_pos'] - current_inversion['start_pos'] + 1
                current_inversion['end_marker'] = next_inv['end_marker']
            else:
                merged_inversions.append(current_inversion)
                current_inversion = next_inv.to_dict()
        
        merged_inversions.append(current_inversion)
        inversions_df = pd.DataFrame(merged_inversions)
        
        # Renumber inversion IDs
        inversions_df['inversion_id'] = range(1, len(inversions_df) + 1)
    
    logger.info(f"  Detected {len(inversions_df)} nucleotide-level inversion blocks")
    
    # Calculate summary statistics
    total_markers = len(comparison_df)
    changed_markers = comparison_df['sequence_changed'].sum()
    reversed_markers = comparison_df['sequences_reversed'].sum()
    
    logger.info(f"  Total markers: {total_markers:,}")
    logger.info(f"  Markers with changed sequences: {changed_markers:,}")
    logger.info(f"  Markers with reversed sequences: {reversed_markers:,}")
    
    return {
        'inversions': inversions_df,
        'comparison': comparison_df,
        'stats': {
            'total_markers': total_markers,
            'changed_markers': changed_markers,
            'reversed_markers': reversed_markers
        }
    }

def detect_inversions_improved(comparison_df):
    """
    Detect inversions between two genome arrangements (improved version)
    
    Args:
        comparison_df (DataFrame): Comparison dataframe with original and inverted positions
        
    Returns:
        dict: Dictionary with inversion detection results
    """
    logger.info("Detecting inversions (improved algorithm)...")
    
    # Sort by original position to ensure proper order
    comparison_df = comparison_df.sort_values('original_position').reset_index(drop=True)
    
    # Initialize inversions list
    inversions = []
    
    if len(comparison_df) < 2:
        return {
            'inversions': pd.DataFrame(inversions),
            'comparison': comparison_df
        }
    
    # Detect inversions by finding regions where inverted_position decreases
    i = 0
    inversion_count = 0
    
    while i < len(comparison_df) - 1:
        # Check if this is the start of an inversion (decreasing sequence)
        if comparison_df.iloc[i]['inverted_position'] > comparison_df.iloc[i + 1]['inverted_position']:
            inversion_count += 1
            start_pos = i + 1  # Convert to 1-based indexing
            
            # Find the end of this inversion
            j = i + 1
            while (j < len(comparison_df) - 1 and 
                   comparison_df.iloc[j]['inverted_position'] <= comparison_df.iloc[j - 1]['inverted_position']):
                j += 1
            end_pos = j  # Convert to 1-based indexing
            
            # Record this inversion if it's large enough
            inversion_size = end_pos - start_pos + 1
            if inversion_size >= 2:  # Minimum size filter
                inversions.append({
                    'inversion_id': inversion_count,
                    'start_pos': start_pos,
                    'end_pos': end_pos,
                    'size': inversion_size,
                    'start_marker': comparison_df.iloc[start_pos - 1]['marker_id'],  # Convert back to 0-based
                    'end_marker': comparison_df.iloc[end_pos - 1]['marker_id'],
                    'type': 'inversion'
                })
            
            i = j
        else:
            i += 1
    
    inversions_df = pd.DataFrame(inversions)
    logger.info(f"  Detected {len(inversions_df)} inversion blocks")
    
    return {
        'inversions': inversions_df,
        'comparison': comparison_df
    }

def calculate_inversion_frequency_improved(inversions, comparison_df):
    """
    Count inversion events per marker (improved version)
    
    Args:
        inversions (DataFrame): Inversion detection results
        comparison_df (DataFrame): Comparison dataframe
        
    Returns:
        DataFrame: DataFrame with inversion frequencies
    """
    logger.info("Calculating inversion frequencies (improved)...")
    
    if len(inversions) == 0:
        return pd.DataFrame({
            'marker_id': [],
            'frequency': [],
            'total_events': []
        })
    
    # Get all markers involved in inversions
    involved_markers = []
    for _, inversion in inversions.iterrows():
        start_pos = int(inversion['start_pos']) - 1  # Convert to 0-based
        end_pos = int(inversion['end_pos']) - 1      # Convert to 0-based
        
        # Get markers in this range
        markers_in_range = comparison_df.iloc[start_pos:end_pos + 1]['marker_id'].tolist()
        involved_markers.extend(markers_in_range)
    
    # Count frequencies
    if involved_markers:
        frequency_df = pd.DataFrame({'marker_id': involved_markers})
        frequency_df = frequency_df.groupby('marker_id').size().reset_index(name='frequency')
        frequency_df['total_events'] = len(inversions)
        frequency_df = frequency_df.sort_values('frequency', ascending=False).reset_index(drop=True)
    else:
        frequency_df = pd.DataFrame({
            'marker_id': [],
            'frequency': [],
            'total_events': []
        })
    
    logger.info(f"  {len(frequency_df)} markers involved in inversions")
    
    return frequency_df

def create_detailed_analysis_improved(inversions, markers_df, config):
    """
    Create comprehensive inversion analysis (improved version)
    
    Args:
        inversions (DataFrame): Detected inversions
        markers_df (DataFrame): Markers dataframe
        config (dict): Analysis configuration
        
    Returns:
        DataFrame: DataFrame with detailed analysis
    """
    logger.info("Creating detailed analysis (improved)...")
    
    if len(inversions) == 0:
        return pd.DataFrame({
            'inversion_id': [],
            'chromosome': [],
            'start_bp': [],
            'end_bp': [],
            'size_markers': [],
            'size_mb': [],
            'num_genes': [],
            'inversion_rate_per_mb': [],
            'type': [],
            'inversion_events': []
        })
    
    # Create mapping from positions to genomic coordinates
    position_to_coords = markers_df.copy()
    position_to_coords['linear_position'] = range(1, len(position_to_coords) + 1)
    
    # Enhance inversions with genomic coordinates
    detailed_analysis = []
    
    for _, inversion in inversions.iterrows():
        start_pos = int(inversion['start_pos'])
        end_pos = int(inversion['end_pos'])
        
        # Get start position info
        start_info = position_to_coords[position_to_coords['linear_position'] == start_pos]
        end_info = position_to_coords[position_to_coords['linear_position'] == end_pos]
        
        if len(start_info) > 0 and len(end_info) > 0:
            start_info = start_info.iloc[0]
            end_info = end_info.iloc[0]
            
            chromosome = start_info['chromosome'] if pd.notna(start_info['chromosome']) else end_info['chromosome']
            size_markers = int(inversion['size'])
            size_mb = size_markers * config['marker_size_bp'] / 1e6
            inversion_rate_per_mb = 1 / size_mb if size_mb > 0 else 0
            
            detailed_analysis.append({
                'inversion_id': int(inversion['inversion_id']),
                'chromosome': chromosome,
                'start_bp': int(start_info['start_bp']),
                'end_bp': int(end_info['end_bp']),
                'size_markers': size_markers,
                'size_mb': size_mb,
                'num_genes': size_markers,  # Assuming 1 gene per marker
                'inversion_rate_per_mb': inversion_rate_per_mb,
                'type': inversion['type'],
                'inversion_events': 1 if config['detect_inversion_events'] else 0
            })
    
    detailed_df = pd.DataFrame(detailed_analysis)
    logger.info(f"  Detailed analysis for {len(detailed_df)} inversions")
    
    return detailed_df

def create_comparison_dotplot(comparison_df, output_path, title="Genome Comparison"):
    """
    Create comparison dotplot
    
    Args:
        comparison_df (DataFrame): Comparison dataframe
        output_path (str): Output path for plot
        title (str): Plot title
    """
    logger.info("Creating comparison dotplot...")
    
    # Handle different comparison formats
    if 'inverted_position' in comparison_df.columns:
        # Traditional position-based comparison
        x_col = 'original_position'
        y_col = 'inverted_position'
        plot_data = comparison_df[comparison_df[y_col].notna()]
    else:
        # Sequence-based comparison - create synthetic positions
        plot_data = comparison_df.copy()
        plot_data['synthetic_inverted_pos'] = plot_data['original_position']
        
        # Highlight inverted regions differently
        if 'sequences_reversed' in plot_data.columns:
            # Create offset for reversed sequences to show inversion
            reversed_mask = plot_data['sequences_reversed']
            plot_data.loc[reversed_mask, 'synthetic_inverted_pos'] = (
                plot_data.loc[reversed_mask, 'original_position'].max() - 
                plot_data.loc[reversed_mask, 'original_position'] + 
                plot_data.loc[reversed_mask, 'original_position'].min()
            )
        
        x_col = 'original_position'
        y_col = 'synthetic_inverted_pos'
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    
    if 'sequences_reversed' in plot_data.columns:
        # Color code by inversion status
        normal_data = plot_data[~plot_data['sequences_reversed']]
        inverted_data = plot_data[plot_data['sequences_reversed']]
        
        if len(normal_data) > 0:
            plt.scatter(normal_data[x_col], normal_data[y_col], 
                       color='steelblue', s=0.8, alpha=0.6, label='Normal')
        
        if len(inverted_data) > 0:
            plt.scatter(inverted_data[x_col], inverted_data[y_col], 
                       color='red', s=1.2, alpha=0.8, label='Inverted')
        
        plt.legend()
    else:
        plt.scatter(plot_data[x_col], plot_data[y_col], 
                   color='steelblue', s=0.8, alpha=0.6)
    
    # Add diagonal reference line for normal data
    if 'sequences_reversed' not in plot_data.columns or len(plot_data[~plot_data.get('sequences_reversed', False)]) > 0:
        max_val = max(plot_data[x_col].max(), plot_data[y_col].max())
        plt.plot([1, max_val], [1, max_val], color='gray', linestyle='--', alpha=0.7, label='Expected (no inversion)')
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.figtext(0.5, 0.92, f"Markers compared: {len(comparison_df):,}", 
               ha='center', fontsize=12)
    
    if 'sequences_reversed' in comparison_df.columns:
        inverted_count = comparison_df['sequences_reversed'].sum()
        plt.figtext(0.5, 0.88, f"Inverted markers: {inverted_count:,}", 
                   ha='center', fontsize=11, color='red')
    
    # plt.xlabel("Original Genome Position", fontsize=11)
    # plt.ylabel("Rearranged Genome Position", fontsize=11)

    plt.xlabel("Dioctria rufipes", fontsize=11)
    plt.ylabel("Dioctria linearis", fontsize=11)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Comparison dotplot saved to: {output_path}")

def create_inversion_heatmap(frequency_df, markers_df, output_path):
    """
    Create inversion heatmap
    
    Args:
        frequency_df (DataFrame): Frequency dataframe
        markers_df (DataFrame): Markers dataframe
        output_path (str): Output path for heatmap
    """
    logger.info("Creating inversion heatmap...")
    
    if len(frequency_df) == 0:
        logger.info("No inversions detected, skipping heatmap")
        return None
    
    # Join frequency with position data
    heatmap_data = markers_df.merge(frequency_df, on='marker_id', how='left')
    heatmap_data['frequency'] = heatmap_data['frequency'].fillna(0)
    
    # Create frequency categories
    def freq_category(freq):
        if freq == 0:
            return "0"
        elif freq == 1:
            return "1"
        elif freq == 2:
            return "2"
        else:
            return "3+"
    
    heatmap_data['freq_category'] = heatmap_data['frequency'].apply(freq_category)
    
    # Create heatmap
    plt.figure(figsize=(14, 6))
    
    # Create a pivot table for the heatmap
    chromosomes = heatmap_data['chromosome'].unique()
    colors = {'0': 'white', '1': 'lightblue', '2': 'blue', '3+': 'darkblue'}
    
    for i, chromosome in enumerate(chromosomes):
        chrom_data = heatmap_data[heatmap_data['chromosome'] == chromosome]
        
        # Create color map
        color_values = [colors[cat] for cat in chrom_data['freq_category']]
        
        plt.scatter(chrom_data['linear_position'], [i] * len(chrom_data), 
                   c=color_values, s=0.5, alpha=0.8)
    
    plt.yticks(range(len(chromosomes)), chromosomes, fontsize=8)
    plt.xlabel("Linear Genome Position", fontsize=11)
    plt.ylabel("Chromosome", fontsize=11)
    plt.title("Inversion Frequency Heatmap", fontsize=14, fontweight='bold')
    
    # Create custom legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=colors[cat], markersize=8, label=cat)
                      for cat in ['0', '1', '2', '3+']]
    plt.legend(handles=legend_elements, title='Inversion\nFrequency', loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Inversion heatmap saved to: {output_path}")

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

def run_genome_inversion_analysis(config=None):
    """
    Main function to run genome inversion analysis
    
    Args:
        config (dict): Analysis configuration
        
    Returns:
        dict: Results of the analysis
    """
    if config is None:
        config = ANALYSIS_CONFIG
    
    logger.info("=== GENOME INVERSION ANALYZER ===")
    logger.info("Analyzing inversions between original and inverted genomes\n")
    
    # Step 1: Process both FASTA files WITH SEQUENCE CONTENT
    logger.info("Step 1: Processing FASTA files with sequence analysis")
    original_result = process_fasta_to_markers_with_sequences(config['original_fasta_path'], config['marker_size_bp'])
    inverted_result = process_fasta_to_markers_with_sequences(config['inverted_fasta_path'], config['marker_size_bp'])
    
    original_markers = original_result['markers_df']['marker_id'].tolist()
    
    # Check if marker order file exists (from simulation)
    use_marker_file = (config.get('use_marker_order_file', False) and 
                      config.get('marker_order_file') and 
                      os.path.exists(config['marker_order_file']))
    
    if use_marker_file:
        logger.info(f"Using marker order file: {config['marker_order_file']}")
        
        # Read the marker order from simulation
        marker_order_df = pd.read_csv(config['marker_order_file'])
        
        # Create comparison directly from the marker order file
        comparison_df = marker_order_df[['original_position', 'marker_id', 'inverted_position']].copy()
        comparison_df = comparison_df.sort_values('original_position').reset_index(drop=True)
        
        logger.info(f"  Loaded marker order with {len(comparison_df):,} markers")
        
        # Quick check for inversions
        position_diffs = comparison_df['inverted_position'].diff()
        decreasing_count = (position_diffs < 0).sum()
        logger.info(f"  Found {decreasing_count} decreasing positions in marker order")
        
    else:
        logger.info("Independent mode: Sequence-based analysis")
        logger.info("Comparing actual sequence content between genomes...")
        
        # Use sequence-based inversion detection
        detection_result = detect_sequence_based_inversions(
            original_result['markers_df'], 
            inverted_result['markers_df']
        )
        
        comparison_df = detection_result['comparison']
        
        logger.info("✓ Completed sequence-based comparison")
        logger.info(f"  Found {detection_result['stats']['changed_markers']} markers with sequence changes")
        logger.info(f"  Found {detection_result['stats']['reversed_markers']} markers with reversed sequences")
    
    logger.info(f"  Generated comparison for {len(comparison_df):,} markers")
    
    # Step 2: Detect inversions using sequence-based algorithm
    logger.info("\nStep 2: Detecting inversions")
    if not use_marker_file:
        # Already detected in sequence-based analysis above
        pass  # detection_result already set
    else:
        detection_result = detect_inversions_improved(comparison_df)
    
    # Step 3: Calculate frequencies
    logger.info("\nStep 3: Calculating frequencies")
    frequency_result = calculate_inversion_frequency_improved(
        detection_result['inversions'], 
        comparison_df
    )
    
    # Step 4: Create detailed analysis
    logger.info("\nStep 4: Creating detailed analysis")
    detailed_analysis = create_detailed_analysis_improved(
        detection_result['inversions'],
        original_result['markers_df'],
        config
    )
    
    # Step 5: Create visualizations
    logger.info("\nStep 5: Creating visualizations")
    create_comparison_dotplot(
        comparison_df,
        config['comparison_dotplot_png'],
        "Original vs Inverted Genome Comparison"
    )
    
    create_inversion_heatmap(
        frequency_result,
        original_result['markers_df'],
        config['inversion_heatmap_png']
    )
    
    # Step 6: Save CSV outputs
    logger.info("\nStep 6: Saving analysis results")
    
    # Inversion summary
    if len(detailed_analysis) > 0:
        detailed_analysis.to_csv(config['inversion_summary_csv'], index=False)
        logger.info(f"Inversion summary saved to: {config['inversion_summary_csv']}")
        
        # Frequency analysis
        frequency_result.to_csv(config['inversion_frequency_csv'], index=False)
        logger.info(f"Frequency analysis saved to: {config['inversion_frequency_csv']}")
        
        # Detailed analysis with all metrics
        comprehensive_analysis = detailed_analysis.copy()
        comprehensive_analysis['total_inversions'] = len(detailed_analysis)
        comprehensive_analysis['total_markers'] = len(original_markers)
        comprehensive_analysis['genome_size_mb'] = original_result['markers_df']['size_mb'].sum()
        comprehensive_analysis['inversion_density'] = len(detailed_analysis) / comprehensive_analysis['genome_size_mb']
        
        comprehensive_analysis.to_csv(config['detailed_analysis_csv'], index=False)
        logger.info(f"Detailed analysis saved to: {config['detailed_analysis_csv']}")
    else:
        logger.info("No inversions detected - creating empty output files")
        
        # Create empty DataFrames with proper column structure
        empty_summary = pd.DataFrame({
            'inversion_id': [],
            'chromosome': [],
            'start_bp': [],
            'end_bp': [],
            'size_markers': [],
            'size_mb': [],
            'num_genes': [],
            'inversion_rate_per_mb': [],
            'type': [],
            'inversion_events': []
        })
        
        empty_frequency = pd.DataFrame({
            'marker_id': [],
            'frequency': [],
            'total_events': []
        })
        
        # Write empty files safely
        try:
            empty_summary.to_csv(config['inversion_summary_csv'], index=False)
            logger.info(f"Empty summary saved to: {config['inversion_summary_csv']}")
        except Exception as e:
            logger.error(f"Error saving summary: {e}")
        
        try:
            empty_frequency.to_csv(config['inversion_frequency_csv'], index=False)
            logger.info(f"Empty frequency saved to: {config['inversion_frequency_csv']}")
        except Exception as e:
            logger.error(f"Error saving frequency: {e}")
        
        try:
            empty_summary.to_csv(config['detailed_analysis_csv'], index=False)
            logger.info(f"Empty detailed analysis saved to: {config['detailed_analysis_csv']}")
        except Exception as e:
            logger.error(f"Error saving detailed analysis: {e}")
    
    # Step 7: Generate summary report
    logger.info("\n=== ANALYSIS COMPLETED ===")
    logger.info("Summary Statistics:")
    logger.info(f"  - Total inversions detected: {len(detection_result['inversions'])}")
    logger.info(f"  - Total markers analyzed: {len(original_markers):,}")
    logger.info(f"  - Markers involved in inversions: {len(frequency_result)}")
    
    if len(detailed_analysis) > 0:
        avg_size = detailed_analysis['size_markers'].mean()
        logger.info(f"  - Average inversion size: {avg_size:.2f} markers")
        logger.info(f"  - Largest inversion: {detailed_analysis['size_markers'].max()} markers")
        logger.info(f"  - Smallest inversion: {detailed_analysis['size_markers'].min()} markers")
        logger.info(f"  - Total inverted sequence: {detailed_analysis['size_mb'].sum():.2f} Mb")
    else:
        logger.info("  - Average inversion size: 0 markers")
    
    logger.info("\nFiles created:")
    logger.info(f"  - Comparison dotplot: {config['comparison_dotplot_png']}")
    if len(frequency_result) > 0:
        logger.info(f"  - Inversion heatmap: {config['inversion_heatmap_png']}")
    logger.info(f"  - Inversion summary: {config['inversion_summary_csv']}")
    logger.info(f"  - Frequency analysis: {config['inversion_frequency_csv']}")
    logger.info(f"  - Detailed analysis: {config['detailed_analysis_csv']}")
    
    return {
        'config': config,
        'original_data': original_result,
        'detection_result': detection_result,
        'frequency_result': frequency_result,
        'detailed_analysis': detailed_analysis,
        'comparison_df': comparison_df
    }

################################################################################
# MAIN EXECUTION
################################################################################

if __name__ == "__main__":
    # Execute the analysis
    analysis_results = run_genome_inversion_analysis(ANALYSIS_CONFIG)
    
    # Print detailed summary
    print("\n" + "=" * 50)
    print("DETAILED INVERSION ANALYSIS SUMMARY")
    print("=" * 50)
    
    if len(analysis_results['detailed_analysis']) > 0:
        print("Inversion Details:")
        print(analysis_results['detailed_analysis'].to_string(index=False))
        
        print("\nFrequency Analysis:")
        print(analysis_results['frequency_result'].head(10).to_string(index=False))
    else:
        print("No inversions detected in the comparison.")