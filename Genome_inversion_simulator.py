#!/usr/bin/env python3
"""
MODULE 1: GENOME INVERSION SIMULATOR (Python Version)
Purpose: Load FASTA, simulate inversions, create dotplot, save inverted FASTA
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import os
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

################################################################################
# CONFIGURATION PARAMETERS
################################################################################

CONFIG = {
    # Input/Output paths
    'input_fasta_path': 'GCA_910594885.2_idBibMarc1.2_genomic.fna',
    'inverted_fasta_output': 'Bibio_marci_inverted.fasta',
    'dotplot_output': 'manual_inversion_dotplot.png',
    
    # Simulation parameters
    'simulate_inversions': True,
    'num_inversions': 1,
    'marker_size_bp': 5000,  # Size of each marker block in base pairs
    
    # NEW: Inversion type controls
    'inversion_type': 'both',  # 'nucleotide', 'marker', or 'both'
    'nucleotide_level_inversion': True,   # Reverse DNA sequences within markers
    'marker_level_inversion': True,       # Rearrange marker order
    
    # Optional features
    'use_busco': False,
    'save_intermediate_states': True
}

################################################################################
# CORE FUNCTIONS
################################################################################

def load_and_process_genome(fasta_path, marker_size=1000):
    """
    Load and process FASTA file into genome markers
    
    Args:
        fasta_path (str): Path to input FASTA file
        marker_size (int): Size of each marker block in base pairs
        
    Returns:
        tuple: (genome_data, genome_df)
    """
    logger.info(f"Loading genome from: {fasta_path}")
    
    # Read FASTA sequences
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
    
    logger.info(f"Loaded {len(genome_df)} scaffolds/chromosomes")
    logger.info(f"Total genome size: {genome_df['chromosome_size_bp'].sum():,} bp")
    
    # Create markers (1 marker per specified bp block)
    markers_list = []
    
    for _, row in genome_df.iterrows():
        chromosome = row['chromosome']
        chr_size = row['chromosome_size_bp']
        sequence = row['sequence']
        
        marker_count = chr_size // marker_size
        
        if marker_count > 0:  # Skip very small scaffolds
            for marker_idx in range(1, marker_count + 1):
                start_bp = (marker_idx - 1) * marker_size + 1
                end_bp = marker_idx * marker_size
                marker_id = f"{chromosome}_{marker_idx:06d}"
                
                markers_list.append({
                    'chromosome': chromosome,
                    'marker_index': marker_idx,
                    'start_bp': start_bp,
                    'end_bp': end_bp,
                    'marker_id': marker_id,
                    'sequence': sequence
                })
    
    markers_df = pd.DataFrame(markers_list)
    
    # Add cumulative positions across all chromosomes
    markers_df = markers_df.sort_values(['chromosome', 'marker_index']).reset_index(drop=True)
    
    # Calculate linear positions
    chr_offsets = {}
    current_offset = 0
    
    for chromosome in markers_df['chromosome'].unique():
        chr_offsets[chromosome] = current_offset
        chr_markers = markers_df[markers_df['chromosome'] == chromosome]
        current_offset += len(chr_markers) * marker_size
    
    markers_df['chr_offset'] = markers_df['chromosome'].map(chr_offsets)
    markers_df['linear_start'] = markers_df['chr_offset'] + (markers_df['marker_index'] - 1) * marker_size
    markers_df['linear_end'] = markers_df['linear_start'] + marker_size
    
    logger.info(f"Created {len(markers_df):,} markers")
    
    return genome_df, markers_df

def simulate_single_inversion(genome_vec):
    """
    Simulate a single inversion event
    
    Args:
        genome_vec (list): Vector of gene/marker identifiers
        
    Returns:
        dict: Dictionary with new genome and inversion coordinates
    """
    n = len(genome_vec)
    if n < 2:
        raise ValueError("Genome too small for inversion")
    
    # Choose random inversion boundaries
    start_idx = random.randint(0, n - 2)  # 0-based indexing
    end_idx = random.randint(start_idx + 1, n - 1)
    
    # Perform inversion
    new_genome = genome_vec.copy()
    new_genome[start_idx:end_idx + 1] = new_genome[start_idx:end_idx + 1][::-1]
    
    return {
        'new_genome': new_genome,
        'inversion_start': start_idx + 1,  # Convert to 1-based for output
        'inversion_end': end_idx + 1,
        'inversion_size': end_idx - start_idx + 1
    }

def simulate_multiple_inversions(genome_vec, num_inversions=1):
    """
    Simulate multiple inversion events
    
    Args:
        genome_vec (list): Vector of gene/marker identifiers
        num_inversions (int): Number of inversions to simulate
        
    Returns:
        dict: Dictionary with final genome, inversion log, and history
    """
    logger.info(f"Simulating {num_inversions} inversion(s)")
    
    current_genome = genome_vec.copy()
    inversions_log = []
    genome_history = [current_genome.copy()]  # Store original
    
    for i in range(num_inversions):
        result = simulate_single_inversion(current_genome)
        
        # Log this inversion (convert to 1-based indexing for consistency with R)
        inversions_log.append({
            'inversion_id': i + 1,
            'start_idx': result['inversion_start'],
            'end_idx': result['inversion_end'],
            'size': result['inversion_size'],
            'original_start_marker': genome_vec[result['inversion_start'] - 1],  # Original marker at this position
            'original_end_marker': genome_vec[result['inversion_end'] - 1]
        })
        
        current_genome = result['new_genome']
        genome_history.append(current_genome.copy())
        
        logger.info(f"  Inversion {i + 1}: {result['inversion_start']} to {result['inversion_end']} "
                   f"({result['inversion_size']} markers)")
    
    return {
        'original_genome': genome_vec,
        'final_genome': current_genome,
        'inversions_log': pd.DataFrame(inversions_log),
        'genome_history': genome_history
    }

def create_inversion_dotplot(original_genome, inverted_genome, output_path, 
                           title="Genome Inversion Dotplot"):
    """
    Create dotplot comparing original and inverted genomes
    
    Args:
        original_genome (list): Original gene order
        inverted_genome (list): Inverted gene order
        output_path (str): Path to save plot
        title (str): Plot title
    """
    logger.info("Creating dotplot...")
    
    # Create plotting data
    inverted_positions = []
    for marker in original_genome:
        try:
            inverted_positions.append(inverted_genome.index(marker) + 1)  # 1-based indexing
        except ValueError:
            inverted_positions.append(np.nan)  # Marker not found
    
    plot_data = pd.DataFrame({
        'original_index': range(1, len(original_genome) + 1),
        'inverted_index': inverted_positions,
        'marker_id': original_genome
    })
    
    # Remove NaN values
    plot_data = plot_data.dropna()
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    plt.scatter(plot_data['original_index'], plot_data['inverted_index'], 
               color='darkred', s=1.5, alpha=0.7)
    
    # Add diagonal reference line
    max_val = max(plot_data['original_index'].max(), plot_data['inverted_index'].max())
    plt.plot([1, max_val], [1, max_val], color='gray', linestyle='--', alpha=0.7)
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.subtitle = f"Total markers: {len(original_genome):,}"
    plt.xlabel("Original Genome Position (Index)", fontsize=11)
    plt.ylabel("Inverted Genome Position (Index)", fontsize=11)
    
    # Add subtitle manually
    plt.figtext(0.5, 0.92, f"Total markers: {len(original_genome):,}", 
               ha='center', fontsize=12)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Dotplot saved to: {output_path}")

def save_inverted_fasta(genome_data, original_genome, inverted_genome, 
                       markers_df, output_path, inversions_log, config):
    """
    Save inverted genome to FASTA format with configurable inversion types
    
    Args:
        genome_data (DataFrame): Original genome data with sequences
        original_genome (list): Original marker order
        inverted_genome (list): New order of markers after inversions
        markers_df (DataFrame): Genome dataframe with marker information
        output_path (str): Path to save inverted FASTA
        inversions_log (DataFrame): Log of inversions performed
        config (dict): Configuration including inversion type settings
        
    Returns:
        dict: Summary information about the reconstruction
    """
    logger.info("Saving inverted genome to FASTA with configurable inversion reconstruction...")
    logger.info(f"Inversion type: {config.get('inversion_type', 'both')}")
    logger.info(f"Nucleotide level: {config.get('nucleotide_level_inversion', True)}")
    logger.info(f"Marker level: {config.get('marker_level_inversion', True)}")
    logger.info("Inversion details:")
    print(inversions_log.to_string(index=False))
    
    # Start with original sequences
    reconstructed_sequences = {}
    for _, row in genome_data.iterrows():
        reconstructed_sequences[row['chromosome']] = row['sequence']
    
    # STEP 1: Apply marker-level inversions if enabled
    if config.get('marker_level_inversion', True):
        logger.info("Applying MARKER-LEVEL inversions (rearranging blocks)...")
        
        # Reconstruct sequences based on inverted marker order
        for chrom in reconstructed_sequences.keys():
            # Get markers for this chromosome in original order
            chrom_markers_orig = markers_df[markers_df['chromosome'] == chrom].copy()
            chrom_markers_orig = chrom_markers_orig.sort_values('start_bp').reset_index(drop=True)
            
            if len(chrom_markers_orig) == 0:
                continue
                
            # Get markers in inverted order for this chromosome
            inverted_order_df = pd.DataFrame({
                'new_order': range(len(inverted_genome)),
                'marker_id': inverted_genome
            })
            
            chrom_markers_inv = chrom_markers_orig.merge(inverted_order_df, on='marker_id', how='inner')
            
            if len(chrom_markers_inv) > 0:
                # Sort by new order to get rearranged sequence
                chrom_markers_inv = chrom_markers_inv.sort_values('new_order')
                
                # Reconstruct this chromosome's sequence from rearranged markers
                original_seq = reconstructed_sequences[chrom]
                new_seq = ""
                
                logger.info(f"  Rearranging {len(chrom_markers_inv)} markers for {chrom}")
                
                for _, marker in chrom_markers_inv.iterrows():
                    start_bp = max(0, marker['start_bp'] - 1)  # Convert to 0-based
                    end_bp = min(marker['end_bp'], len(original_seq))
                    
                    if start_bp < end_bp and start_bp < len(original_seq):
                        segment = original_seq[start_bp:end_bp]
                        new_seq += segment
                
                if len(new_seq) > 0:
                    reconstructed_sequences[chrom] = new_seq
                    logger.info(f"    Rearranged {chrom}: {len(new_seq):,} bp")
    
    # STEP 2: Apply nucleotide-level inversions if enabled
    if config.get('nucleotide_level_inversion', True):
        logger.info("Applying NUCLEOTIDE-LEVEL inversions (reversing sequences)...")
        
        # Apply each inversion to reverse sequences within regions
        for _, inversion in inversions_log.iterrows():
            logger.info(f"Applying nucleotide inversion {inversion['inversion_id']}:")
            logger.info(f"  Start marker index: {inversion['start_idx']}")
            logger.info(f"  End marker index: {inversion['end_idx']}")
            logger.info(f"  Size: {inversion['size']} markers")
            
            # Get the chromosomes and positions for start and end markers
            start_marker_info = markers_df[markers_df['marker_id'] == inversion['original_start_marker']]
            end_marker_info = markers_df[markers_df['marker_id'] == inversion['original_end_marker']]
            
            if len(start_marker_info) == 0 or len(end_marker_info) == 0:
                logger.warning("  Warning: Could not find marker positions, skipping inversion")
                continue
            
            start_marker_info = start_marker_info.iloc[0]
            end_marker_info = end_marker_info.iloc[0]
            
            logger.info(f"  Start marker: {start_marker_info['marker_id']}")
            logger.info(f"  End marker: {end_marker_info['marker_id']}")
            
            # Single chromosome inversion - reverse nucleotide sequence
            if start_marker_info['chromosome'] == end_marker_info['chromosome']:
                chrom = start_marker_info['chromosome']
                current_seq = reconstructed_sequences[chrom]
                
                # Convert marker indices to base pair positions (convert to 0-based)
                start_bp = max(0, start_marker_info['start_bp'] - 1)
                end_bp = min(end_marker_info['end_bp'], len(current_seq))
                
                logger.info(f"  Reversing nucleotides in {chrom} from bp {start_bp + 1} to {end_bp}")
                
                if start_bp < end_bp:
                    # Extract segments
                    before_segment = current_seq[:start_bp] if start_bp > 0 else ""
                    inversion_segment = current_seq[start_bp:end_bp]
                    after_segment = current_seq[end_bp:] if end_bp < len(current_seq) else ""
                    
                    # Reverse the inversion segment (nucleotide level)
                    reversed_segment = inversion_segment[::-1]
                    
                    # Reconstruct the sequence
                    new_seq = before_segment + reversed_segment + after_segment
                    reconstructed_sequences[chrom] = new_seq
                    
                    logger.info(f"    Reversed {len(inversion_segment):,} nucleotides in {chrom}")
                    logger.info(f"    New sequence length: {len(new_seq):,} bp")
    
    # Save reconstructed sequences
    records = []
    for chrom, seq in reconstructed_sequences.items():
        record = SeqRecord(Seq(seq), id=chrom, description="")
        records.append(record)
    
    SeqIO.write(records, output_path, "fasta")
    
    # Save marker order file for verification
    marker_order_file = output_path.replace('.fasta', '_marker_order.csv')
    
    # Map original positions to inverted positions
    original_to_inverted_pos = []
    for marker in original_genome:
        try:
            original_to_inverted_pos.append(inverted_genome.index(marker) + 1)  # 1-based
        except ValueError:
            original_to_inverted_pos.append(np.nan)
    
    inverted_marker_df = pd.DataFrame({
        'original_position': range(1, len(original_genome) + 1),
        'marker_id': original_genome,
        'inverted_position': original_to_inverted_pos
    })
    
    # Join with marker information
    marker_info = markers_df[['marker_id', 'chromosome', 'start_bp', 'end_bp']]
    inverted_marker_df = inverted_marker_df.merge(marker_info, on='marker_id', how='left')
    inverted_marker_df = inverted_marker_df.sort_values('original_position')
    
    # Verify the inversion is captured
    position_diffs = inverted_marker_df['inverted_position'].diff()
    inversions_detected = (position_diffs < 0).sum()
    
    logger.info(f"Detected {inversions_detected} decreasing positions in marker order")
    
    inverted_marker_df.to_csv(marker_order_file, index=False)
    logger.info(f"Marker order saved to: {marker_order_file}")
    
    # Save inversion log
    inversion_log_file = output_path.replace('.fasta', '_inversion_log.csv')
    inversions_log.to_csv(inversion_log_file, index=False)
    logger.info(f"Inversion log saved to: {inversion_log_file}")
    
    logger.info(f"Inverted FASTA with configurable reconstruction saved to: {output_path}")
    
    # Calculate and report size differences
    original_total_size = genome_data['chromosome_size_bp'].sum()
    reconstructed_total_size = sum(len(seq) for seq in reconstructed_sequences.values())
    
    logger.info(f"Original genome size: {original_total_size:,} bp")
    logger.info(f"Reconstructed genome size: {reconstructed_total_size:,} bp")
    
    if abs(original_total_size - reconstructed_total_size) < 1000:
        logger.info("✓ Sequence reconstruction successful (size difference < 1kb)")
    else:
        logger.warning("⚠ Warning: Significant size difference in reconstruction")
    
    return {
        'marker_order_file': marker_order_file,
        'inversion_log_file': inversion_log_file,
        'inversions_in_order': inversions_detected,
        'original_size': original_total_size,
        'reconstructed_size': reconstructed_total_size,
        'inversion_type': config.get('inversion_type', 'both')
    }

################################################################################
# MAIN EXECUTION FUNCTION
################################################################################

def run_genome_inversion_simulation(config=None):
    """
    Main function to run the genome inversion simulation
    
    Args:
        config (dict): Configuration dictionary with parameters
        
    Returns:
        dict: Results of the simulation
    """
    if config is None:
        config = CONFIG
    
    logger.info("=== GENOME INVERSION SIMULATOR ===")
    logger.info(f"Starting simulation with {config['num_inversions']} inversion(s)\n")
    
    # Step 1: Load and process genome
    genome_data, markers_df = load_and_process_genome(
        config['input_fasta_path'], 
        config['marker_size_bp']
    )
    
    # Step 2: Simulate inversions (if enabled)
    if config['simulate_inversions']:
        marker_vector = markers_df['marker_id'].tolist()
        
        simulation_result = simulate_multiple_inversions(
            marker_vector, 
            config['num_inversions']
        )
        
        # Step 3: Create dotplot
        create_inversion_dotplot(
            simulation_result['original_genome'],
            simulation_result['final_genome'],
            config['dotplot_output'],
            f"Manual Inversion Dotplot - {config['num_inversions']} Inversion(s)"
        )
        
        # Step 4: Save inverted FASTA
        save_result = save_inverted_fasta(
            genome_data,
            simulation_result['original_genome'],
            simulation_result['final_genome'],
            markers_df,
            config['inverted_fasta_output'],
            simulation_result['inversions_log'],
            config  # Pass config for inversion type settings
        )
        
        logger.info(f"Marker order verification: found {save_result['inversions_in_order']} inversions in saved order")
        
        # Step 5: Return results for potential further analysis
        logger.info("\n=== SIMULATION COMPLETED ===")
        logger.info("Files created:")
        logger.info(f"  - Dotplot: {config['dotplot_output']}")
        logger.info(f"  - Inverted FASTA: {config['inverted_fasta_output']}")
        
        return {
            'config': config,
            'genome_data': genome_data,
            'markers_df': markers_df,
            'simulation_result': simulation_result,
            'save_result': save_result
        }
    else:
        logger.info("Simulation disabled. Only genome loading performed.")
        return {
            'config': config,
            'genome_data': genome_data,
            'markers_df': markers_df
        }

################################################################################
# MAIN EXECUTION
################################################################################

if __name__ == "__main__":
    # Set random seed for reproducibility (optional)
    random.seed(42)
    np.random.seed(42)
    
    # Execute the simulation
    results = run_genome_inversion_simulation(CONFIG)
    
    # Print summary
    if CONFIG['simulate_inversions']:
        print("\nInversion Summary:")
        print(results['simulation_result']['inversions_log'].to_string(index=False))