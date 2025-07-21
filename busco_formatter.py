#!/usr/bin/env python3
"""
BUSCO Data Formatter for Syngraph
Converts BUSCO full_table.tsv files to syngraph format
"""

import pandas as pd
import os
import sys

def format_busco_for_syngraph(busco_file, output_file, taxon_name):
    """
    Format BUSCO full_table.tsv for syngraph input
    
    Args:
        busco_file (str): Path to BUSCO full_table.tsv
        output_file (str): Output file path
        taxon_name (str): Taxon name for the file
    """
    print(f"Processing {busco_file} for taxon {taxon_name}")
    
    # Read BUSCO table
    busco_data = []
    with open(busco_file, 'r') as f:
        for line in f:
            if not line.startswith('#') and line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 5 and parts[1] == 'Complete':
                    busco_data.append({
                        'busco_id': parts[0],
                        'status': parts[1],
                        'sequence': parts[2],
                        'start': parts[3],
                        'end': parts[4]
                    })
    
    # Create DataFrame
    df = pd.DataFrame(busco_data)
    
    # Filter for complete BUSCOs only
    df = df[df['status'] == 'Complete']
    
    # Remove rows with 'N/A' coordinates
    df = df[(df['start'] != 'N/A') & (df['end'] != 'N/A')]
    
    # Convert coordinates to integers
    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])
    
    # Sort by sequence and position
    df = df.sort_values(['sequence', 'start'])
    
    # Create syngraph format: BUSCO_ID, Sequence, Gene_Start, Gene_End
    output_df = df[['busco_id', 'sequence', 'start', 'end']]
    
    # Save to file
    output_df.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"  Saved {len(output_df)} complete BUSCO genes to {output_file}")
    return len(output_df)

def main():
    """Main function to process BUSCO files"""
    
    # Configuration - adjust these paths
    busco_files = {
        'Bibio_marci': 'busco-data/Bibio_marci.tsv',
        'Dilophus_febrilis': 'busco-data/Dilophus_febrilis.tsv',
        'Dioctria_linearis': 'busco-data/Dioctria_linearisis.tsv',
        'Dioctria_rufipes': 'busco-data/Dioctria_rufipes.tsv',
        'Plecia_longiforceps': 'busco-data/Plecia_longiforceps.tsv',
    }
    
    output_dir = 'syngraph_input'
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== BUSCO Data Formatter for Syngraph ===")
    
    total_genes = 0
    for taxon_name, busco_file in busco_files.items():
        if os.path.exists(busco_file):
            output_file = os.path.join(output_dir, f"{taxon_name}.tsv")
            genes_count = format_busco_for_syngraph(busco_file, output_file, taxon_name)
            total_genes += genes_count
        else:
            print(f"Warning: {busco_file} not found!")
    
    print(f"\nTotal genes processed: {total_genes}")
    print(f"Output files created in: {output_dir}/")
    print("\nNext steps:")
    print("1. Create a newick tree file with your taxa")
    print("2. Run syngraph build")
    print("3. Run syngraph infer")

if __name__ == "__main__":
    main()