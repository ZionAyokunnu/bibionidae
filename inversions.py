#!/usr/bin/env python3
"""
Simple runner for Genome Inversion Analyzer
"""

import sys
from pathlib import Path
import argparse

# Add the current directory to Python path 
sys.path.insert(0, str(Path(__file__).parent))

from genome_inversion_analyser import run_complete_enhanced_analysis_with_hybrid
from genome_inversion_analyser.config import (
    COMPLETE_ENHANCED_CONFIG, 
    ENHANCED_HYBRID_CONFIG, 
    FAST_HYBRID_CONFIG
)

def main():
    parser = argparse.ArgumentParser(description='Run Genome Inversion Analysis')
    parser.add_argument('--mode', choices=['fast', 'hybrid', 'complete'], 
                       default='complete', help='Analysis mode')
    args = parser.parse_args()
    
    # Choose config based on mode
    if args.mode == 'fast':
        config = FAST_HYBRID_CONFIG.copy()
        print("🚀 Using FAST mode (Biopython only, minimal features)")
    elif args.mode == 'hybrid':
        config = ENHANCED_HYBRID_CONFIG.copy()
        print("⚡ Using HYBRID mode (Minimap2 + Biopython)")
    else:
        config = COMPLETE_ENHANCED_CONFIG.copy()
        print("🔬 Using COMPLETE mode (All features)")

    # Update config with your file paths
    config.update({
        'first_fasta_path': 'GCA_910594885.2_idBibMarc1.2_genomic.fna',
        'second_fasta_path': 'GCA_958336335.1_idDilFebr1.1_genomic.fna', 
        'first_busco_path': 'Bibio_marci/full_table.tsv',
        'second_busco_path': 'Dilophus_febrilis/full_table.tsv',
        'base_output_dir': f'v4/{args.mode}_results'  # Different output dirs per mode
    })
    
    # Verify files exist
    required_files = [
        config['first_fasta_path'],
        config['second_fasta_path'], 
        config['first_busco_path'],
        config['second_busco_path']
    ]
    
    for file_path in required_files:
        if not Path(file_path).exists():
            print(f"❌ File not found: {file_path}")
            return None
    
    try:
        print("🧬 Starting Genome Inversion Analysis...")
        results = run_complete_enhanced_analysis_with_hybrid(config)
        print("✅ Analysis Complete!")
        return results
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        return None

if __name__ == "__main__":
    main()