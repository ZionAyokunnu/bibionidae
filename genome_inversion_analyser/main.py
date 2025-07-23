#!/usr/bin/env python3
"""
Main entry point for the genome inversion analysis pipeline.
Handles CLI argument parsing and workflow orchestration.
"""

import sys
import argparse
import logging
from pathlib import Path

from .config import (
    ENHANCED_HYBRID_CONFIG, 
    FAST_HYBRID_CONFIG, 
    COMPLETE_ENHANCED_CONFIG
)
from .logger import setup_logger, get_logger

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Enhanced Genome Synteny and Inversion Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m genome_inversion_analyser_v5 --mode fast
  python -m genome_inversion_analyser_v5 --mode hybrid --threads 8
  python -m genome_inversion_analyser_v5 --mode complete --output results/
        """
    )
    
    parser.add_argument(
        '--mode', 
        choices=['fast', 'hybrid', 'complete'],
        default='hybrid',
        help='Analysis mode: fast (Biopython only), hybrid (Minimap2+Biopython), complete (all features)'
    )
    
    parser.add_argument(
        '--first-genome', 
        type=str,
        help='Path to first genome FASTA file'
    )
    
    parser.add_argument(
        '--second-genome', 
        type=str,
        help='Path to second genome FASTA file'
    )
    
    parser.add_argument(
        '--first-busco', 
        type=str,
        help='Path to first genome BUSCO full_table.tsv'
    )
    
    parser.add_argument(
        '--second-busco', 
        type=str,
        help='Path to second genome BUSCO full_table.tsv'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '--threads', '-t',
        type=int,
        help='Number of threads to use for parallel processing'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--log-file',
        type=str,
        help='Path to log file'
    )
    
    return parser.parse_args()

def select_config(mode):
    """Select configuration based on mode."""
    configs = {
        'fast': FAST_HYBRID_CONFIG,
        'hybrid': ENHANCED_HYBRID_CONFIG,
        'complete': COMPLETE_ENHANCED_CONFIG
    }
    return configs[mode]

def override_config_from_args(config, args):
    """Override configuration values from command line arguments."""
    if args.first_genome:
        config.update({'first_fasta_path': args.first_genome})
    if args.second_genome:
        config.update({'second_fasta_path': args.second_genome})
    if args.first_busco:
        config.update({'first_busco_path': args.first_busco})
    if args.second_busco:
        config.update({'second_busco_path': args.second_busco})
    if args.output:
        config.update({'base_output_dir': args.output})
    if args.threads:
        config.update({'minimap2_threads': args.threads})

def run_enhanced_analysis(config):
    """
    Main analysis workflow - placeholder for Phase 2 implementation.
    
    Args:
        config: AnalysisConfig object with analysis parameters
        
    Returns:
        Dictionary with analysis results
    """
    logger = get_logger()
    
    logger.info("Starting enhanced genome analysis...")
    logger.info(f"Input genomes: {config.get('first_fasta_path')} vs {config.get('second_fasta_path')}")
    logger.info(f"BUSCO tables: {config.get('first_busco_path')} vs {config.get('second_busco_path')}")
    logger.info(f"Alignment strategy: {config.get('alignment_strategy')}")
    
    # Placeholder results structure
    results = {
        'status': 'success',
        'config': config.to_dict(),
        'synteny_blocks': [],
        'inversions': [],
        'rearrangements': [],
        'quality_metrics': {},
        'statistics': {}
    }
    
    logger.info("Analysis workflow completed (Phase 1 placeholder)")
    return results

def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logger(level=log_level, log_file=args.log_file)
    
    # Select and configure analysis
    config = select_config(args.mode)
    override_config_from_args(config, args)
    
    logger.section_header(f"ENHANCED GENOME INVERSION ANALYSIS - {args.mode.upper()} MODE")
    
    logger.info(f"Configuration: {args.mode} mode")
    logger.info(f"Alignment strategy: {config.get('alignment_strategy')}")
    logger.info(f"Output directory: {config.get('base_output_dir')}")
    
    try:
        # Run the main workflow
        results = run_enhanced_analysis(config)
        logger.info("Analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()