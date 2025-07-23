# from genome_inversion_analyser.config import FAST_HYBRID_CONFIG
# from genome_inversion_analyser.workflow import run_complete_analysis

# if __name__ == "__main__":
#     config = FAST_HYBRID_CONFIG
#     run_complete_analysis(config)



# =============================================================================
# Main Entry Point (main.py)
# =============================================================================

"""
Main entry point for the genome inversion analysis pipeline.
Handles CLI argument parsing and workflow orchestration.
"""

import sys
import argparse
from pathlib import Path

from genome_inversion_analyser.config import (
    ENHANCED_HYBRID_CONFIG, 
    FAST_HYBRID_CONFIG, 
    COMPLETE_ENHANCED_CONFIG
)
from genome_inversion_analyser.logger import setup_logger, get_logger

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Enhanced Genome Synteny and Inversion Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m genome_inversion_analyser --mode fast
  python -m genome_inversion_analyser --mode hybrid --threads 8
  python -m genome_inversion_analyser --mode complete --output results/
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
        # Import and run the main workflow
        # (This will be implemented in Phase 2)
        logger.info("Phase 1 infrastructure setup complete!")
        logger.info("Next: Implement workflow orchestration in Phase 2")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

