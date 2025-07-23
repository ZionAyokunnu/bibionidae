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

# Import all the analysis functions that were in the unmodularized version
# These need to be implemented in separate modules or imported from core.py
from .core import (
    create_output_directory,
    enhanced_parse_busco_table,
    assess_assembly_quality,
    enhanced_filter_busco_genes,
    setup_hybrid_sequence_aligner,
    extract_enhanced_busco_sequences,
    generate_cache_key,
    load_cached_alignment_results,
    run_hybrid_alignment_analysis,
    cache_alignment_results,
    analyze_enhanced_synteny_blocks,
    analyze_enhanced_chromosome_rearrangements,
    analyze_enhanced_inversions,
    save_enhanced_results,
    create_enhanced_visualizations,
    generate_comprehensive_report
)

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
    Main analysis workflow
    
    Args:
        config: Configuration dictionary with analysis parameters
        
    Returns:
        Dictionary with analysis results
    """
    logger = get_logger()
    
    logger.info("=" * 80)
    logger.info("ENHANCED INTEGRATED SYNTENY AND INVERSION ANALYZER")
    logger.info("With Hybrid Alignment System (Minimap2 + Biopython)")
    logger.info("=" * 80)
    
    # Create output directory
    output_dir = create_output_directory(config)
    logger.info(f"Output directory: {output_dir}")
    
    try:
        # Phase 1: Enhanced BUSCO Processing and Quality Assessment
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 1: ENHANCED BUSCO PROCESSING")
        logger.info("=" * 50)
        
        # Parse BUSCO data with enhanced validation
        logger.info("Step 1.1: Parsing BUSCO tables")
        first_busco_raw = enhanced_parse_busco_table(config['first_busco_path'], config)
        second_busco_raw = enhanced_parse_busco_table(config['second_busco_path'], config)
        
        # Assess assembly quality
        logger.info("Step 1.2: Assessing assembly quality")
        first_quality = assess_assembly_quality(config['first_fasta_path'], first_busco_raw, config)
        second_quality = assess_assembly_quality(config['second_fasta_path'], second_busco_raw, config)
        
        # Enhanced BUSCO filtering
        logger.info("Step 1.3: Enhanced BUSCO filtering")
        first_busco_filtered = enhanced_filter_busco_genes(first_busco_raw, config, first_quality)
        second_busco_filtered = enhanced_filter_busco_genes(second_busco_raw, config, second_quality)
        
        # Enhanced sequence extraction
        logger.info("Step 1.4: Enhanced sequence extraction")
        aligner = setup_hybrid_sequence_aligner(config)
        first_busco_seqs = extract_enhanced_busco_sequences(first_busco_filtered, config['first_fasta_path'], config)
        second_busco_seqs = extract_enhanced_busco_sequences(second_busco_filtered, config['second_fasta_path'], config)
        
        # Phase 2: Enhanced Ortholog Mapping with Hybrid Alignment
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 2: ENHANCED ORTHOLOG MAPPING (HYBRID ALIGNMENT)")
        logger.info("=" * 50)
        
        logger.info("Step 2.1: Creating enhanced ortholog mapping with hybrid alignment")
        
        # Check cache first
        cache_key = generate_cache_key(first_busco_seqs, second_busco_seqs, config)
        cached_result = load_cached_alignment_results(cache_key, config)
        
        if cached_result:
            logger.info("  Using cached alignment results")
            ortholog_df, paralog_df = cached_result
        else:
            # Run hybrid alignment analysis
            ortholog_df, paralog_df = run_hybrid_alignment_analysis(
                first_busco_seqs, second_busco_seqs, config
            )
            
            # Cache results
            cache_alignment_results((ortholog_df, paralog_df), cache_key, config)
        
        # Phase 3: Enhanced Synteny Analysis
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 3: ENHANCED SYNTENY ANALYSIS")
        logger.info("=" * 50)
        
        logger.info("Step 3.1: Analyzing enhanced synteny blocks")
        synteny_df, mapping_df = analyze_enhanced_synteny_blocks(ortholog_df, config)
        
        # Phase 4: Enhanced Rearrangement Analysis
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 4: ENHANCED REARRANGEMENT ANALYSIS")
        logger.info("=" * 50)
        
        logger.info("Step 4.1: Analyzing chromosome rearrangements")
        rearrangement_df = analyze_enhanced_chromosome_rearrangements(ortholog_df, config)
        
        # Phase 5: Enhanced Inversion Analysis
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 5: ENHANCED INVERSION ANALYSIS")
        logger.info("=" * 50)
        
        logger.info("Step 5.1: Analyzing inversions")
        inversion_df = analyze_enhanced_inversions(synteny_df, ortholog_df, config)
        
        # Phase 6: Results Integration and Reporting
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 6: RESULTS INTEGRATION AND REPORTING")
        logger.info("=" * 50)
        
        # Save all results
        logger.info("Step 6.1: Saving analysis results")
        save_enhanced_results(output_dir, {
            'ortholog_df': ortholog_df,
            'paralog_df': paralog_df,
            'synteny_df': synteny_df,
            'mapping_df': mapping_df,
            'rearrangement_df': rearrangement_df,
            'inversion_df': inversion_df,
            'first_quality': first_quality,
            'second_quality': second_quality
        }, config)
        
        # Generate comprehensive visualizations
        logger.info("Step 6.2: Creating enhanced visualizations")
        create_enhanced_visualizations(output_dir, {
            'ortholog_df': ortholog_df,
            'synteny_df': synteny_df,
            'rearrangement_df': rearrangement_df,
            'inversion_df': inversion_df
        }, config)
        
        # Generate final report
        logger.info("Step 6.3: Generating comprehensive report")
        generate_comprehensive_report(output_dir, {
            'ortholog_df': ortholog_df,
            'paralog_df': paralog_df,
            'synteny_df': synteny_df,
            'rearrangement_df': rearrangement_df,
            'inversion_df': inversion_df,
            'first_quality': first_quality,
            'second_quality': second_quality,
            'config': config
        })
        
        logger.info("\n" + "=" * 80)
        logger.info("ENHANCED ANALYSIS WITH HYBRID ALIGNMENT COMPLETED SUCCESSFULLY")
        logger.info("=" * 80)
        
        return {
            'status': 'success',
            'config': config,
            'synteny_blocks': synteny_df.to_dict('records') if not synteny_df.empty else [],
            'inversions': inversion_df.to_dict('records') if not inversion_df.empty else [],
            'rearrangements': rearrangement_df.to_dict('records') if not rearrangement_df.empty else [],
            'quality_metrics': {
                'first_quality': first_quality,
                'second_quality': second_quality
            },
            'statistics': {
                'ortholog_pairs': len(ortholog_df),
                'synteny_blocks_count': len(synteny_df),
                'inversions_count': len(inversion_df),
                'rearrangements_count': len(rearrangement_df)
            }
        }
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise

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