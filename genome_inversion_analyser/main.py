#!/usr/bin/env python3
"""
Enhanced main.py with Registry System Integration
"""

import sys
import random
import logging
from pathlib import Path
import pandas as pd
import numpy as np

# Import all modules including new registry
from genome_inversion_analyser.config import (
    ENHANCED_HYBRID_CONFIG,
    FAST_HYBRID_CONFIG,
    COMPLETE_ENHANCED_CONFIG
)

from genome_inversion_analyser.utils import (
    create_output_directory,
    generate_cache_key,
    cache_alignment_results,
    load_cached_alignment_results
)

from genome_inversion_analyser.core import (
    assess_assembly_quality,
    enhanced_parse_busco_table,
    enhanced_filter_busco_genes,
    extract_enhanced_busco_sequences,
    run_hybrid_alignment_analysis,
    setup_hybrid_sequence_aligner,
    analyze_enhanced_synteny_blocks,
    analyze_enhanced_chromosome_rearrangements,
    analyze_enhanced_inversions
)

from genome_inversion_analyser.visualization import (
    create_enhanced_visualizations
)

# Import new registry system
from genome_inversion_analyser.registry import (
    FileRegistry,
    AnalysisResultsExporter
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def run_complete_enhanced_analysis_with_registry(config=None):
    """Enhanced analysis with integrated registry system"""
    if config is None:
        config = ENHANCED_HYBRID_CONFIG
    
    logger.info("=" * 80)
    logger.info("ENHANCED GENOME INVERSION ANALYZER WITH REGISTRY SYSTEM")
    logger.info("=" * 80)
    
    # Create output directory and registry
    output_dir = create_output_directory(config)
    registry = FileRegistry(output_dir, project_name="genome_inversion_analysis")
    exporter = AnalysisResultsExporter(registry)
    
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Registry initialized: {registry.registry_file}")
    
    try:
        # Phase 1: Enhanced BUSCO Processing and Quality Assessment
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 1: ENHANCED BUSCO PROCESSING WITH REGISTRY")
        logger.info("=" * 50)
        
        # Parse BUSCO data
        logger.info("Step 1.1: Parsing BUSCO tables")
        first_busco_raw = enhanced_parse_busco_table(config['first_busco_path'], config)
        second_busco_raw = enhanced_parse_busco_table(config['second_busco_path'], config)
        
        # Register raw BUSCO data
        registry.register_file(
            'first_busco_raw',
            first_busco_raw,
            'csv',
            'Raw BUSCO results for first genome',
            metadata={'genome': 'first', 'source': config['first_busco_path']}
        )
        
        registry.register_file(
            'second_busco_raw',
            second_busco_raw,
            'csv',
            'Raw BUSCO results for second genome',
            metadata={'genome': 'second', 'source': config['second_busco_path']}
        )
        
        # Assess assembly quality
        logger.info("Step 1.2: Assessing assembly quality")
        first_quality = assess_assembly_quality(config['first_fasta_path'], first_busco_raw, config)
        second_quality = assess_assembly_quality(config['second_fasta_path'], second_busco_raw, config)
        
        # Register quality assessments
        registry.register_file(
            'first_quality_metrics',
            first_quality,
            'json',
            'Assembly quality metrics for first genome',
            dependencies=['first_busco_raw'],
            metadata={'genome': 'first', 'quality_class': first_quality['quality_class']}
        )
        
        registry.register_file(
            'second_quality_metrics',
            second_quality,
            'json',
            'Assembly quality metrics for second genome',
            dependencies=['second_busco_raw'],
            metadata={'genome': 'second', 'quality_class': second_quality['quality_class']}
        )
        
        # Enhanced BUSCO filtering
        logger.info("Step 1.3: Enhanced BUSCO filtering")
        first_busco_filtered = enhanced_filter_busco_genes(first_busco_raw, config, first_quality)
        second_busco_filtered = enhanced_filter_busco_genes(second_busco_raw, config, second_quality)
        
        # Register filtered BUSCO data
        registry.register_file(
            'first_busco_filtered',
            first_busco_filtered,
            'csv',
            'Filtered BUSCO genes for first genome',
            dependencies=['first_busco_raw', 'first_quality_metrics'],
            metadata={'filter_stats': {'initial': len(first_busco_raw), 'filtered': len(first_busco_filtered)}}
        )
        
        registry.register_file(
            'second_busco_filtered',
            second_busco_filtered,
            'csv',
            'Filtered BUSCO genes for second genome',
            dependencies=['second_busco_raw', 'second_quality_metrics'],
            metadata={'filter_stats': {'initial': len(second_busco_raw), 'filtered': len(second_busco_filtered)}}
        )
        
        # Enhanced sequence extraction
        logger.info("Step 1.4: Enhanced sequence extraction")
        aligner = setup_hybrid_sequence_aligner(config)
        first_busco_seqs = extract_enhanced_busco_sequences(first_busco_filtered, config['first_fasta_path'], config)
        second_busco_seqs = extract_enhanced_busco_sequences(second_busco_filtered, config['second_fasta_path'], config)
        
        # Register sequence data
        registry.register_file(
            'first_busco_sequences',
            first_busco_seqs,
            'csv',
            'Extracted BUSCO sequences for first genome',
            dependencies=['first_busco_filtered'],
            metadata={'extraction_stats': getattr(first_busco_seqs, 'attrs', {}).get('extraction_stats', {})}
        )
        
        registry.register_file(
            'second_busco_sequences',
            second_busco_seqs,
            'csv',
            'Extracted BUSCO sequences for second genome',
            dependencies=['second_busco_filtered'],
            metadata={'extraction_stats': getattr(second_busco_seqs, 'attrs', {}).get('extraction_stats', {})}
        )
        
        # Phase 2: Enhanced Ortholog Mapping with Registry
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 2: ORTHOLOG MAPPING WITH CACHING")
        logger.info("=" * 50)
        
        logger.info("Step 2.1: Creating enhanced ortholog mapping")
        
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
        
        # Register ortholog results
        registry.register_file(
            'ortholog_pairs',
            ortholog_df,
            'csv',
            'Orthologous gene pairs between genomes',
            dependencies=['first_busco_sequences', 'second_busco_sequences'],
            metadata={
                'alignment_method': config.get('alignment_strategy', 'unknown'),
                'total_pairs': len(ortholog_df),
                'average_similarity': float(ortholog_df['similarity'].mean()) if len(ortholog_df) > 0 else 0
            }
        )
        
        if not paralog_df.empty:
            registry.register_file(
                'paralog_relationships',
                paralog_df,
                'csv',
                'Paralogous gene relationships detected',
                dependencies=['ortholog_pairs'],
                metadata={'total_paralogs': len(paralog_df)}
            )
        
        # Export standardized formats
        logger.info("Step 2.2: Exporting standardized formats")
        exporter.export_busco_coordinates(ortholog_df)
        
        # Phase 3: Synteny and Inversion Analysis
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 3: SYNTENY AND INVERSION ANALYSIS")
        logger.info("=" * 50)
        
        # Synteny analysis
        logger.info("Step 3.1: Analyzing synteny blocks")
        synteny_df, mapping_df = analyze_enhanced_synteny_blocks(ortholog_df, config)
        
        registry.register_file(
            'synteny_blocks',
            synteny_df,
            'csv',
            'Detected synteny blocks between genomes',
            dependencies=['ortholog_pairs'],
            metadata={
                'total_blocks': len(synteny_df),
                'average_block_size': float(synteny_df['block_size'].mean()) if len(synteny_df) > 0 else 0
            }
        )
        
        registry.register_file(
            'chromosome_mappings',
            mapping_df,
            'csv',
            'Chromosome-to-chromosome mapping statistics',
            dependencies=['ortholog_pairs'],
            metadata={'total_mappings': len(mapping_df)}
        )
        
        # Export synteny in standardized format
        exporter.export_synteny_blocks(synteny_df)
        
        # Rearrangement analysis
        logger.info("Step 3.2: Analyzing chromosome rearrangements")
        rearrangement_df = analyze_enhanced_chromosome_rearrangements(ortholog_df, config)
        
        registry.register_file(
            'chromosome_rearrangements',
            rearrangement_df,
            'csv',
            'Detected chromosome rearrangements (splits/fusions)',
            dependencies=['ortholog_pairs'],
            metadata={'total_rearrangements': len(rearrangement_df)}
        )
        
        # Inversion analysis
        logger.info("Step 3.3: Analyzing inversions")
        inversion_df = analyze_enhanced_inversions(synteny_df, ortholog_df, config)
        
        registry.register_file(
            'inversion_events',
            inversion_df,
            'csv',
            'Detected inversion events',
            dependencies=['synteny_blocks', 'ortholog_pairs'],
            metadata={
                'total_inversions': len(inversion_df),
                'average_inversion_size': float(inversion_df['size_genes'].mean()) if len(inversion_df) > 0 else 0
            }
        )
        
        # Export inversions in standardized format
        if not inversion_df.empty:
            exporter.export_inversion_regions(inversion_df)
        
        # Phase 4: Results Integration and Export
        logger.info("\n" + "=" * 50)
        logger.info("PHASE 4: RESULTS INTEGRATION AND EXPORT")
        logger.info("=" * 50)
        
        # Create comprehensive results dictionary
        results_dict = {
            'ortholog_df': ortholog_df,
            'paralog_df': paralog_df,
            'synteny_df': synteny_df,
            'mapping_df': mapping_df,
            'rearrangement_df': rearrangement_df,
            'inversion_df': inversion_df,
            'first_quality': first_quality,
            'second_quality': second_quality,
            'config': config
        }
        
        # Save legacy results (for compatibility)
        logger.info("Step 4.1: Saving legacy format results")
        save_enhanced_results(output_dir, results_dict, config)
        
        # Export analysis summary
        logger.info("Step 4.2: Exporting analysis summary")
        exporter.export_analysis_summary(results_dict)
        
        # Generate comprehensive visualizations
        logger.info("Step 4.3: Creating enhanced visualizations")
        create_enhanced_visualizations(output_dir, results_dict, config)
        
        # Register visualization outputs
        plots_dir = output_dir / 'plots'
        if plots_dir.exists():
            for plot_file in plots_dir.glob('*.png'):
                registry.register_file(
                    f'plot_{plot_file.stem}',
                    {'path': str(plot_file)},
                    'json',
                    f'Visualization: {plot_file.stem}',
                    dependencies=['ortholog_pairs', 'synteny_blocks', 'inversion_events'],
                    metadata={'plot_type': 'png', 'visualization': True}
                )
        
        # Generate final report with registry integration
        logger.info("Step 4.4: Generating comprehensive report")
        generate_comprehensive_report(output_dir, results_dict, registry)
        
        # Export file manifest
        logger.info("Step 4.5: Exporting file manifest")
        manifest_path = registry.export_manifest()
        
        # Verify registry integrity
        logger.info("Step 4.6: Verifying file integrity")
        integrity_report = {}
        for file_id in registry.list_files():
            integrity_report[file_id] = registry.verify_integrity(file_id)
        
        failed_files = [fid for fid, valid in integrity_report.items() if not valid]
        if failed_files:
            logger.warning(f"Integrity check failed for: {failed_files}")
        else:
            logger.info("✅ All files passed integrity check")
        
        logger.info("\n" + "=" * 80)
        logger.info("ENHANCED ANALYSIS WITH REGISTRY COMPLETED SUCCESSFULLY")
        logger.info("=" * 80)
        
        # Enhanced return with registry information
        return {
            'ortholog_df': ortholog_df,
            'paralog_df': paralog_df,
            'synteny_df': synteny_df,
            'mapping_df': mapping_df,
            'rearrangement_df': rearrangement_df,
            'inversion_df': inversion_df,
            'first_quality': first_quality,
            'second_quality': second_quality,
            'output_dir': output_dir,
            'config': config,
            'registry': registry,
            'file_manifest': str(manifest_path),
            'integrity_report': integrity_report
        }
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        if config.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
        raise


def save_enhanced_results(output_dir, results, config):
    """Legacy results saving function (for compatibility)"""
    data_dir = output_dir / 'data'
    
    # Save main results
    results['ortholog_df'].to_csv(data_dir / Path(config['synteny_analysis_csv']).name, index=False)
    results['inversion_df'].to_csv(data_dir / Path(config['inversion_summary_csv']).name, index=False)
    results['rearrangement_df'].to_csv(data_dir / Path(config['chromosome_rearrangements_csv']).name, index=False)
    
    # Save paralog data if available
    if 'paralog_df' in results and not results['paralog_df'].empty:
        results['paralog_df'].to_csv(data_dir / Path(config['paralog_analysis_csv']).name, index=False)
    
    # Save quality reports
    quality_data = []
    for genome, quality_info in [('first', results['first_quality']), ('second', results['second_quality'])]:
        quality_record = {'genome': genome}
        quality_record.update(quality_info['metrics'])
        quality_record['quality_score'] = quality_info['quality_score']
        quality_record['quality_class'] = quality_info['quality_class']
        quality_data.append(quality_record)
    
    pd.DataFrame(quality_data).to_csv(data_dir / Path(config['quality_report_csv']).name, index=False)
    
    logger.info(f"  Legacy results saved to {data_dir}")


def generate_comprehensive_report(output_dir, results, registry=None):
    """Generate comprehensive analysis report with registry integration"""
    reports_dir = output_dir / 'reports'
    reports_dir.mkdir(exist_ok=True)
    
    # Generate enhanced report with registry information
    if registry:
        report_path = reports_dir / 'analysis_report_with_registry.md'
        
        with open(report_path, 'w') as f:
            f.write("# Genome Inversion Analysis Report\n\n")
            f.write(f"Generated: {pd.Timestamp.now()}\n\n")
            
            # Registry summary
            f.write("## File Registry Summary\n\n")
            f.write(f"- Total registered files: {len(registry.list_files())}\n")
            f.write(f"- Registry location: {registry.registry_file}\n")
            f.write(f"- Manifest: {output_dir}/file_manifest.json\n\n")
            
            # File type breakdown
            f.write("### Files by Type\n\n")
            file_types = {}
            for file_id in registry.list_files():
                file_info = registry.get_file_info(file_id)
                file_type = file_info['type']
                file_types[file_type] = file_types.get(file_type, 0) + 1
            
            for file_type, count in sorted(file_types.items()):
                f.write(f"- {file_type}: {count} files\n")
            
            f.write("\n### Standardized Exports Available\n\n")
            f.write("- BED format: BUSCO coordinates for genome browsers\n")
            f.write("- BEDPE format: Inversion regions for structural variant analysis\n")
            f.write("- JSON format: Machine-readable analysis summaries\n")
            f.write("- GFF format: Feature annotations (if applicable)\n\n")
            
            # Analysis summary
            f.write("## Analysis Results Summary\n\n")
            f.write(f"- Ortholog pairs: {len(results['ortholog_df'])}\n")
            f.write(f"- Synteny blocks: {len(results['synteny_df'])}\n")
            f.write(f"- Inversions detected: {len(results['inversion_df'])}\n")
            f.write(f"- Chromosome rearrangements: {len(results['rearrangement_df'])}\n\n")
            
            # Data provenance
            f.write("## Data Provenance\n\n")
            f.write("This analysis maintains complete data provenance through the registry system:\n\n")
            
            # Show dependency graph for key files
            key_files = ['ortholog_pairs', 'synteny_blocks', 'inversion_events']
            for file_id in key_files:
                if file_id in registry.list_files():
                    deps = registry.get_dependencies(file_id)
                    f.write(f"- **{file_id}**: depends on {deps}\n")
        
        logger.info(f"Enhanced report with registry info: {report_path}")
    
    # Standard report generation
    logger.info(f"Comprehensive report generated in {reports_dir}")


# Update the main execution block to use the new registry system
if __name__ == "__main__":
    # Set random seed for reproducible results
    random.seed(42)
    np.random.seed(42)
    
    # Configuration selection
    if len(sys.argv) > 1:
        if sys.argv[1] == '--fast':
            config = FAST_HYBRID_CONFIG
            logger.info("Starting Fast Analysis with Registry System")
        elif sys.argv[1] == '--hybrid':
            config = ENHANCED_HYBRID_CONFIG
            logger.info("Starting Hybrid Analysis with Registry System")
        elif sys.argv[1] == '--complete':
            config = COMPLETE_ENHANCED_CONFIG
            logger.info("Starting Complete Analysis with Registry System")
        else:
            logger.error(f"Unknown option: {sys.argv[1]}")
            logger.info("Available options: --fast, --hybrid, --complete")
            sys.exit(1)
    else:
        # Default to complete configuration
        config = COMPLETE_ENHANCED_CONFIG
        logger.info("Starting Complete Analysis with Registry System (default)")
        logger.info("Available options: --fast, --hybrid, --complete")
    
    try:
        # Run analysis with registry system
        results = run_complete_enhanced_analysis_with_registry(config)
        
        # Print comprehensive summary with registry info
        print("\n" + "=" * 80)
        config_name = {
            FAST_HYBRID_CONFIG: "FAST",
            ENHANCED_HYBRID_CONFIG: "HYBRID", 
            COMPLETE_ENHANCED_CONFIG: "COMPLETE"
        }.get(config, "UNKNOWN")
        print(f"{config_name} ANALYSIS WITH REGISTRY SYSTEM - SUMMARY")
        print("=" * 80)
        
        # Standard analysis summary (existing code)
        print(f"\nConfiguration Details:")
        strategy = config.get('alignment_strategy', 'unknown')
        print(f"  Alignment strategy: {strategy}")
        
        print(f"\nOrtholog Analysis:")
        print(f"  Total ortholog pairs: {len(results['ortholog_df'])}")
        if len(results['ortholog_df']) > 0:
            print(f"  Average similarity: {results['ortholog_df']['similarity'].mean():.3f}")
            print(f"  Average confidence: {results['ortholog_df']['confidence'].mean():.3f}")
        
        print(f"\nSynteny Analysis:")
        print(f"  Synteny blocks found: {len(results['synteny_df'])}")
        
        print(f"\nInversion Analysis:")
        print(f"  Inversion regions: {len(results['inversion_df'])}")
        
        # Registry system summary
        print(f"\n🗂️ REGISTRY SYSTEM SUMMARY:")
        registry = results['registry']
        print(f"  📁 Total files registered: {len(registry.list_files())}")
        print(f"  📊 File types: {len(set(info['type'] for info in registry.registry['files'].values()))}")
        print(f"  ✅ Integrity check: {'PASSED' if all(results['integrity_report'].values()) else 'FAILED'}")
        print(f"  📋 Manifest: {results['file_manifest']}")
        
        # Standardized exports
        print(f"\n📤 STANDARDIZED EXPORTS:")
        bed_files = registry.list_files('bed')
        json_files = registry.list_files('json')
        print(f"  🧬 BED format files: {len(bed_files)}")
        print(f"  📄 JSON format files: {len(json_files)}")
        
        if bed_files:
            print(f"    → Use with genome browsers (IGV, UCSC, etc.)")
        if json_files:
            print(f"    → Machine-readable for downstream analysis")
        
        print(f"\n📂 OUTPUT STRUCTURE:")
        print(f"  Base directory: {results['output_dir']}")
        print(f"  ├── data/           # CSV files (legacy format)")
        print(f"  ├── exports/        # Standardized formats")
        print(f"  │   ├── bed/        # Genome browser files")
        print(f"  │   ├── json/       # Machine-readable data")
        print(f"  │   └── gff/        # Feature annotations")
        print(f"  ├── plots/          # Visualizations")
        print(f"  ├── cache/          # Analysis cache")
        print(f"  ├── reports/        # Generated reports")
        print(f"  ├── registry.json   # File registry")
        print(f"  └── file_manifest.json  # Complete manifest")
        
        print(f"\n🚀 NEXT STEPS:")
        print(f"  • Load BED files into genome browsers for visualization")
        print(f"  • Use JSON exports for downstream phylogenetic analysis")
        print(f"  • Check file_manifest.json for complete data provenance")
        print(f"  • Registry enables reproducible analysis and data sharing")
        
        print(f"\n🎯 REGISTRY BENEFITS:")
        print(f"  ✅ Complete data provenance tracking")
        print(f"  ✅ Standardized output formats") 
        print(f"  ✅ Dependency management")
        print(f"  ✅ File integrity verification")
        print(f"  ✅ Easy data sharing and collaboration")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        print(f"\nError: {str(e)}")
        if config.get('enable_debug_output', False):
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 80)