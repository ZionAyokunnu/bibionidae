"""
Configuration settings for the Genome Inversion Analyzer
Contains all configuration dictionaries and parameter sets
"""

################################################################################
# HYBRID CONFIGURATION SYSTEM
################################################################################

# Enhanced configuration with hybrid alignment settings
ENHANCED_HYBRID_CONFIG = {
    # ==== INPUT FILES ====
    'first_fasta_path': 'template.fna',
    'second_fasta_path': 'template.fna',
    'first_busco_path': 'template/full_table.tsv',
    'second_busco_path': 'template/full_table.tsv',
    
    # ==== OUTPUT FILES ====
    'base_output_dir': 'enhanced_results',
    'synteny_analysis_csv': 'enhanced_synteny_analysis.csv',
    'inversion_summary_csv': 'enhanced_inversion_summary.csv',
    'chromosome_rearrangements_csv': 'enhanced_chromosome_rearrangements.csv',
    'paralog_analysis_csv': 'paralog_analysis.csv',
    'quality_report_csv': 'assembly_quality_report.csv',
    
    # ==== HYBRID ALIGNMENT CONFIGURATION ====
    'alignment_strategy': 'minimap2',  # 'biopython', 'minimap2', or 'hybrid'
    
    # Sequence length thresholds for alignment method selection
    'short_sequence_threshold': 100,      # BP - use Biopython below this
    'long_sequence_threshold': 500,      # BP - use Minimap2 above this
    'buffer_zone_method': 'dual',         # 'biopython', 'minimap2', or 'dual'
    
    # Minimap2 specific settings
    'minimap2_path': '/Users/za7/Documents/minimap2/minimap2',
    'minimap2_preset': 'sr',            # Short read preset for BUSCO-sized sequences
    'minimap2_kmer_size': 17,             # Smaller k-mer for better sensitivity
    'minimap2_threads': 3,                # CPU cores to use
    'minimap2_min_score': 100,            # Minimum alignment score
    'minimap2_min_identity': 0.7,         # Minimum identity percentage
    'minimap2_extra_flags': '-c --cs',    # Additional flags for detailed output
    
    # Biopython alignment settings (for short sequences)
    'biopython_match_score': 2,
    'biopython_mismatch_score': -1,
    'biopython_gap_open_score': -2,
    'biopython_gap_extend_score': -0.5,
    'biopython_mode': 'local',            # 'local' or 'global'
    'biopython_batch_size': 100,          # Process in batches for progress tracking
    
    # Score normalization and confidence calculation
    'normalize_scores': True,             # Convert all scores to 0-1 scale
    'confidence_weighting': {
        'identity': 0.4,                  # Weight for sequence identity
        'coverage': 0.3,                  # Weight for alignment coverage  
        'length_ratio': 0.2,              # Weight for length similarity
        'score_quality': 0.1              # Weight for raw alignment quality
    },
    
    # Ortholog mapping strategy
    'use_reciprocal_best_hits': True,     # Apply RBH filtering
    'identity_gap_threshold': 0.02,       # 2% - for resolving ties
    'max_paralogs_per_busco': 3,          # Limit paralogs to reduce noise
    'require_bidirectional_match': True,   # Both directions must be best hits
    
    # Quality control and validation
    'validate_alignment_consistency': True, # Check for coordinate mismatches
    'flag_discrepant_alignments': True,   # Mark suspicious alignments for review
    'alignment_length_ratio_min': 0.7,    # Min alignment_length/gene_length
    'cross_validate_buffer_zone': True,   # Run both methods on buffer zone genes
    
    # Performance and debugging
    'enable_parallel_alignment': True,    # Use multiprocessing for Biopython
    'alignment_cache_enabled': True,      # Cache results to avoid recomputation
    'temp_file_cleanup': True,            # Clean up minimap2 temp files
    'detailed_alignment_logging': False,   # Log every alignment (verbose)
    'progress_reporting_interval': 50,    # Report progress every N alignments
    
    # Fallback and error handling
    'fallback_to_simple_similarity': True, # Use difflib if alignment fails
    'skip_failed_alignments': False,      # Whether to skip or retry failed alignments
    'max_alignment_retries': 2,           # Retry failed alignments
    'timeout_per_alignment': 30,          # Seconds before timing out alignment
    
    # ==== STANDARD CONFIGURATION PARAMETERS ====
    'base_similarity_threshold': 0.5,
    'high_quality_similarity_threshold': 0.8,
    'medium_quality_similarity_threshold': 0.6, 
    'low_quality_similarity_threshold': 0.3,
    'fragmented_assembly_similarity_threshold': 0.2,
    
    'base_min_busco_length': 150,
    'high_quality_min_length': 200,
    'medium_quality_min_length': 150,
    'low_quality_min_length': 100,
    'fragmented_assembly_min_length': 50,
    
    'base_min_genes_per_chromosome': 3,
    'base_synteny_correlation_threshold': 0.5,
    'relaxed_correlation_threshold': 0.3,
    'strict_correlation_threshold': 0.8,
    
    'base_min_synteny_block_size': 3,
    'micro_synteny_block_size': 1,
    'base_max_gap_in_synteny': 1000000,
    'adaptive_max_gap_multiplier': 2.0,
    
    'base_min_inversion_size': 2,
    'micro_inversion_size': 1,
    'strand_consistency_threshold': 0.6,
    'inversion_confidence_threshold': 0.7,
    
    # Quality assessment thresholds
    'high_quality_busco_threshold': 0.95,
    'medium_quality_busco_threshold': 0.85,
    'low_quality_busco_threshold': 0.70,
    'high_quality_n50_threshold': 10000000,
    'medium_quality_n50_threshold': 1000000,
    
    # Boolean flags for enhanced features
    'enable_adaptive_thresholds': True,
    'enable_paralog_detection': True,
    'enable_gene_boundary_validation': True,
    'enable_strand_validation': True,
    'enable_translation_check': False,     # Disabled for speed
    'enable_exclusion_warnings': True,
    'enable_chromosome_bounds_check': True,
    'enable_duplicate_handling': True,
    
    'enable_paralog_ortholog_mapping': True,
    'use_dynamic_similarity_threshold': True,
    'enable_gene_model_validation': False, # Disabled for speed
    'enable_ortholog_confidence_scoring': True,
    'enable_many_to_many_mapping': False,  # Use RBH instead
    
    'enable_small_synteny_blocks': True,
    'enable_translocation_detection': True,
    'enable_synteny_confidence_scoring': True,
    
    'enable_single_gene_inversions': True,
    'enable_micro_inversions': True,
    'use_content_based_inversion': False,  # Disabled for speed
    'enable_inversion_confidence': True,
    
    'enable_assembly_quality_assessment': True,
    'enable_statistical_validation': False, # Disabled for speed
    'enable_debug_output': True,
    
    # Plotting parameters
    'plot_width': 15,
    'plot_height': 10,
    'dpi': 300,
    'color_palette': 'viridis',
    'font_size': 12,
    'figure_format': 'png',
    
    # Visualization parameters
    'generate_dotplots': True,            # Generate synteny dotplots
    'dotplot_show_labels': False,         # Show gene labels on plots
    'dotplot_by': 'busco',               # 'gene', 'nucleotide', 'busco'
    'dotplot_size': (12, 10),            # Figure size for dotplots
    'synteny_color': '#1f77b4',          # Blue for syntenic regions
    'inversion_color': '#d62728',        # Red for inverted regions
    'confidence_alpha': True,            # Use alpha transparency for confidence
    'show_synteny_blocks': True,         # Overlay synteny block lines
    'show_breakpoints': True,            # Mark rearrangement breakpoints
    
    # Reporting parameters
    'report_format': 'txt',              # 'txt', 'md', 'html'
    'include_detailed_stats': True,      # Include detailed statistics
    'include_methodology': True,         # Include methods section
    'include_quality_metrics': True,     # Include assembly quality details
    'generate_summary_plots': True,      # Generate summary visualizations
}

# Fast configuration for quick testing (Biopython only, minimal features)
FAST_HYBRID_CONFIG = {
    **ENHANCED_HYBRID_CONFIG,  # Inherit all settings
    
    # Override for speed
    'alignment_strategy': 'biopython',     # Skip minimap2 for simplicity
    'short_sequence_threshold': 99999,     # Force all sequences through Biopython
    'enable_parallel_alignment': True,     # Use multiprocessing
    'biopython_batch_size': 50,           # Smaller batches
    'detailed_alignment_logging': False,
    'progress_reporting_interval': 25,
    
    # Simplified features
    'use_reciprocal_best_hits': False,    # Simple best hit
    'validate_alignment_consistency': False,
    'cross_validate_buffer_zone': False,
    'enable_translation_check': False,
    'enable_gene_model_validation': False,
    'use_content_based_inversion': False,
    'enable_statistical_validation': False,
    
    # Relaxed thresholds for speed
    'base_similarity_threshold': 0.4,
    'alignment_length_ratio_min': 0.6,
    'timeout_per_alignment': 15,
    
    'base_output_dir': 'v4/fast_results',
    'synteny_analysis_csv': 'v4/fast_synteny_analysis.csv',
    'inversion_summary_csv': 'v4/fast_inversion_summary.csv',
    'chromosome_rearrangements_csv': 'v4/fast_chromosome_rearrangements.csv',
    'quality_report_csv': 'v4/fast_quality_report.csv'
}

# Complete configuration (original, all features)
COMPLETE_ENHANCED_CONFIG = {
    **ENHANCED_HYBRID_CONFIG,
    # CHANGE THIS LINE - use 'hybrid' instead of 'biopython'
    'alignment_strategy': 'hybrid',     # ← This was your problem!
    
    # Keep all the enhanced features
    'enable_translation_check': True,
    'enable_gene_model_validation': True,
    'use_content_based_inversion': True,
    'enable_statistical_validation': True,
    'enable_assembly_quality_assessment': True,
}

PUBLICATION_CONFIG = {
    # Synteny visualization
    'use_synteny_plotter': True,
    'create_curved_ribbons': True,
    'create_chord_diagrams': True,
    
    # Phylogenetic integration  
    'use_existing_diptera_tree': True,
    'diptera_tree_path': '/path/to/sanger_diptera.nwk',
    'prune_to_species': True,
    'annotate_with_inversions': True,
    
    # Annotation metrics
    'inversion_metrics': ['raw_count', 'normalized_rate', 'per_mb'],
    'edge_annotation_style': 'heatmap',  # 'heatmap', 'width', 'color'
    
    # Output control
    'create_publication_plots': True,
    'plot_formats': ['png', 'pdf', 'svg'],
    'high_dpi': 300,
    
    # Charlotte's R-based synteny plotter settings (WORKING VERSION)
    'external_tools': {
        'synteny_plotter': 'genome_inversion_analyser/external_tools/synteny_plotter'  # Directory path
    },
    
    # Synteny visualization settings
    'synteny_visualization': {
        'enabled': True,
        'create_curved_ribbons': True,
        'create_straight_links': False,
        'options': {
            'filter_threshold': 5,  # Minimum BUSCOs per chromosome (R script -f parameter)
            'gap': 6,              # Gap between chromosomal sets (R script -g parameter)
            'alpha': 10,           # Transparency percentage (R script -alpha parameter)
            'output_format': 'png'
        }
    },
    
    # Phylogenetic tree annotation
    'tree_annotation': {
        'enabled': True,
        'source_tree_path': 'diptera_clean_20species.newick',
        'prune_to_target_species': True,
        'annotation_metrics': {
            'inversion_count': True,
            'inversion_rate_per_mb': True,
            'normalized_inversion_score': True
        },
        'visualization': {
            'edge_annotation_style': 'heatmap',  # 'heatmap', 'width', 'color', 'labels'
            'node_support_values': True,
            'output_formats': ['png', 'pdf', 'newick'],
            'tree_layout': 'rectangular',  # 'rectangular', 'circular'
            'dpi': 300
        }
    },
    
    # Publication suite control
    'publication_suite': {
        'enabled': True,
        'create_all_plots': True,
        'high_quality_output': True,
        'create_supplementary_data': True
    }
}