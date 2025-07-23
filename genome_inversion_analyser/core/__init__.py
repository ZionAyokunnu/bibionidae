"""Core analysis modules for genome inversion analyzer"""

from .quality_assessment import (
    assess_assembly_quality,
    calculate_quality_score,
    classify_assembly_quality,
    suggest_parameter_adjustments
)

from .busco_processor import (
    enhanced_parse_busco_table,
    detect_and_annotate_paralogs,
    enhanced_filter_busco_genes,
    validate_gene_boundaries,
    validate_translation,
    extract_enhanced_busco_sequences
)

from .alignment import (
    partition_sequences_by_length,
    create_minimap2_fasta,
    run_minimap2_alignment,
    parse_minimap2_paf,
    run_biopython_alignment_batch,
    run_simple_biopython_alignment,
    run_parallel_biopython_alignment,
    normalize_alignment_scores,
    apply_reciprocal_best_hit_filtering,
    convert_alignment_results_to_ortholog_pairs,
    run_hybrid_alignment_analysis,
    check_minimap2_available,
    select_best_buffer_results,
    setup_hybrid_sequence_aligner
)

from .synteny_analyzer import (
    analyze_enhanced_synteny_blocks,
    classify_enhanced_synteny_type,
    calculate_synteny_confidence,
    analyze_enhanced_chromosome_rearrangements,
    calculate_rearrangement_confidence,
    analyze_enhanced_inversions,
    identify_enhanced_inversions,
    calculate_inversion_confidence,
    detect_single_gene_inversions
)

__all__ = [
    # Quality assessment
    'assess_assembly_quality',
    'calculate_quality_score', 
    'classify_assembly_quality',
    'suggest_parameter_adjustments',
    # BUSCO processing
    'enhanced_parse_busco_table',
    'detect_and_annotate_paralogs',
    'enhanced_filter_busco_genes',
    'validate_gene_boundaries',
    'validate_translation',
    'extract_enhanced_busco_sequences',
    # Alignment system
    'partition_sequences_by_length',
    'create_minimap2_fasta',
    'run_minimap2_alignment',
    'parse_minimap2_paf',
    'run_biopython_alignment_batch',
    'run_simple_biopython_alignment',
    'run_parallel_biopython_alignment',
    'normalize_alignment_scores',
    'apply_reciprocal_best_hit_filtering',
    'convert_alignment_results_to_ortholog_pairs',
    'run_hybrid_alignment_analysis',
    'check_minimap2_available',
    'select_best_buffer_results',
    'setup_hybrid_sequence_aligner',
    # Synteny analysis
    'analyze_enhanced_synteny_blocks',
    'classify_enhanced_synteny_type',
    'calculate_synteny_confidence',
    'analyze_enhanced_chromosome_rearrangements',
    'calculate_rearrangement_confidence',
    'analyze_enhanced_inversions',
    'identify_enhanced_inversions',
    'calculate_inversion_confidence',
    'detect_single_gene_inversions'
]