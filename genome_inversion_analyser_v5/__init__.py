# =============================================================================
# Genome Inversion Analyser Package (__init__.py)
# =============================================================================

"""
Enhanced Genome Synteny and Inversion Analysis Package

A comprehensive toolkit for analyzing chromosomal rearrangements, synteny blocks,
and inversions between genomes using BUSCO genes as anchors.

Features:
- Hybrid alignment system (Minimap2 + Biopython)
- Assembly quality assessment  
- Reciprocal best hit filtering
- Parallel processing and caching
- Comprehensive visualization
- Publication-ready reports
"""

__version__ = "2.0.0"
__author__ = "Genome Analysis Team"
__description__ = "Comprehensive genome synteny and inversion analysis toolkit"

# Import main configuration classes
from .config import (
    AnalysisConfig,
    ENHANCED_HYBRID_CONFIG,
    FAST_HYBRID_CONFIG, 
    COMPLETE_ENHANCED_CONFIG
)

# Import logging utilities
from .logger import setup_logger, get_logger

# Import main analysis modules
from .analysis import (
    SyntenyAnalyzer,
    InversionDetector,
    RearrangementAnalyzer,
    QualityAssessor,
    StatisticalValidator
)

# Import quality assessment
# from .quality import AssemblyQualityAnalyzer

# Import visualization components
from .visualization import (
    SyntenyDotplot,
    SummaryPlotter,
    QualityPlotter,
    AnalysisDashboard
)

# Import reporting
from .reporting import ReportGenerator, StatisticsFormatter

# Main analysis function for easy use
def analyze_genomes(first_fasta, second_fasta, first_busco, second_busco, 
                   config=None, output_dir='results'):
    """
    Main convenience function for genome analysis.
    
    Args:
        first_fasta: Path to first genome FASTA file
        second_fasta: Path to second genome FASTA file  
        first_busco: Path to first genome BUSCO table
        second_busco: Path to second genome BUSCO table
        config: Optional configuration (uses ENHANCED_HYBRID_CONFIG if None)
        output_dir: Output directory path
        
    Returns:
        Dictionary with analysis results
    """
    # Use default config if none provided
    if config is None:
        config_dict = ENHANCED_HYBRID_CONFIG.to_dict()
        config_dict.update({
            'first_fasta_path': first_fasta,
            'second_fasta_path': second_fasta,
            'first_busco_path': first_busco,
            'second_busco_path': second_busco,
            'base_output_dir': output_dir
        })
        config = AnalysisConfig(config_dict)
    
    # Setup logger
    logger = setup_logger(log_file=f"{output_dir}/analysis.log")
    logger.section_header("Genome Inversion Analysis Pipeline v2.0")
    
    # Import main analysis function (to avoid circular imports)
    # from .main_analysis import run_enhanced_analysis
    
    # Run analysis
    # results = run_enhanced_analysis(config)
    
    logger.info("Analysis completed successfully!")
    return results

# Package metadata
__all__ = [
    # Configuration
    'AnalysisConfig',
    'ENHANCED_HYBRID_CONFIG',
    'FAST_HYBRID_CONFIG',
    'COMPLETE_ENHANCED_CONFIG',
    
    # Logging
    'setup_logger',
    'get_logger',
    
    # Analysis modules
    'SyntenyAnalyzer',
    'InversionDetector',
    'RearrangementAnalyzer',
    'QualityAssessor',
    'StatisticalValidator',
    
    # Quality assessment
    # 'AssemblyQualityAnalyzer',
    
    # Visualization
    'SyntenyDotplot',
    'SummaryPlotter',
    'QualityPlotter',
    'AnalysisDashboard',
    
    # Reporting
    'ReportGenerator',
    'StatisticsFormatter',
    
    # Main function
    'analyze_genomes'
]