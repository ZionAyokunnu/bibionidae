# =============================================================================
# Analysis Workflow Integration (workflow.py)
# =============================================================================

"""
Main workflow integration for the analysis system.
Orchestrates synteny, inversion, rearrangement, and quality analysis.
"""

import pandas as pd
from typing import Dict, Tuple, Any

# Additional imports for enhanced/hybrid workflow
from ..config import ENHANCED_HYBRID_CONFIG
from ..utils import create_output_directory, generate_cache_key
from ..io_utils.busco_parser import BuscoParser
from ..io_utils.sequence_extractor import SequenceExtractor
from ..alignment.workflow import run_hybrid_alignment_analysis
from .synteny_analyzer import SyntenyAnalyzer
from .inversion_detector import InversionDetector
from .rearrangement_analyzer import RearrangementAnalyzer
from .quality_assessor import QualityAssessor
from .statistical_validator import StatisticalValidator
from ..visualization.dashboard import create_enhanced_visualizations
from ..reporting.report_generator import generate_comprehensive_report
from ..logger import get_logger

logger = get_logger()

class AnalysisWorkflow:
    """
    Main workflow for genomic analysis.
    Orchestrates synteny, inversion, rearrangement, and quality analysis.
    """
    
    def __init__(self, config):
        """
        Initialize analysis workflow.
        
        Args:
            config: Configuration object
        """
        self.config = config
        
        # Initialize analysis components
        self.synteny_analyzer = SyntenyAnalyzer(config)
        self.inversion_detector = InversionDetector(config)
        self.rearrangement_analyzer = RearrangementAnalyzer(config)
        self.quality_assessor = QualityAssessor(config)
        self.statistical_validator = StatisticalValidator(config)
        
        logger.info("Analysis workflow initialized")
    
    def run_complete_analysis(self, ortholog_df: pd.DataFrame,
                            first_fasta_path: str = None,
                            second_fasta_path: str = None,
                            first_busco_df: pd.DataFrame = None,
                            second_busco_df: pd.DataFrame = None) -> Dict[str, Any]:
        """
        Run complete genomic analysis pipeline.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            first_fasta_path: Optional path to first genome FASTA
            second_fasta_path: Optional path to second genome FASTA
            first_busco_df: Optional first genome BUSCO data
            second_busco_df: Optional second genome BUSCO data
            
        Returns:
            Dictionary with all analysis results
        """
        logger.info("Running complete genomic analysis pipeline...")
        
        results = {}
        
        # Step 1: Quality assessment (if FASTA and BUSCO data provided)
        if first_fasta_path and first_busco_df is not None:
            logger.info("Step 1: Assessing first genome quality...")
            results['first_quality'] = self.quality_assessor.assess_assembly_quality(
                first_fasta_path, first_busco_df
            )
        
        if second_fasta_path and second_busco_df is not None:
            logger.info("Step 1: Assessing second genome quality...")
            results['second_quality'] = self.quality_assessor.assess_assembly_quality(
                second_fasta_path, second_busco_df
            )
        
        # Quality comparison if both available
        if 'first_quality' in results and 'second_quality' in results:
            results['quality_comparison'] = self.quality_assessor.compare_assembly_qualities(
                results['first_quality'], results['second_quality']
            )
        
        # Step 2: Synteny analysis
        logger.info("Step 2: Analyzing synteny blocks...")
        synteny_df, mapping_df = self.synteny_analyzer.analyze_synteny_blocks(ortholog_df)
        results['synteny_df'] = synteny_df
        results['mapping_df'] = mapping_df
        
        # Step 3: Inversion detection
        logger.info("Step 3: Detecting inversions...")
        inversion_df = self.inversion_detector.detect_inversions(synteny_df, ortholog_df)
        results['inversion_df'] = inversion_df
        
        # Step 4: Rearrangement analysis
        logger.info("Step 4: Analyzing chromosome rearrangements...")
        rearrangement_df = self.rearrangement_analyzer.analyze_chromosome_rearrangements(ortholog_df)
        results['rearrangement_df'] = rearrangement_df
        
        # Step 5: Statistical validation (if enabled)
        if self.config.get('enable_statistical_validation', False):
            logger.info("Step 5: Performing statistical validation...")
            results['synteny_validation'] = self.statistical_validator.validate_synteny_results(
                results['synteny_df'], ortholog_df
            )
        return results

    def run_complete_enhanced_analysis_with_hybrid(self, config=None):
        """
        Run the complete enhanced analysis workflow using the hybrid alignment strategy.
        This orchestrates BUSCO parsing, sequence extraction, hybrid alignment, synteny/inversion/rearrangement analysis,
        quality assessment, statistical validation, visualization, and reporting.
        Args:
            config: Optional configuration object. If None, uses ENHANCED_HYBRID_CONFIG.
        Returns:
            Dictionary with all analysis results and generated outputs.
        """
        # Use default config if not provided
        if config is None:
            config = ENHANCED_HYBRID_CONFIG

        # Step 1: Prepare output directory and cache key
        output_dir = config.get('base_output_dir', 'results')
        create_output_directory(output_dir)
        cache_key = generate_cache_key(config)

        # Step 2: Parse BUSCO tables
        logger.info("Parsing BUSCO tables...")
        first_busco_parser = BuscoParser(config.get('first_busco_path'))
        first_busco_df = first_busco_parser.parse_busco_table(config)
        second_busco_parser = BuscoParser(config.get('second_busco_path'))
        second_busco_df = second_busco_parser.parse_busco_table(config)

        # Step 3: Extract BUSCO gene sequences
        logger.info("Extracting BUSCO gene sequences...")
        seq_extractor_1 = SequenceExtractor(config.get('first_fasta_path'))
        busco_seqs_1 = seq_extractor_1.extract_busco_sequences(first_busco_df, config)
        seq_extractor_2 = SequenceExtractor(config.get('second_fasta_path'))
        busco_seqs_2 = seq_extractor_2.extract_busco_sequences(second_busco_df, config)

        # Step 4: Run hybrid alignment analysis
        logger.info("Running hybrid alignment analysis...")
        alignment_results, sequence_pairs = run_hybrid_alignment_analysis(
            busco_seqs_1, busco_seqs_2, config
        )

        # Step 5: Orthology processing (score, filter, convert)
        logger.info("Processing orthologs...")
        from ..orthology.workflow import OrthologyWorkflow
        orthology_workflow = OrthologyWorkflow(config)
        ortholog_df, paralog_df = orthology_workflow.process_alignment_results(
            alignment_results, sequence_pairs
        )

        # Step 6: Synteny analysis
        logger.info("Analyzing synteny blocks...")
        synteny_df, mapping_df = self.synteny_analyzer.analyze_synteny_blocks(ortholog_df)

        # Step 7: Inversion detection
        logger.info("Detecting inversions...")
        inversion_df = self.inversion_detector.detect_inversions(synteny_df, ortholog_df)

        # Step 8: Rearrangement analysis
        logger.info("Analyzing chromosome rearrangements...")
        rearrangement_df = self.rearrangement_analyzer.analyze_chromosome_rearrangements(ortholog_df)

        # Step 9: Quality assessment
        logger.info("Assessing assembly quality...")
        first_quality = self.quality_assessor.assess_assembly_quality(
            config.get('first_fasta_path'), first_busco_df
        )
        second_quality = self.quality_assessor.assess_assembly_quality(
            config.get('second_fasta_path'), second_busco_df
        )
        quality_comparison = self.quality_assessor.compare_assembly_qualities(
            first_quality, second_quality
        )

        # Step 10: Statistical validation (if enabled)
        synteny_validation = None
        inversion_validation = None
        if self.config.get('enable_statistical_validation', False):
            logger.info("Performing statistical validation...")
            synteny_validation = self.statistical_validator.validate_synteny_results(
                synteny_df, ortholog_df
            )
            inversion_validation = self.statistical_validator.validate_inversion_results(
                inversion_df
            )

        # Step 11: Visualization
        logger.info("Generating visualizations...")
        create_enhanced_visualizations(
            ortholog_df=ortholog_df,
            synteny_df=synteny_df,
            inversion_df=inversion_df,
            rearrangement_df=rearrangement_df,
            quality_report={'first': first_quality, 'second': second_quality, 'comparison': quality_comparison},
            output_dir=output_dir,
            config=config
        )

        # Step 12: Reporting
        logger.info("Generating comprehensive report...")
        report_path = generate_comprehensive_report(
            ortholog_df=ortholog_df,
            synteny_df=synteny_df,
            inversion_df=inversion_df,
            rearrangement_df=rearrangement_df,
            paralog_df=paralog_df,
            quality_report={'first': first_quality, 'second': second_quality, 'comparison': quality_comparison},
            synteny_validation=synteny_validation,
            inversion_validation=inversion_validation,
            output_dir=output_dir,
            config=config
        )

        # Collate and return all results
        results = {
            'ortholog_df': ortholog_df,
            'paralog_df': paralog_df,
            'synteny_df': synteny_df,
            'mapping_df': mapping_df,
            'inversion_df': inversion_df,
            'rearrangement_df': rearrangement_df,
            'first_quality': first_quality,
            'second_quality': second_quality,
            'quality_comparison': quality_comparison,
            'synteny_validation': synteny_validation,
            'inversion_validation': inversion_validation,
            'visualizations': output_dir,
            'report_path': report_path
        }
        logger.info("Enhanced hybrid analysis complete.")
        return results