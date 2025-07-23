# =============================================================================
# Analysis Workflow Integration (workflow.py)
# =============================================================================

"""
Main workflow integration for the analysis system.
Orchestrates synteny, inversion, rearrangement, and quality analysis.
"""

import pandas as pd
from typing import Dict, Tuple, Any

from .synteny_analyzer import SyntenyAnalyzer
from .inversion_detector import InversionDetector
from .rearrangement_analyzer import RearrangementAnalyzer
from .quality_assessor import QualityAssessor
from .statistical_validator import StatisticalValidator
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
                synteny_df, ortholog_df
            )
            results['inversion_validation'] = self.statistical_validator.validate_inversion_results(
                inversion_df
            )
        
        # Step 6: Generate comprehensive statistics
        logger.info("Step 6: Generating comprehensive statistics...")
        results['statistics'] = self._generate_comprehensive_statistics(results)
        
        # Log final summary
        self._log_analysis_summary(results)
        
        return results
    
    def _generate_comprehensive_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate comprehensive statistics for all analysis components."""
        stats = {}
        
        # Synteny statistics
        if 'synteny_df' in results and 'mapping_df' in results:
            stats['synteny'] = self.synteny_analyzer.get_synteny_statistics(
                results['synteny_df'], results['mapping_df']
            )
        
        # Inversion statistics
        if 'inversion_df' in results:
            stats['inversions'] = self.inversion_detector.get_inversion_statistics(
                results['inversion_df']
            )
        
        # Rearrangement statistics
        if 'rearrangement_df' in results:
            stats['rearrangements'] = self.rearrangement_analyzer.get_rearrangement_statistics(
                results['rearrangement_df']
            )
        
        # Quality statistics
        if 'first_quality' in results:
            stats['first_genome_quality'] = self.quality_assessor.get_quality_statistics(
                results['first_quality']
            )
        
        if 'second_quality' in results:
            stats['second_genome_quality'] = self.quality_assessor.get_quality_statistics(
                results['second_quality']
            )
        
        # Overall summary statistics
        stats['summary'] = self._calculate_summary_statistics(results)
        
        return stats
    
    def _calculate_summary_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Calculate high-level summary statistics."""
        summary = {}
        
        # Count major findings
        summary['total_synteny_blocks'] = len(results.get('synteny_df', []))
        summary['total_inversions'] = len(results.get('inversion_df', []))
        summary['total_rearrangements'] = len(results.get('rearrangement_df', []))
        
        # Calculate structural variation rates - need to pass ortholog count
        # Try to infer ortholog count from synteny blocks or passed data
        total_orthologs = 0
        if 'synteny_df' in results and len(results['synteny_df']) > 0:
            total_orthologs = results['synteny_df']['block_size'].sum()
        elif 'mapping_df' in results and len(results['mapping_df']) > 0:
            total_orthologs = results['mapping_df']['gene_count'].sum()
        
        if total_orthologs > 0:
            summary['inversion_rate'] = summary['total_inversions'] / total_orthologs
            summary['rearrangement_rate'] = summary['total_rearrangements'] / total_orthologs
        
        # Quality indicators
        if 'first_quality' in results and 'second_quality' in results:
            summary['average_quality_score'] = (
                results['first_quality']['quality_score'] + 
                results['second_quality']['quality_score']
            ) / 2
            
            summary['quality_classes'] = [
                results['first_quality']['quality_class'],
                results['second_quality']['quality_class']
            ]
        
        # Validation status
        if 'synteny_validation' in results:
            summary['synteny_validated'] = results['synteny_validation'].get('overall_validation', {}).get('validated', False)
        
        if 'inversion_validation' in results:
            summary['inversions_validated'] = results['inversion_validation'].get('overall_validation', {}).get('validated', False)
        
        return summary
    
    def _log_analysis_summary(self, results: Dict[str, Any]):
        """Log comprehensive analysis summary."""
        logger.info("Analysis pipeline completed successfully!")
        
        # Log major findings
        logger.info("Major findings:")
        logger.info(f"  Synteny blocks: {len(results.get('synteny_df', []))}")
        logger.info(f"  Inversions detected: {len(results.get('inversion_df', []))}")
        logger.info(f"  Rearrangements detected: {len(results.get('rearrangement_df', []))}")
        
        # Log quality assessment
        if 'first_quality' in results and 'second_quality' in results:
            q1 = results['first_quality']
            q2 = results['second_quality']
            logger.info("Assembly quality:")
            logger.info(f"  First genome: {q1['quality_class']} (score: {q1['quality_score']:.3f})")
            logger.info(f"  Second genome: {q2['quality_class']} (score: {q2['quality_score']:.3f})")
        
        # Log validation status
        validation_status = []
        if 'synteny_validation' in results:
            synteny_valid = results['synteny_validation'].get('overall_validation', {}).get('validated', False)
            validation_status.append(f"Synteny: {'✓' if synteny_valid else '✗'}")
        
        if 'inversion_validation' in results:
            inversion_valid = results['inversion_validation'].get('overall_validation', {}).get('validated', False)
            validation_status.append(f"Inversions: {'✓' if inversion_valid else '✗'}")
        
        if validation_status:
            logger.info(f"Statistical validation: {', '.join(validation_status)}")
    
    def get_analysis_configuration(self) -> Dict[str, Any]:
        """Get current analysis configuration for reporting."""
        return {
            'synteny_parameters': {
                'min_genes_per_block': getattr(self.synteny_analyzer, 'min_genes_per_block', 3),
                'min_synteny_block_size': getattr(self.synteny_analyzer, 'min_synteny_block_size', 3),
                'correlation_threshold': getattr(self.synteny_analyzer, 'correlation_threshold', 0.5),
                'strand_consistency_threshold': getattr(self.synteny_analyzer, 'strand_consistency_threshold', 0.6)
            },
            'inversion_parameters': {
                'min_inversion_size': self.inversion_detector.min_inversion_size,
                'enable_single_gene': self.inversion_detector.enable_single_gene,
                'enable_micro_inversions': self.inversion_detector.enable_micro_inversions,
                'confidence_threshold': self.inversion_detector.confidence_threshold
            },
            'rearrangement_parameters': {
                'enable_translocation_detection': self.rearrangement_analyzer.enable_translocation_detection,
                'min_genes_for_rearrangement': self.rearrangement_analyzer.min_genes_for_rearrangement,
                'rearrangement_confidence_threshold': self.rearrangement_analyzer.rearrangement_confidence_threshold
            },
            'quality_parameters': {
                'enable_assessment': self.quality_assessor.enable_assessment,
                'enable_adaptive_thresholds': self.quality_assessor.enable_adaptive_thresholds,
                'high_quality_busco_threshold': self.quality_assessor.high_quality_busco_threshold,
                'high_quality_n50_threshold': self.quality_assessor.high_quality_n50_threshold
            },
            'validation_parameters': {
                'enable_validation': self.statistical_validator.enable_validation,
                'confidence_level': self.statistical_validator.confidence_level,
                'bootstrap_iterations': self.statistical_validator.bootstrap_iterations,
                'significance_threshold': self.statistical_validator.significance_threshold
            }
        }