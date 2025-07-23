
# =============================================================================
# 7. Orthology Workflow Integration (orthology/workflow.py)
# =============================================================================

"""
Main workflow integration for the orthology processing system.
Orchestrates the complete pipeline from alignment results to ortholog identification.
"""

import pandas as pd
from typing import List, Tuple, Dict, Any

from ..alignment.alignment_result import AlignmentResult
from .scorer import AlignmentScorer
from .rbh_filter import ReciprocalBestHitFilter
from .converter import OrthologConverter
from .paralog_handler import ParalogHandler
from .confidence_calculator import ConfidenceCalculator
from ..logger import get_logger

logger = get_logger()

class OrthologyWorkflow:
    """
    Main workflow for orthology processing.
    Orchestrates scoring, filtering, conversion, and paralog analysis.
    """
    
    def __init__(self, config):
        """
        Initialize orthology workflow.
        
        Args:
            config: Configuration object
        """
        self.config = config
        
        # Initialize components
        self.scorer = AlignmentScorer(config)
        self.rbh_filter = ReciprocalBestHitFilter(config)
        self.converter = OrthologConverter(config)
        self.paralog_handler = ParalogHandler(config)
        self.confidence_calculator = ConfidenceCalculator(config)
        
        logger.info("Orthology workflow initialized")
    
    def process_alignment_results(self, alignment_results: List[AlignmentResult],
                                sequence_pairs: List[Dict]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process alignment results through the complete orthology pipeline.
        
        Args:
            alignment_results: List of AlignmentResult objects
            sequence_pairs: Original sequence pairs with metadata
            
        Returns:
            Tuple of (ortholog_df, paralog_df)
        """
        logger.info("Processing alignment results through orthology pipeline...")
        
        # Step 1: Score normalization and confidence calculation
        logger.info("  Step 1: Normalizing scores and calculating confidence...")
        normalized_results = self.scorer.normalize_alignment_scores(alignment_results)
        
        # Step 2: Enhanced confidence calculation
        logger.info("  Step 2: Calculating enhanced confidence scores...")
        enhanced_results = self._enhance_confidence_scores(normalized_results)
        
        # Step 3: Reciprocal best hit filtering
        logger.info("  Step 3: Applying reciprocal best hit filtering...")
        filtered_results = self.rbh_filter.apply_rbh_filtering(enhanced_results)
        
        # Step 4: Convert to ortholog pairs
        logger.info("  Step 4: Converting to ortholog pairs...")
        ortholog_df = self.converter.convert_to_ortholog_pairs(filtered_results, sequence_pairs)
        
        # Step 5: Paralog analysis
        logger.info("  Step 5: Analyzing paralog relationships...")
        paralog_df = self.paralog_handler.analyze_ortholog_paralog_relationships(ortholog_df)
        
        # Step 6: Final validation and statistics
        logger.info("  Step 6: Final validation and statistics...")
        self._validate_and_report_results(ortholog_df, paralog_df, alignment_results, filtered_results)
        
        return ortholog_df, paralog_df
    
    def _enhance_confidence_scores(self, results: List[AlignmentResult]) -> List[AlignmentResult]:
        """Enhance confidence scores using advanced calculation."""
        enhanced_results = []
        
        for result in results:
            # Calculate enhanced confidence
            enhanced_confidence = self.confidence_calculator.calculate_ortholog_confidence(result)
            
            # Update result
            result.confidence = enhanced_confidence
            enhanced_results.append(result)
        
        return enhanced_results
    
    def _validate_and_report_results(self, ortholog_df: pd.DataFrame, paralog_df: pd.DataFrame,
                                   original_results: List[AlignmentResult], 
                                   filtered_results: List[AlignmentResult]):
        """Validate results and generate comprehensive statistics."""
        
        # Basic validation
        if len(ortholog_df) == 0:
            logger.warning("  No ortholog pairs identified - check thresholds and input quality")
        
        # Generate comprehensive statistics
        logger.info("  Final orthology processing statistics:")
        logger.info(f"    Original alignment results: {len(original_results)}")
        logger.info(f"    After RBH filtering: {len(filtered_results)}")
        logger.info(f"    Final ortholog pairs: {len(ortholog_df)}")
        logger.info(f"    Complex paralog relationships: {len(paralog_df)}")
        
        if len(ortholog_df) > 0:
            # Quality metrics
            avg_similarity = ortholog_df['similarity'].mean()
            avg_confidence = ortholog_df['confidence'].mean()
            avg_identity = ortholog_df['identity'].mean()
            
            logger.info(f"    Average similarity: {avg_similarity:.3f}")
            logger.info(f"    Average confidence: {avg_confidence:.3f}")
            logger.info(f"    Average identity: {avg_identity:.3f}")
            
            # Coverage statistics
            unique_buscos = ortholog_df['busco_id'].nunique()
            chromosome_pairs = len(ortholog_df.groupby(['first_chr', 'second_chr']))
            
            logger.info(f"    Unique BUSCOs covered: {unique_buscos}")
            logger.info(f"    Chromosome pairs involved: {chromosome_pairs}")
            
            # Method distribution
            if 'alignment_method' in ortholog_df.columns:
                method_dist = ortholog_df['alignment_method'].value_counts()
                logger.info(f"    Final method distribution: {method_dist.to_dict()}")
    
    def get_comprehensive_statistics(self, ortholog_df: pd.DataFrame = None, 
                                   paralog_df: pd.DataFrame = None,
                                   alignment_results: List[AlignmentResult] = None) -> Dict[str, Any]:
        """
        Get comprehensive statistics for the entire orthology processing pipeline.
        
        Args:
            ortholog_df: Final ortholog DataFrame
            paralog_df: Final paralog DataFrame  
            alignment_results: Original alignment results
            
        Returns:
            Dictionary with comprehensive statistics
        """
        stats = {}
        
        # Scoring statistics
        if alignment_results:
            stats['scoring'] = self.scorer.get_score_statistics(alignment_results)
        
        # Filtering statistics
        if alignment_results and ortholog_df is not None:
            # Approximate filtered results count from ortholog count
            filtered_count = len(ortholog_df)  # This is an approximation
            stats['filtering'] = self.rbh_filter.get_filtering_statistics(
                alignment_results, alignment_results[:filtered_count]
            )
        
        # Conversion statistics
        if ortholog_df is not None:
            stats['conversion'] = self.converter.get_conversion_statistics(ortholog_df)
        
        # Paralog statistics
        if paralog_df is not None:
            stats['paralogs'] = self.paralog_handler.get_paralog_statistics(
                ortholog_df=ortholog_df, paralog_df=paralog_df
            )
        
        # Confidence statistics
        if alignment_results:
            stats['confidence'] = self.confidence_calculator.get_confidence_statistics(alignment_results)
        
        return stats