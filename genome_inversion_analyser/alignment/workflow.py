
# =============================================================================
# Main Workflow Integration (alignment/workflow.py)
# =============================================================================

"""
Main workflow integration for the alignment system.
Provides high-level interface for the alignment pipeline.
"""

from typing import Tuple, Dict
import pandas as pd

from ..logger import get_logger
from .hybrid_alignment import HybridAlignmentEngine

logger = get_logger()

class AlignmentWorkflow:
    """
    High-level workflow interface for the alignment system.
    Simplifies integration with the main analysis pipeline.
    """
    
    def __init__(self, config):
        """
        Initialize alignment workflow.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.engine = HybridAlignmentEngine(config)
        
        logger.info("Alignment workflow initialized")
    
    def run_alignment_analysis(self, first_busco_df: pd.DataFrame, 
                             second_busco_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run complete alignment analysis.
        
        Args:
            first_busco_df: First genome BUSCO data with sequences
            second_busco_df: Second genome BUSCO data with sequences
            
        Returns:
            Tuple of (ortholog_df, paralog_df)
        """
        logger.info("Starting alignment analysis workflow...")
        
        # Validate input data
        self._validate_input_data(first_busco_df, second_busco_df)
        
        # Run alignment pipeline
        ortholog_df, paralog_df = self.engine.run_hybrid_alignment_analysis(
            first_busco_df, second_busco_df
        )
        
        # Validate output data
        self._validate_output_data(ortholog_df, paralog_df)
        
        logger.info("Alignment analysis workflow completed successfully")
        return ortholog_df, paralog_df
    
    def _validate_input_data(self, first_busco_df: pd.DataFrame, second_busco_df: pd.DataFrame):
        """Validate input BUSCO dataframes."""
        required_columns = ['busco_id', 'gene_sequence', 'gene_length']
        
        for df_name, df in [('first', first_busco_df), ('second', second_busco_df)]:
            missing_cols = [col for col in required_columns if col not in df.columns]
            if missing_cols:
                raise ValueError(f"{df_name} BUSCO dataframe missing columns: {missing_cols}")
            
            if len(df) == 0:
                raise ValueError(f"{df_name} BUSCO dataframe is empty")
        
        logger.info(f"  Input validation passed")
        logger.info(f"    First dataset: {len(first_busco_df)} sequences")
        logger.info(f"    Second dataset: {len(second_busco_df)} sequences")
    
    def _validate_output_data(self, ortholog_df: pd.DataFrame, paralog_df: pd.DataFrame):
        """Validate output dataframes."""
        if len(ortholog_df) == 0:
            logger.warning("No ortholog pairs found - check input data and thresholds")
        else:
            logger.info(f"  Output validation passed: {len(ortholog_df)} ortholog pairs")
        
        # Check for required columns in ortholog dataframe
        expected_cols = ['busco_id', 'similarity', 'confidence']
        missing_cols = [col for col in expected_cols if col not in ortholog_df.columns]
        if missing_cols:
            logger.warning(f"Output dataframe missing expected columns: {missing_cols}")
    
    def get_alignment_statistics(self) -> Dict:
        """Get comprehensive alignment statistics."""
        # Get statistics from engine components
        stats = {}
        
        if hasattr(self.engine.partitioner, 'partitions'):
            stats['partitioning'] = self.engine.partitioner.get_partition_statistics()
        
        return stats