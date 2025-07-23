
# =============================================================================
# Integration Module (workflow.py)
# =============================================================================

"""
Workflow integration for I/O operations.
Orchestrates the complete I/O pipeline from input to output.
"""

from pathlib import Path
from typing import Dict, Tuple, Any

from .fasta_loader import FastaLoader
from .busco_parser import BuscoParser
from .sequence_extractor import SequenceExtractor
from .output_writer import OutputWriter
from ..logger import get_logger

logger = get_logger()

class IOWorkflow:
    """
    Orchestrates the complete I/O workflow for genome analysis.
    Handles loading, parsing, extraction, and output operations.
    """
    
    def __init__(self, config):
        """
        Initialize I/O workflow.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.output_writer = OutputWriter(Path(config.get('base_output_dir', 'results')))
    
    def load_and_process_genomes(self) -> Tuple[Dict, Dict]:
        """
        Load and process both genome datasets.
        
        Returns:
            Tuple of (first_genome_data, second_genome_data)
        """
        logger.info("Loading and processing genome datasets...")
        
        # Process first genome
        logger.info("Processing first genome...")
        first_data = self._process_single_genome(
            self.config.get('first_fasta_path'),
            self.config.get('first_busco_path'),
            'first'
        )
        
        # Process second genome
        logger.info("Processing second genome...")
        second_data = self._process_single_genome(
            self.config.get('second_fasta_path'),
            self.config.get('second_busco_path'),
            'second'
        )
        
        logger.info("Genome processing completed successfully")
        return first_data, second_data
    
    def _process_single_genome(self, fasta_path: str, busco_path: str, 
                             genome_name: str) -> Dict[str, Any]:
        """
        Process a single genome dataset.
        
        Args:
            fasta_path: Path to genome FASTA file
            busco_path: Path to BUSCO table
            genome_name: Name for this genome dataset
            
        Returns:
            Dictionary with processed genome data
        """
        # Load FASTA sequences
        fasta_loader = FastaLoader(fasta_path)
        sequences = fasta_loader.load_sequences()
        genome_stats = fasta_loader.get_genome_stats()
        
        # Parse BUSCO table
        busco_parser = BuscoParser(busco_path)
        busco_raw_df = busco_parser.parse_busco_table(self.config)
        
        # Filter BUSCO genes
        busco_filtered_df = busco_parser.filter_busco_genes(busco_raw_df, self.config)
        
        # Extract sequences
        sequence_extractor = SequenceExtractor(fasta_loader)
        busco_sequences_df = sequence_extractor.extract_busco_sequences(busco_filtered_df, self.config)
        
        # Save intermediate results
        self.output_writer.save_busco_results(
            busco_raw_df, 
            f'{genome_name}_busco_raw.csv',
            busco_parser.get_parsing_statistics()
        )
        
        self.output_writer.save_busco_results(
            busco_filtered_df,
            f'{genome_name}_busco_filtered.csv'
        )
        
        self.output_writer.save_sequence_data(
            busco_sequences_df,
            f'{genome_name}_sequences.csv',
            sequence_extractor.get_extraction_statistics()
        )
        
        return {
            'fasta_loader': fasta_loader,
            'sequences': sequences,
            'genome_stats': genome_stats,
            'busco_raw_df': busco_raw_df,
            'busco_filtered_df': busco_filtered_df,
            'busco_sequences_df': busco_sequences_df,
            'parsing_stats': busco_parser.get_parsing_statistics(),
            'extraction_stats': sequence_extractor.get_extraction_statistics()
        }
    
    def save_final_results(self, results_dict: Dict, quality_data: Dict = None) -> Dict[str, Path]:
        """
        Save all final analysis results.
        
        Args:
            results_dict: Dictionary with analysis results
            quality_data: Optional quality assessment data
            
        Returns:
            Dictionary mapping result types to file paths
        """
        logger.info("Saving final analysis results...")
        
        # Save main analysis results
        saved_files = self.output_writer.save_analysis_results(results_dict, self.config)
        
        # Save quality report if available
        if quality_data:
            quality_list = []
            for genome_name, quality_info in quality_data.items():
                quality_record = {'genome': genome_name}
                quality_record.update(quality_info.get('metrics', {}))
                quality_record['quality_score'] = quality_info.get('quality_score', 0.0)
                quality_record['quality_class'] = quality_info.get('quality_class', 'unknown')
                quality_list.append(quality_record)
            
            quality_path = self.output_writer.save_quality_report(quality_list, self.config)
            saved_files['quality_report'] = quality_path
        
        # Generate summary report
        summary_path = self.output_writer.generate_summary_report(results_dict, self.config)
        saved_files['summary_report'] = summary_path
        
        logger.info(f"All results saved to: {self.output_writer.output_dir}")
        return saved_files