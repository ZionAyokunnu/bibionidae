
# =============================================================================
# Output Writer (output_writer.py)
# =============================================================================

"""
Comprehensive output writer for analysis results.
Handles CSV, JSON, and formatted report generation.
"""

import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime

from ..logger import get_logger

logger = get_logger()

class OutputWriter:
    """
    Comprehensive output writer for all analysis results.
    Supports multiple formats and includes metadata preservation.
    """
    
    def __init__(self, output_dir: Path):
        """
        Initialize output writer.
        
        Args:
            output_dir: Base output directory
        """
        self.output_dir = Path(output_dir)
        self.data_dir = self.output_dir / 'data'
        self.reports_dir = self.output_dir / 'reports'
        
        # Ensure directories exist
        self.data_dir.mkdir(exist_ok=True)
        self.reports_dir.mkdir(exist_ok=True)
        
        logger.info(f"Output writer initialized for: {self.output_dir}")
    
    def save_busco_results(self, busco_df: pd.DataFrame, filename: str, 
                          parsing_stats: Optional[Dict] = None) -> Path:
        """
        Save BUSCO parsing results with metadata.
        
        Args:
            busco_df: DataFrame with BUSCO data
            filename: Output filename
            parsing_stats: Optional parsing statistics
            
        Returns:
            Path to saved file
        """
        output_path = self.data_dir / filename
        
        # Save main data
        busco_df.to_csv(output_path, index=False)
        
        # Save metadata if available
        if parsing_stats:
            metadata_path = output_path.with_suffix('.metadata.json')
            metadata = {
                'created_at': datetime.now().isoformat(),
                'record_count': len(busco_df),
                'parsing_statistics': parsing_stats,
                'columns': list(busco_df.columns)
            }
            
            with open(metadata_path, 'w') as f:
                json.dump(metadata, f, indent=2)
        
        logger.info(f"  BUSCO results saved: {output_path}")
        return output_path
    
    def save_sequence_data(self, sequence_df: pd.DataFrame, filename: str,
                          extraction_stats: Optional[Dict] = None) -> Path:
        """
        Save extracted sequence data.
        
        Args:
            sequence_df: DataFrame with sequence data
            filename: Output filename
            extraction_stats: Optional extraction statistics
            
        Returns:
            Path to saved file
        """
        output_path = self.data_dir / filename
        
        # Save main data (excluding sequences for size)
        save_df = sequence_df.drop(columns=['gene_sequence'], errors='ignore')
        save_df.to_csv(output_path, index=False)
        
        # Save sequences separately in FASTA format
        fasta_path = output_path.with_suffix('.fasta')
        self._save_sequences_as_fasta(sequence_df, fasta_path)
        
        # Save metadata
        if extraction_stats:
            metadata_path = output_path.with_suffix('.metadata.json')
            metadata = {
                'created_at': datetime.now().isoformat(),
                'record_count': len(sequence_df),
                'extraction_statistics': extraction_stats,
                'fasta_file': fasta_path.name
            }
            
            with open(metadata_path, 'w') as f:
                json.dump(metadata, f, indent=2)
        
        logger.info(f"  Sequence data saved: {output_path}")
        logger.info(f"  Sequences saved: {fasta_path}")
        return output_path
    
    def save_analysis_results(self, results_dict: Dict[str, pd.DataFrame], 
                            config) -> Dict[str, Path]:
        """
        Save all analysis results.
        
        Args:
            results_dict: Dictionary of analysis results
            config: Configuration object
            
        Returns:
            Dictionary mapping result types to saved file paths
        """
        logger.info("Saving analysis results...")
        
        saved_files = {}
        
        # Define output mappings
        output_mappings = {
            'ortholog_df': Path(config.get('synteny_analysis_csv', 'synteny_analysis.csv')).name,
            'inversion_df': Path(config.get('inversion_summary_csv', 'inversion_summary.csv')).name,
            'rearrangement_df': Path(config.get('chromosome_rearrangements_csv', 'chromosome_rearrangements.csv')).name,
            'synteny_df': 'synteny_blocks.csv',
            'mapping_df': 'chromosome_mappings.csv',
            'paralog_df': Path(config.get('paralog_analysis_csv', 'paralog_analysis.csv')).name
        }
        
        # Save each result type
        for result_key, filename in output_mappings.items():
            if result_key in results_dict and not results_dict[result_key].empty:
                output_path = self.data_dir / filename
                results_dict[result_key].to_csv(output_path, index=False)
                saved_files[result_key] = output_path
                logger.info(f"  {result_key} saved: {output_path}")
        
        return saved_files
    
    def save_quality_report(self, quality_data: List[Dict], config) -> Path:
        """
        Save assembly quality report.
        
        Args:
            quality_data: List of quality assessment dictionaries
            config: Configuration object
            
        Returns:
            Path to saved quality report
        """
        output_path = self.data_dir / Path(config.get('quality_report_csv', 'assembly_quality_report.csv')).name
        
        quality_df = pd.DataFrame(quality_data)
        quality_df.to_csv(output_path, index=False)
        
        logger.info(f"  Quality report saved: {output_path}")
        return output_path
    
    def generate_summary_report(self, results_dict: Dict, config) -> Path:
        """
        Generate comprehensive text summary report.
        
        Args:
            results_dict: Dictionary with all analysis results
            config: Configuration object
            
        Returns:
            Path to generated report
        """
        report_path = self.reports_dir / 'analysis_summary.txt'
        
        with open(report_path, 'w') as f:
            f.write("ENHANCED GENOME SYNTENY AND INVERSION ANALYSIS REPORT\n")
            f.write("=" * 60 + "\n\n")
            
            # Analysis configuration
            f.write("ANALYSIS CONFIGURATION\n")
            f.write("-" * 25 + "\n")
            f.write(f"Alignment Strategy: {config.get('alignment_strategy', 'unknown')}\n")
            f.write(f"Enable Parallel Processing: {config.get('enable_parallel_alignment', False)}\n")
            f.write(f"Use Reciprocal Best Hits: {config.get('use_reciprocal_best_hits', False)}\n")
            f.write(f"Enable Paralog Detection: {config.get('enable_paralog_detection', False)}\n")
            f.write("\n")
            
            # Input data summary
            f.write("INPUT DATA SUMMARY\n")
            f.write("-" * 20 + "\n")
            f.write(f"First genome: {config.get('first_fasta_path', 'N/A')}\n")
            f.write(f"Second genome: {config.get('second_fasta_path', 'N/A')}\n")
            f.write(f"First BUSCO: {config.get('first_busco_path', 'N/A')}\n")
            f.write(f"Second BUSCO: {config.get('second_busco_path', 'N/A')}\n")
            f.write("\n")
            
            # Results summary
            self._write_results_summary(f, results_dict)
            
            # Quality assessment
            if 'first_quality' in results_dict and 'second_quality' in results_dict:
                self._write_quality_summary(f, results_dict)
            
            f.write(f"\nReport generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        logger.info(f"  Summary report saved: {report_path}")
        return report_path
    
    def _save_sequences_as_fasta(self, sequence_df: pd.DataFrame, fasta_path: Path):
        """Save sequences in FASTA format."""
        if 'gene_sequence' not in sequence_df.columns:
            return
        
        with open(fasta_path, 'w') as f:
            for _, row in sequence_df.iterrows():
                header = f">{row['busco_id']}_{row['sequence_id']}_{row['gene_start']}_{row['gene_end']}_{row['strand']}"
                f.write(f"{header}\n")
                f.write(f"{row['gene_sequence']}\n")
    
    def _write_results_summary(self, f, results_dict: Dict):
        """Write results summary to report file."""
        f.write("ANALYSIS RESULTS SUMMARY\n")
        f.write("-" * 25 + "\n")
        
        # Ortholog analysis
        if 'ortholog_df' in results_dict:
            ortholog_df = results_dict['ortholog_df']
            f.write(f"Ortholog pairs found: {len(ortholog_df)}\n")
            
            if len(ortholog_df) > 0:
                if 'similarity' in ortholog_df.columns:
                    f.write(f"Average similarity: {ortholog_df['similarity'].mean():.3f}\n")
                if 'confidence' in ortholog_df.columns:
                    f.write(f"Average confidence: {ortholog_df['confidence'].mean():.3f}\n")
                if 'alignment_method' in ortholog_df.columns:
                    method_counts = ortholog_df['alignment_method'].value_counts()
                    f.write("Alignment methods used:\n")
                    for method, count in method_counts.items():
                        f.write(f"  {method}: {count} ({count/len(ortholog_df)*100:.1f}%)\n")
        
        # Synteny analysis
        if 'synteny_df' in results_dict:
            synteny_df = results_dict['synteny_df']
            f.write(f"Synteny blocks found: {len(synteny_df)}\n")
            
            if len(synteny_df) > 0 and 'block_size' in synteny_df.columns:
                f.write(f"Average block size: {synteny_df['block_size'].mean():.1f} genes\n")
        
        # Inversion analysis
        if 'inversion_df' in results_dict:
            inversion_df = results_dict['inversion_df']
            f.write(f"Inversions detected: {len(inversion_df)}\n")
            
            if len(inversion_df) > 0 and 'size_genes' in inversion_df.columns:
                f.write(f"Average inversion size: {inversion_df['size_genes'].mean():.1f} genes\n")
        
        # Rearrangement analysis
        if 'rearrangement_df' in results_dict:
            rearrangement_df = results_dict['rearrangement_df']
            f.write(f"Chromosome rearrangements: {len(rearrangement_df)}\n")
        
        f.write("\n")
    
    def _write_quality_summary(self, f, results_dict: Dict):
        """Write quality assessment summary."""
        f.write("ASSEMBLY QUALITY ASSESSMENT\n")
        f.write("-" * 30 + "\n")
        
        first_quality = results_dict['first_quality']
        second_quality = results_dict['second_quality']
        
        f.write(f"First genome quality: {first_quality['quality_class']} (score: {first_quality['quality_score']:.3f})\n")
        f.write(f"Second genome quality: {second_quality['quality_class']} (score: {second_quality['quality_score']:.3f})\n")
        
        # Detailed metrics if available
        for genome_name, quality_info in [('First', first_quality), ('Second', second_quality)]:
            if 'metrics' in quality_info:
                metrics = quality_info['metrics']
                f.write(f"\n{genome_name} genome metrics:\n")
                
                if 'n50' in metrics:
                    f.write(f"  N50: {metrics['n50']:,} bp\n")
                if 'n_contigs' in metrics:
                    f.write(f"  Contigs: {metrics['n_contigs']:,}\n")
                if 'busco_completeness' in metrics:
                    f.write(f"  BUSCO completeness: {metrics['busco_completeness']:.1%}\n")
        
        f.write("\n")