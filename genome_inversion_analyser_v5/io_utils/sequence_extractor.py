
# =============================================================================
# Sequence Extractor (sequence_extractor.py)
# =============================================================================

"""
Enhanced sequence extraction with comprehensive validation and error handling.
Handles gene boundary validation, strand processing, and translation checks.
"""

import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
from Bio.Seq import Seq

from ..logger import get_logger
from .fasta_loader import FastaLoader

logger = get_logger()

class SequenceExtractor:
    """
    Enhanced sequence extractor with comprehensive validation and error handling.
    Extracts BUSCO sequences with proper strand handling and validation.
    """
    
    def __init__(self, fasta_loader: FastaLoader):
        """
        Initialize sequence extractor.
        
        Args:
            fasta_loader: Initialized FastaLoader instance
        """
        self.fasta_loader = fasta_loader
        self.extraction_stats = {}
    
    def extract_busco_sequences(self, busco_df: pd.DataFrame, config) -> pd.DataFrame:
        """
        Extract BUSCO sequences with comprehensive validation.
        
        Args:
            busco_df: DataFrame with BUSCO records
            config: Configuration object
            
        Returns:
            DataFrame with extracted sequences
        """
        logger.info(f"Extracting BUSCO sequences with enhanced validation...")
        
        # Ensure genome sequences are loaded
        if not self.fasta_loader._indexed:
            self.fasta_loader.load_sequences()
        
        busco_with_seqs = []
        extraction_stats = {
            'total_attempted': len(busco_df),
            'successful': 0,
            'failed': 0,
            'warnings': 0
        }
        
        # Process each gene
        for idx, gene in busco_df.iterrows():
            try:
                extracted_info, validation_info = self._extract_and_validate_gene(gene, config)
                
                if extracted_info and validation_info['valid']:
                    # Create enhanced gene record
                    gene_record = self._create_gene_record(gene, extracted_info, validation_info, config)
                    busco_with_seqs.append(gene_record)
                    extraction_stats['successful'] += 1
                    
                    if validation_info.get('warnings'):
                        extraction_stats['warnings'] += 1
                else:
                    extraction_stats['failed'] += 1
                    if config.get('enable_debug_output', False) and extraction_stats['failed'] <= 10:
                        errors = '; '.join(validation_info.get('errors', ['Unknown error']))
                        logger.warning(f"  Failed to extract {gene['busco_id']}: {errors}")
                        
            except Exception as e:
                extraction_stats['failed'] += 1
                if config.get('enable_debug_output', False) and extraction_stats['failed'] <= 10:
                    logger.warning(f"  Extraction error for {gene['busco_id']}: {e}")
        
        # Report statistics
        self._report_extraction_statistics(extraction_stats, config)
        
        busco_seq_df = pd.DataFrame(busco_with_seqs)
        self.extraction_stats = extraction_stats
        
        return busco_seq_df
    
    def _extract_and_validate_gene(self, gene_info: pd.Series, config) -> Tuple[Optional[Dict], Dict]:
        """
        Extract and validate a single gene sequence.
        
        Args:
            gene_info: Series with gene information
            config: Configuration object
            
        Returns:
            Tuple of (extracted_info or None, validation_info)
        """
        seq_id = gene_info['sequence']
        start = int(gene_info['gene_start'])
        end = int(gene_info['gene_end'])
        strand = gene_info['strand']
        busco_id = gene_info['busco_id']
        
        validation_info = {
            'valid': False,
            'warnings': [],
            'errors': []
        }
        
        # Validate coordinates
        is_valid, error_msg = self.fasta_loader.validate_coordinates(seq_id, start, end)
        if not is_valid:
            validation_info['errors'].append(error_msg)
            return None, validation_info
        
        # Extract sequence
        gene_seq = self.fasta_loader.extract_subsequence(seq_id, start, end, strand)
        if gene_seq is None:
            validation_info['errors'].append("Failed to extract sequence")
            return None, validation_info
        
        # Additional validations
        if config.get('enable_gene_boundary_validation', False):
            boundary_warnings = self._validate_gene_boundaries(gene_seq, busco_id, config)
            validation_info['warnings'].extend(boundary_warnings)
        
        # Translation check if enabled
        if config.get('enable_translation_check', False):
            translation_info = self._validate_translation(gene_seq, busco_id)
            validation_info.update(translation_info)
        
        # Final validation
        if len(gene_seq) < 10:
            validation_info['warnings'].append(f"Very short gene sequence ({len(gene_seq)} bp)")
        
        validation_info['valid'] = len(validation_info['errors']) == 0
        
        extracted_info = {
            'sequence': gene_seq,
            'validated_start': start,
            'validated_end': end,
            'validated_strand': strand,
            'length': len(gene_seq),
            'chromosome_length': self.fasta_loader.get_sequence_length(seq_id)
        }
        
        return extracted_info, validation_info
    
    def _validate_gene_boundaries(self, gene_seq: str, busco_id: str, config) -> List[str]:
        """
        Validate gene boundaries and sequence quality.
        
        Args:
            gene_seq: Gene sequence
            busco_id: BUSCO identifier
            config: Configuration object
            
        Returns:
            List of warning messages
        """
        warnings = []
        
        # Check for unusual composition
        gc_content = (gene_seq.count('G') + gene_seq.count('C')) / len(gene_seq)
        if gc_content < 0.2 or gc_content > 0.8:
            warnings.append(f"Unusual GC content: {gc_content:.2f}")
        
        # Check for excessive N's
        n_content = gene_seq.count('N') / len(gene_seq)
        if n_content > 0.1:
            warnings.append(f"High N content: {n_content:.2f}")
        
        return warnings
    
    def _validate_translation(self, gene_seq: str, busco_id: str) -> Dict:
        """
        Validate gene translation for CDS quality.
        
        Args:
            gene_seq: Gene sequence
            busco_id: BUSCO identifier
            
        Returns:
            Dictionary with translation validation info
        """
        validation_info = {'translation_warnings': [], 'translation_errors': []}
        
        try:
            # Check if sequence length is multiple of 3
            if len(gene_seq) % 3 != 0:
                validation_info['translation_warnings'].append("Sequence length not multiple of 3")
            
            # Attempt translation
            protein_seq = str(Seq(gene_seq).translate())
            
            # Check for premature stop codons
            internal_stops = protein_seq[:-1].count('*')
            if internal_stops > 0:
                validation_info['translation_warnings'].append(f"{internal_stops} internal stop codons")
            
            # Check for start codon
            if len(gene_seq) >= 3:
                start_codon = gene_seq[:3].upper()
                if start_codon not in ['ATG', 'GTG', 'TTG']:
                    validation_info['translation_warnings'].append(f"Non-standard start codon: {start_codon}")
            
            # Check for reasonable protein length
            if len(protein_seq) < 50:
                validation_info['translation_warnings'].append(f"Short protein sequence ({len(protein_seq)} aa)")
                
        except Exception as e:
            validation_info['translation_errors'].append(f"Translation failed: {e}")
        
        return validation_info
    
    def _create_gene_record(self, gene: pd.Series, extracted_info: Dict, 
                          validation_info: Dict, config) -> Dict:
        """
        Create enhanced gene record with all information.
        
        Args:
            gene: Original gene information
            extracted_info: Extracted sequence information
            validation_info: Validation results
            config: Configuration object
            
        Returns:
            Complete gene record dictionary
        """
        gene_record = {
            'busco_id': gene['busco_id'],
            'sequence_id': gene['sequence'],
            'gene_start': extracted_info['validated_start'],
            'gene_end': extracted_info['validated_end'],
            'strand': extracted_info['validated_strand'],
            'gene_sequence': extracted_info['sequence'],
            'gene_length': extracted_info['length'],
            'status': gene['status'],
            'is_paralog': gene.get('is_paralog', False),
            'paralog_rank': gene.get('paralog_rank', 1),
            'paralog_count': gene.get('paralog_count', 1),
            'chromosome_length': extracted_info['chromosome_length']
        }
        
        # Add validation information if debug enabled
        if config.get('enable_debug_output', False):
            gene_record.update({
                'validation_warnings': ';'.join(validation_info.get('warnings', [])),
                'translation_warnings': ';'.join(validation_info.get('translation_warnings', [])),
                'extraction_method': 'enhanced'
            })
        
        return gene_record
    
    def _report_extraction_statistics(self, extraction_stats: Dict, config):
        """Report extraction statistics."""
        success_rate = extraction_stats['successful'] / extraction_stats['total_attempted'] * 100
        
        logger.info(f"  Extraction completed: {extraction_stats['successful']}/{extraction_stats['total_attempted']} ({success_rate:.1f}% success)")
        
        if extraction_stats['warnings'] > 0:
            logger.info(f"  {extraction_stats['warnings']} genes extracted with warnings")
        
        if extraction_stats['failed'] > 10 and config.get('enable_debug_output', False):
            logger.warning(f"  ... and {extraction_stats['failed'] - 10} more extraction failures")
    
    def get_extraction_statistics(self) -> Dict:
        """Get extraction statistics."""
        return self.extraction_stats.copy()