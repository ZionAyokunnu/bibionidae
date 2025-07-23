# =============================================================================
# FASTA Loader (fasta_loader.py)
# =============================================================================

"""
FASTA file loading and genome sequence management.
Handles large genome files efficiently with optional indexing.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Iterator, Tuple
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gzip
import numpy as np

from ..logger import get_logger
from ..utils import standardize_sequence_id

logger = get_logger()

class FastaLoader:
    """
    Efficient FASTA file loader with sequence standardization and validation.
    Supports both regular and gzipped FASTA files.
    """
    
    def __init__(self, fasta_path: str, index_sequences: bool = True):
        """
        Initialize FASTA loader.
        
        Args:
            fasta_path: Path to FASTA file
            index_sequences: Whether to create an index for fast random access
        """
        self.fasta_path = Path(fasta_path)
        self.index_sequences = index_sequences
        self.sequences = {}
        self.sequence_lengths = {}
        self.sequence_stats = {}
        self._indexed = False
        
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_path}")
        
        logger.info(f"Initializing FASTA loader for: {self.fasta_path}")
    
    def load_sequences(self, standardize_ids: bool = True) -> Dict[str, str]:
        """
        Load all sequences into memory.
        
        Args:
            standardize_ids: Whether to standardize sequence IDs
            
        Returns:
            Dictionary mapping sequence IDs to sequences
        """
        logger.info(f"Loading sequences from {self.fasta_path}")
        
        sequences = {}
        sequence_lengths = {}
        
        try:
            # Handle gzipped files
            if self.fasta_path.suffix.lower() == '.gz':
                file_handle = gzip.open(self.fasta_path, 'rt')
            else:
                file_handle = open(self.fasta_path, 'r')
            
            with file_handle:
                for record in SeqIO.parse(file_handle, "fasta"):
                    # Standardize sequence ID if requested
                    seq_id = record.id
                    if standardize_ids:
                        seq_id = standardize_sequence_id(record.id)
                        if seq_id is None:
                            logger.warning(f"Skipping sequence with invalid ID: {record.id}")
                            continue
                    
                    # Store sequence and length
                    sequences[seq_id] = str(record.seq)
                    sequence_lengths[seq_id] = len(record.seq)
                    
        except Exception as e:
            logger.error(f"Failed to load FASTA file {self.fasta_path}: {e}")
            raise
        
        # Calculate statistics
        self._calculate_genome_stats(sequence_lengths)
        
        self.sequences = sequences
        self.sequence_lengths = sequence_lengths
        self._indexed = True
        
        logger.info(f"Loaded {len(sequences)} sequences")
        logger.info(f"Total genome size: {sum(sequence_lengths.values()):,} bp")
        
        return sequences
    
    def get_sequence(self, seq_id: str) -> Optional[str]:
        """
        Get a specific sequence by ID.
        
        Args:
            seq_id: Sequence identifier
            
        Returns:
            Sequence string or None if not found
        """
        if not self._indexed:
            self.load_sequences()
        
        return self.sequences.get(seq_id)
    
    def get_sequence_length(self, seq_id: str) -> Optional[int]:
        """
        Get length of a specific sequence.
        
        Args:
            seq_id: Sequence identifier
            
        Returns:
            Sequence length or None if not found
        """
        if not self._indexed:
            self.load_sequences()
        
        return self.sequence_lengths.get(seq_id)
    
    def get_sequence_names(self) -> List[str]:
        """
        Get list of all sequence names.
        
        Returns:
            List of sequence identifiers
        """
        if not self._indexed:
            self.load_sequences()
        
        return list(self.sequences.keys())
    
    def get_genome_stats(self) -> Dict[str, any]:
        """
        Get comprehensive genome statistics.
        
        Returns:
            Dictionary with genome statistics
        """
        if not self._indexed:
            self.load_sequences()
        
        return self.sequence_stats.copy()
    
    def extract_subsequence(self, seq_id: str, start: int, end: int, 
                          strand: str = '+') -> Optional[str]:
        """
        Extract subsequence from a chromosome.
        
        Args:
            seq_id: Sequence identifier
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            strand: Strand orientation ('+' or '-')
            
        Returns:
            Extracted subsequence or None if invalid
        """
        if not self._indexed:
            self.load_sequences()
        
        if seq_id not in self.sequences:
            return None
        
        full_seq = self.sequences[seq_id]
        
        # Convert to 0-based indexing
        start_idx = max(0, start - 1)
        end_idx = min(len(full_seq), end)
        
        if start_idx >= end_idx:
            return None
        
        # Extract sequence
        subseq = full_seq[start_idx:end_idx]
        
        # Handle reverse complement for negative strand
        if strand == '-':
            from Bio.Seq import Seq
            subseq = str(Seq(subseq).reverse_complement())
        
        return subseq
    
    def validate_coordinates(self, seq_id: str, start: int, end: int) -> Tuple[bool, str]:
        """
        Validate genomic coordinates.
        
        Args:
            seq_id: Sequence identifier
            start: Start position
            end: End position
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not self._indexed:
            self.load_sequences()
        
        # Check if sequence exists
        if seq_id not in self.sequences:
            return False, f"Sequence '{seq_id}' not found"
        
        seq_length = self.sequence_lengths[seq_id]
        
        # Check coordinate validity
        if start < 1:
            return False, f"Start position {start} is less than 1"
        
        if end > seq_length:
            return False, f"End position {end} exceeds sequence length {seq_length}"
        
        if start > end:
            return False, f"Start position {start} is greater than end position {end}"
        
        return True, ""
    
    def _calculate_genome_stats(self, sequence_lengths: Dict[str, int]) -> None:
        """Calculate comprehensive genome statistics."""
        if not sequence_lengths:
            self.sequence_stats = {}
            return
        
        lengths = list(sequence_lengths.values())
        total_length = sum(lengths)
        n_contigs = len(lengths)
        
        # Sort lengths for N50/N90 calculation
        sorted_lengths = sorted(lengths, reverse=True)
        cumsum = 0
        n50 = n90 = 0
        
        for length in sorted_lengths:
            cumsum += length
            if n50 == 0 and cumsum >= total_length * 0.5:
                n50 = length
            if n90 == 0 and cumsum >= total_length * 0.9:
                n90 = length
                break
        
        self.sequence_stats = {
            'total_length': total_length,
            'n_contigs': n_contigs,
            'n50': n50,
            'n90': n90,
            'max_contig': max(lengths),
            'min_contig': min(lengths),
            'mean_contig_length': total_length / n_contigs,
            'median_contig_length': sorted_lengths[n_contigs // 2] if n_contigs > 0 else 0,
            'gc_content': None  # Could be calculated if needed
        }
    
    def __len__(self) -> int:
        """Return number of sequences."""
        if not self._indexed:
            self.load_sequences()
        return len(self.sequences)
    
    def __contains__(self, seq_id: str) -> bool:
        """Check if sequence ID exists."""
        if not self._indexed:
            self.load_sequences()
        return seq_id in self.sequences