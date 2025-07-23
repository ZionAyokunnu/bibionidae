
# =============================================================================
# I/O Utils Package (__init__.py)
# =============================================================================

"""
I/O utilities for genome analysis pipeline.
Handles FASTA loading, BUSCO parsing, and results output.
"""

from .fasta_loader import FastaLoader
from .busco_parser import BuscoParser, BuscoRecord
from .output_writer import OutputWriter
from .sequence_extractor import SequenceExtractor

__all__ = [
    'FastaLoader', 
    'BuscoParser', 
    'BuscoRecord',
    'OutputWriter',
    'SequenceExtractor'
]
