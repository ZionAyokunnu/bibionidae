
# =============================================================================
# Sequence Partitioner (partition.py)
# =============================================================================

"""
Intelligent sequence partitioning for optimal alignment method selection.
Partitions sequence pairs by length and complexity for hybrid alignment.
"""

from typing import List, Dict, Tuple, Any
from collections import defaultdict

from ..logger import get_logger

logger = get_logger()

class SequencePartitioner:
    """
    Partitions sequence pairs by length and complexity for optimal alignment method selection.
    Implements intelligent decision-making for hybrid alignment strategies.
    """
    
    def __init__(self, config):
        """
        Initialize sequence partitioner.
        
        Args:
            config: Configuration object with partitioning thresholds
        """
        self.config = config
        self.short_threshold = config.get('short_sequence_threshold', 500)
        self.long_threshold = config.get('long_sequence_threshold', 1500)
        self.partitions = {}
        
        logger.info(f"Sequence partitioner initialized:")
        logger.info(f"  Short threshold: ≤{self.short_threshold} bp → Biopython")
        logger.info(f"  Long threshold: ≥{self.long_threshold} bp → Minimap2")
        logger.info(f"  Buffer zone: {self.short_threshold}-{self.long_threshold} bp → {config.get('buffer_zone_method', 'dual')}")
    
    def partition_sequence_pairs(self, sequence_pairs: List[Dict]) -> Dict[str, List[Dict]]:
        """
        Partition sequence pairs based on length and complexity.
        
        Args:
            sequence_pairs: List of sequence pair dictionaries
            
        Returns:
            Dictionary with partitioned sequence pairs
        """
        logger.info(f"Partitioning {len(sequence_pairs)} sequence pairs...")
        
        partitions = {
            'short_pairs': [],      # Use Biopython (fast and accurate for short sequences)
            'long_pairs': [],       # Use Minimap2 (designed for long sequences)
            'buffer_pairs': [],     # Use configurable method or both
            'mixed_pairs': []       # One short, one long - special handling
        }
        
        for pair in sequence_pairs:
            partition_key = self._classify_sequence_pair(pair)
            partitions[partition_key].append(pair)
        
        # Log partition statistics
        self._log_partition_statistics(partitions)
        
        self.partitions = partitions
        return partitions
    
    def _classify_sequence_pair(self, pair: Dict) -> str:
        """
        Classify a sequence pair into the appropriate partition.
        
        Args:
            pair: Sequence pair dictionary
            
        Returns:
            Partition key string
        """
        # Extract sequence lengths
        len1 = pair['first_gene']['gene_length']
        len2 = pair['second_gene']['gene_length']
        
        max_len = max(len1, len2)
        min_len = min(len1, len2)
        
        # Classify based on length thresholds
        if max_len <= self.short_threshold:
            return 'short_pairs'
        elif min_len >= self.long_threshold:
            return 'long_pairs'
        elif self.short_threshold < max_len < self.long_threshold:
            return 'buffer_pairs'
        else:
            # One sequence much longer than the other
            return 'mixed_pairs'
    
    def _log_partition_statistics(self, partitions: Dict[str, List]):
        """Log detailed partition statistics."""
        total_pairs = sum(len(pairs) for pairs in partitions.values())
        
        if self.config.get('detailed_alignment_logging', False):
            logger.info("  Detailed partitioning results:")
            for partition_name, pairs in partitions.items():
                percentage = len(pairs) / total_pairs * 100 if total_pairs > 0 else 0
                logger.info(f"    {partition_name}: {len(pairs)} pairs ({percentage:.1f}%)")
                
                # Log length statistics for each partition
                if pairs and self.config.get('enable_debug_output', False):
                    lengths = []
                    for pair in pairs:
                        lengths.extend([pair['first_gene']['gene_length'], 
                                      pair['second_gene']['gene_length']])
                    
                    if lengths:
                        logger.info(f"      Length range: {min(lengths)}-{max(lengths)} bp")
                        logger.info(f"      Average length: {sum(lengths)/len(lengths):.0f} bp")
        else:
            logger.info(f"  Partitioned into:")
            logger.info(f"    Short pairs (≤{self.short_threshold}bp): {len(partitions['short_pairs'])}")
            logger.info(f"    Long pairs (≥{self.long_threshold}bp): {len(partitions['long_pairs'])}")
            logger.info(f"    Buffer zone pairs: {len(partitions['buffer_pairs'])}")
            logger.info(f"    Mixed length pairs: {len(partitions['mixed_pairs'])}")
    
    def get_partition_statistics(self) -> Dict[str, Any]:
        """
        Get comprehensive partition statistics.
        
        Returns:
            Dictionary with partition statistics
        """
        if not self.partitions:
            return {}
        
        total_pairs = sum(len(pairs) for pairs in self.partitions.values())
        
        stats = {
            'total_pairs': total_pairs,
            'partition_counts': {name: len(pairs) for name, pairs in self.partitions.items()},
            'partition_percentages': {
                name: len(pairs) / total_pairs * 100 if total_pairs > 0 else 0 
                for name, pairs in self.partitions.items()
            },
            'short_threshold': self.short_threshold,
            'long_threshold': self.long_threshold
        }
        
        return stats