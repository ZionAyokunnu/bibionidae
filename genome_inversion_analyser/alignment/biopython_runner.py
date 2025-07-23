
# =============================================================================
# Biopython Alignment Runner (biopython_runner.py)
# =============================================================================

"""
Biopython-based sequence alignment with parallel processing support.
Optimized for short to medium-length sequences with high accuracy requirements.
"""

import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from typing import List, Dict, Optional
from Bio.Align import PairwiseAligner
from difflib import SequenceMatcher

from ..logger import get_logger
from .alignment_result import AlignmentResult

logger = get_logger()

class BioPythonAligner:
    """
    Biopython-based sequence aligner with parallel processing and fallback support.
    Optimized for accuracy on short to medium sequences.
    """
    
    def __init__(self, config):
        """
        Initialize Biopython aligner.
        
        Args:
            config: Configuration object with alignment parameters
        """
        self.config = config
        self.aligner = self._setup_aligner()
        
        logger.info("Biopython aligner initialized")
        logger.info(f"  Mode: {config.get('biopython_mode', 'local')}")
        logger.info(f"  Parallel processing: {config.get('enable_parallel_alignment', True)}")
        logger.info(f"  Batch size: {config.get('biopython_batch_size', 100)}")
    
    def _setup_aligner(self) -> Optional[PairwiseAligner]:
        """Setup Biopython pairwise aligner with configuration."""
        try:
            aligner = PairwiseAligner()
            aligner.match_score = self.config.get('biopython_match_score', 2)
            aligner.mismatch_score = self.config.get('biopython_mismatch_score', -1)
            aligner.open_gap_score = self.config.get('biopython_gap_open_score', -2)
            aligner.extend_gap_score = self.config.get('biopython_gap_extend_score', -0.5)
            aligner.mode = self.config.get('biopython_mode', 'local')
            
            logger.info("  Biopython aligner configured successfully")
            return aligner
        except Exception as e:
            logger.warning(f"  Failed to setup Biopython aligner: {e}")
            return None
    
    def align_sequences(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """
        Align sequence pairs using Biopython with optional parallel processing.
        
        Args:
            sequence_pairs: List of sequence pair dictionaries
            
        Returns:
            List of AlignmentResult objects
        """
        if not sequence_pairs:
            return []
        
        logger.info(f"Running Biopython alignment on {len(sequence_pairs)} sequence pairs...")
        
        # Decide on processing strategy
        if (self.config.get('enable_parallel_alignment', True) and 
            len(sequence_pairs) > 200 and 
            multiprocessing.cpu_count() > 1):
            return self._parallel_alignment(sequence_pairs)
        else:
            return self._sequential_alignment(sequence_pairs)
    
    def _parallel_alignment(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """Run parallel Biopython alignment."""
        batch_size = self.config.get('biopython_batch_size', 100)
        max_workers = min(multiprocessing.cpu_count(), self.config.get('minimap2_threads', 4))
        
        # Split into batches
        batches = [sequence_pairs[i:i + batch_size] for i in range(0, len(sequence_pairs), batch_size)]
        
        logger.info(f"  Using parallel processing with {max_workers} workers, {len(batches)} batches")
        
        all_results = []
        
        # Extract serializable config parameters
        simple_config = self._extract_serializable_config()
        
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = []
                for batch in batches:
                    future = executor.submit(align_sequence_batch, batch, simple_config)
                    futures.append(future)
                
                for i, future in enumerate(futures):
                    try:
                        timeout = simple_config.get('timeout_per_alignment', 30)
                        batch_results = future.result(timeout=timeout)
                        all_results.extend(batch_results)
                        
                        if i % 5 == 0:  # Progress every 5 batches
                            logger.info(f"    Progress: {i+1}/{len(batches)} batches completed")
                            
                    except Exception as e:
                        logger.warning(f"  Batch {i+1} failed: {e}")
                        # Try sequential processing for failed batch
                        try:
                            batch_results = align_sequence_batch(batches[i], simple_config)
                            all_results.extend(batch_results)
                            logger.info(f"    Batch {i+1} recovered with sequential processing")
                        except Exception as e2:
                            logger.error(f"    Batch {i+1} failed completely: {e2}")
                            
        except Exception as e:
            logger.warning(f"  Parallel processing failed: {e}. Falling back to sequential...")
            return self._sequential_alignment(sequence_pairs)
        
        return all_results
    
    def _sequential_alignment(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """Run sequential Biopython alignment."""
        logger.info(f"  Using sequential processing for {len(sequence_pairs)} sequence pairs")
        
        batch_size = self.config.get('biopython_batch_size', 100)
        all_results = []
        
        for i in range(0, len(sequence_pairs), batch_size):
            batch = sequence_pairs[i:i + batch_size]
            batch_results = align_sequence_batch(batch, self.config)
            all_results.extend(batch_results)
            
            # Progress reporting
            if i % (batch_size * 10) == 0 or i + batch_size >= len(sequence_pairs):
                processed = min(i + batch_size, len(sequence_pairs))
                logger.info(f"    Processed {processed}/{len(sequence_pairs)} pairs ({processed/len(sequence_pairs)*100:.1f}%)")
        
        return all_results
    
    def _extract_serializable_config(self) -> Dict:
        """Extract only serializable configuration parameters."""
        return {
            'biopython_match_score': self.config.get('biopython_match_score', 2),
            'biopython_mismatch_score': self.config.get('biopython_mismatch_score', -1),
            'biopython_gap_open_score': self.config.get('biopython_gap_open_score', -2),
            'biopython_gap_extend_score': self.config.get('biopython_gap_extend_score', -0.5),
            'biopython_mode': self.config.get('biopython_mode', 'local'),
            'fallback_to_simple_similarity': self.config.get('fallback_to_simple_similarity', True),
            'timeout_per_alignment': self.config.get('timeout_per_alignment', 30),
            'detailed_alignment_logging': self.config.get('detailed_alignment_logging', False)
        }

# Standalone function for multiprocessing (must be at module level)
def align_sequence_batch(sequence_pairs_batch: List[Dict], config: Dict) -> List[AlignmentResult]:
    """
    Align a batch of sequence pairs using Biopython.
    This function must be at module level for multiprocessing.
    
    Args:
        sequence_pairs_batch: Batch of sequence pairs
        config: Configuration dictionary
        
    Returns:
        List of AlignmentResult objects
    """
    results = []
    
    # Setup aligner
    aligner = None
    use_aligner = False
    
    try:
        aligner = PairwiseAligner()
        aligner.match_score = config.get('biopython_match_score', 2)
        aligner.mismatch_score = config.get('biopython_mismatch_score', -1)
        aligner.open_gap_score = config.get('biopython_gap_open_score', -2)
        aligner.extend_gap_score = config.get('biopython_gap_extend_score', -0.5)
        aligner.mode = config.get('biopython_mode', 'local')
        use_aligner = True
    except Exception:
        use_aligner = False
    
    # Process each pair in the batch
    for pair_idx, pair in enumerate(sequence_pairs_batch):
        try:
            seq1 = pair['first_gene']['gene_sequence']
            seq2 = pair['second_gene']['gene_sequence']
            busco_id = pair['first_gene']['busco_id']
            
            result = None
            
            if use_aligner:
                # Try Biopython alignment
                try:
                    alignments = aligner.align(seq1, seq2)
                    if alignments:
                        best_alignment = alignments[0]
                        
                        # Calculate detailed statistics
                        alignment_stats = _calculate_alignment_statistics(best_alignment, seq1, seq2)
                        
                        result = AlignmentResult.from_biopython(
                            pair_index=pair_idx,
                            busco_id=busco_id,
                            bio_result=alignment_stats
                        )
                        
                except Exception:
                    # Biopython failed, will fall back to difflib
                    pass
            
            # Fallback to difflib if Biopython failed or not available
            if result is None and config.get('fallback_to_simple_similarity', True):
                similarity = SequenceMatcher(None, seq1, seq2).ratio()
                result = AlignmentResult.from_fallback(pair_idx, busco_id, similarity)
            
            if result:
                results.append(result)
                
        except Exception as e:
            # Log error but continue with next pair
            if config.get('detailed_alignment_logging', False):
                logger.warning(f"Failed to process sequence pair {pair_idx}: {e}")
            continue
    
    return results

def _calculate_alignment_statistics(alignment, seq1: str, seq2: str) -> Dict:
    """Calculate detailed alignment statistics from Biopython alignment."""
    try:
        alignment_str = str(alignment)
        alignment_lines = alignment_str.split('\n')
        
        if len(alignment_lines) >= 3:
            seq1_aligned = alignment_lines[0]
            seq2_aligned = alignment_lines[2]
            
            matches = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) 
                         if a == b and a != '-')
            alignment_length = len(seq1_aligned)
            
            # Coverage calculations
            query_coverage = (alignment_length - seq1_aligned.count('-')) / len(seq1)
            target_coverage = (alignment_length - seq2_aligned.count('-')) / len(seq2)
            identity = matches / alignment_length if alignment_length > 0 else 0
        else:
            # Fallback calculation
            matches = int(alignment.score / 2)  # Rough estimate
            alignment_length = max(len(seq1), len(seq2))
            query_coverage = 0.8  # Conservative estimate
            target_coverage = 0.8
            identity = matches / alignment_length if alignment_length > 0 else 0
    
    except Exception:
        # Simple fallback
        matches = int(alignment.score / 2)
        alignment_length = max(len(seq1), len(seq2))
        query_coverage = 0.7
        target_coverage = 0.7
        identity = matches / alignment_length if alignment_length > 0 else 0
    
    return {
        'identity': identity,
        'query_coverage': query_coverage,
        'target_coverage': target_coverage,
        'alignment_length': alignment_length,
        'matches': matches,
        'score': alignment.score
    }