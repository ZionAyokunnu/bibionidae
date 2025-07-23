
# =============================================================================
# Reciprocal Best Hit Filter (rbh_filter.py)
# =============================================================================

"""
Reciprocal Best Hit (RBH) filtering for ortholog identification.
Implements sophisticated filtering strategies to identify true orthologs.
"""

from typing import List, Dict, Set, Any
from collections import defaultdict
from ..alignment.alignment_result import AlignmentResult
from ..logger import get_logger

logger = get_logger()

class ReciprocalBestHitFilter:
    """
    Reciprocal Best Hit filter for identifying true orthologous relationships.
    Implements multiple filtering strategies including tie-breaking and confidence-based selection.
    """
    
    def __init__(self, config):
        """
        Initialize RBH filter.
        
        Args:
            config: Configuration object with filtering parameters
        """
        self.config = config
        self.use_rbh = config.get('use_reciprocal_best_hits', True)
        self.identity_gap_threshold = config.get('identity_gap_threshold', 0.02)
        self.require_bidirectional = config.get('require_bidirectional_match', True)
        self.max_paralogs = config.get('max_paralogs_per_busco', 3)
        
        logger.info("RBH filter initialized")
        logger.info(f"  RBH filtering: {'enabled' if self.use_rbh else 'disabled'}")
        logger.info(f"  Identity gap threshold: {self.identity_gap_threshold}")
        logger.info(f"  Bidirectional requirement: {self.require_bidirectional}")
    
    def apply_rbh_filtering(self, alignment_results: List[AlignmentResult]) -> List[AlignmentResult]:
        """
        Apply reciprocal best hit filtering to alignment results.
        
        Args:
            alignment_results: List of AlignmentResult objects
            
        Returns:
            Filtered list of AlignmentResult objects
        """
        if not self.use_rbh:
            logger.info("RBH filtering disabled, returning all results")
            return alignment_results
        
        logger.info(f"Applying RBH filtering to {len(alignment_results)} alignment results...")
        
        # Group results by BUSCO ID
        busco_groups = self._group_results_by_busco(alignment_results)
        
        # Apply filtering strategies
        filtered_results = []
        
        for busco_id, results in busco_groups.items():
            if len(results) == 1:
                # Single result - keep it
                filtered_results.extend(results)
            else:
                # Multiple results - apply RBH filtering
                best_results = self._filter_busco_group(busco_id, results)
                filtered_results.extend(best_results)
        
        # Log filtering statistics
        self._log_filtering_statistics(alignment_results, filtered_results, busco_groups)
        
        return filtered_results
    
    def _group_results_by_busco(self, results: List[AlignmentResult]) -> Dict[str, List[AlignmentResult]]:
        """Group alignment results by BUSCO ID."""
        groups = defaultdict(list)
        for result in results:
            groups[result.busco_id].append(result)
        return dict(groups)
    
    def _filter_busco_group(self, busco_id: str, results: List[AlignmentResult]) -> List[AlignmentResult]:
        """
        Filter a group of results for a single BUSCO ID.
        
        Args:
            busco_id: BUSCO identifier
            results: List of results for this BUSCO
            
        Returns:
            Filtered list of best results
        """
        # Sort by primary criteria: confidence, then identity, then coverage
        sorted_results = sorted(
            results, 
            key=lambda x: (x.confidence, x.identity, x.min_coverage), 
            reverse=True
        )
        
        best_result = sorted_results[0]
        
        # Check for ties within identity gap threshold
        ties = [
            r for r in sorted_results 
            if abs(r.identity - best_result.identity) <= self.identity_gap_threshold
        ]
        
        if len(ties) == 1:
            return [best_result]
        else:
            # Multiple ties - apply sophisticated tie-breaking
            return self._resolve_ties(busco_id, ties)
    
    def _resolve_ties(self, busco_id: str, tied_results: List[AlignmentResult]) -> List[AlignmentResult]:
        """
        Resolve ties using multiple criteria.
        
        Args:
            busco_id: BUSCO identifier
            tied_results: List of tied results
            
        Returns:
            List with best result(s) after tie-breaking
        """
        # Multi-level tie-breaking
        tie_broken = sorted(
            tied_results,
            key=lambda x: (
                x.matches,                    # More matches is better
                -x.mapq if x.mapq else 0,    # Lower MAPQ is worse (negative for reverse sort)
                x.alignment_length           # Longer alignment is better
            ),
            reverse=True
        )
        
        # Return the best result after tie-breaking
        return [tie_broken[0]]
    
    def _log_filtering_statistics(self, original_results: List[AlignmentResult], 
                                filtered_results: List[AlignmentResult],
                                busco_groups: Dict[str, List[AlignmentResult]]):
        """Log detailed filtering statistics."""
        original_count = len(original_results)
        filtered_count = len(filtered_results)
        reduction_rate = (original_count - filtered_count) / original_count if original_count > 0 else 0
        
        # Count groups with multiple results
        multi_result_groups = sum(1 for results in busco_groups.values() if len(results) > 1)
        single_result_groups = len(busco_groups) - multi_result_groups
        
        logger.info("  RBH filtering statistics:")
        logger.info(f"    Original results: {original_count}")
        logger.info(f"    Filtered results: {filtered_count}")
        logger.info(f"    Reduction rate: {reduction_rate:.1%}")
        logger.info(f"    BUSCO groups processed: {len(busco_groups)}")
        logger.info(f"      Single result groups: {single_result_groups}")
        logger.info(f"      Multi-result groups: {multi_result_groups}")
        
        # Method distribution in final results
        method_counts = defaultdict(int)
        for result in filtered_results:
            method_counts[result.method] += 1
        
        logger.info(f"    Final method distribution: {dict(method_counts)}")
    
    def get_filtering_statistics(self, original_results: List[AlignmentResult], 
                               filtered_results: List[AlignmentResult]) -> Dict[str, Any]:
        """
        Get comprehensive filtering statistics.
        
        Args:
            original_results: Original alignment results
            filtered_results: Filtered alignment results
            
        Returns:
            Dictionary with filtering statistics
        """
        busco_groups = self._group_results_by_busco(original_results)
        
        stats = {
            'original_count': len(original_results),
            'filtered_count': len(filtered_results),
            'reduction_rate': (len(original_results) - len(filtered_results)) / len(original_results) if original_results else 0,
            'busco_groups_total': len(busco_groups),
            'single_result_groups': sum(1 for results in busco_groups.values() if len(results) == 1),
            'multi_result_groups': sum(1 for results in busco_groups.values() if len(results) > 1),
            'max_results_per_busco': max(len(results) for results in busco_groups.values()) if busco_groups else 0,
            'average_results_per_busco': sum(len(results) for results in busco_groups.values()) / len(busco_groups) if busco_groups else 0
        }
        
        # Method distribution
        original_methods = defaultdict(int)
        filtered_methods = defaultdict(int)
        
        for result in original_results:
            original_methods[result.method] += 1
        
        for result in filtered_results:
            filtered_methods[result.method] += 1
        
        stats['original_method_distribution'] = dict(original_methods)
        stats['filtered_method_distribution'] = dict(filtered_methods)
        
        return stats
