
# =============================================================================
# Alignment Scorer (scorer.py)
# =============================================================================

"""
Comprehensive alignment scoring and normalization system.
Normalizes scores across different alignment methods and calculates unified confidence metrics.
"""

import numpy as np
from typing import List, Dict, Any
from ..alignment.alignment_result import AlignmentResult
from ..logger import get_logger

logger = get_logger()

class AlignmentScorer:
    """
    Comprehensive alignment scorer that normalizes results from different methods
    and calculates unified confidence metrics.
    """
    
    def __init__(self, config):
        """
        Initialize alignment scorer.
        
        Args:
            config: Configuration object with scoring parameters
        """
        self.config = config
        self.normalize_scores = config.get('normalize_scores', True)
        self.confidence_weights = config.get('confidence_weighting', {
            'identity': 0.4,
            'coverage': 0.3,
            'length_ratio': 0.2,
            'score_quality': 0.1
        })
        
        logger.info("Alignment scorer initialized")
        logger.info(f"  Score normalization: {'enabled' if self.normalize_scores else 'disabled'}")
        logger.info(f"  Confidence weights: {self.confidence_weights}")
    
    def normalize_alignment_scores(self, alignment_results: List[AlignmentResult]) -> List[AlignmentResult]:
        """
        Normalize alignment scores from different methods to a unified scale.
        
        Args:
            alignment_results: List of AlignmentResult objects
            
        Returns:
            List of AlignmentResult objects with normalized scores
        """
        if not self.normalize_scores:
            return alignment_results
        
        logger.info(f"Normalizing scores for {len(alignment_results)} alignment results...")
        
        normalized_results = []
        
        for result in alignment_results:
            # Calculate unified confidence score
            confidence = self._calculate_unified_confidence(result)
            
            # Calculate combined similarity metric
            similarity = self._calculate_similarity_metric(result)
            
            # Update result with normalized scores
            result.confidence = confidence
            result.similarity = similarity
            
            normalized_results.append(result)
        
        # Log normalization statistics
        self._log_normalization_statistics(normalized_results)
        
        return normalized_results
    
    def _calculate_unified_confidence(self, result: AlignmentResult) -> float:
        """
        Calculate unified confidence score using weighted metrics.
        
        Args:
            result: AlignmentResult object
            
        Returns:
            Unified confidence score (0.0 to 1.0)
        """
        # Identity component
        identity_score = result.identity
        
        # Coverage component
        coverage_score = result.min_coverage
        
        # Length ratio component (if available from sequence data)
        length_ratio = 1.0  # Default - could be calculated from original sequences
        
        # Score quality component (method-dependent normalization)
        score_quality = self._normalize_method_score(result)
        
        # Calculate weighted confidence
        confidence = (
            self.confidence_weights['identity'] * identity_score +
            self.confidence_weights['coverage'] * coverage_score +
            self.confidence_weights['length_ratio'] * length_ratio +
            self.confidence_weights['score_quality'] * score_quality
        )
        
        return max(0.0, min(1.0, confidence))  # Clamp to [0,1]
    
    def _normalize_method_score(self, result: AlignmentResult) -> float:
        """
        Normalize method-specific scores to 0-1 scale.
        
        Args:
            result: AlignmentResult object
            
        Returns:
            Normalized score quality (0.0 to 1.0)
        """
        if result.method == 'minimap2':
            # MAPQ is typically 0-60, higher is better
            mapq = result.mapq if result.mapq is not None else 0
            return min(1.0, mapq / 60.0)
        
        elif result.method == 'biopython':
            # Biopython scores vary widely, use heuristic normalization
            score = result.score if result.score is not None else 0
            # Rough normalization: assume scores range from -100 to +200
            normalized = (score + 100) / 300.0
            return max(0.0, min(1.0, normalized))
        
        else:  # fallback methods
            # Use identity as proxy for score quality
            return result.identity
    
    def _calculate_similarity_metric(self, result: AlignmentResult) -> float:
        """
        Calculate combined similarity metric.
        
        Args:
            result: AlignmentResult object
            
        Returns:
            Combined similarity score (0.0 to 1.0)
        """
        # Combine identity and coverage with emphasis on identity
        return result.identity * 0.7 + result.min_coverage * 0.3
    
    def _log_normalization_statistics(self, results: List[AlignmentResult]):
        """Log statistics about score normalization."""
        if not results:
            return
        
        # Calculate statistics by method
        method_stats = {}
        for result in results:
            method = result.method
            if method not in method_stats:
                method_stats[method] = {
                    'count': 0,
                    'identities': [],
                    'confidences': [],
                    'similarities': []
                }
            
            stats = method_stats[method]
            stats['count'] += 1
            stats['identities'].append(result.identity)
            stats['confidences'].append(result.confidence)
            stats['similarities'].append(result.similarity)
        
        # Log summary statistics
        logger.info("  Score normalization statistics:")
        for method, stats in method_stats.items():
            avg_identity = np.mean(stats['identities'])
            avg_confidence = np.mean(stats['confidences'])
            avg_similarity = np.mean(stats['similarities'])
            
            logger.info(f"    {method} ({stats['count']} results):")
            logger.info(f"      Avg identity: {avg_identity:.3f}")
            logger.info(f"      Avg confidence: {avg_confidence:.3f}")
            logger.info(f"      Avg similarity: {avg_similarity:.3f}")
    
    def get_score_statistics(self, results: List[AlignmentResult]) -> Dict[str, Any]:
        """
        Get comprehensive scoring statistics.
        
        Args:
            results: List of AlignmentResult objects
            
        Returns:
            Dictionary with detailed statistics
        """
        if not results:
            return {}
        
        # Overall statistics
        identities = [r.identity for r in results]
        confidences = [r.confidence for r in results]
        similarities = [r.similarity for r in results]
        
        stats = {
            'total_results': len(results),
            'identity_stats': {
                'mean': np.mean(identities),
                'std': np.std(identities),
                'min': np.min(identities),
                'max': np.max(identities),
                'median': np.median(identities)
            },
            'confidence_stats': {
                'mean': np.mean(confidences),
                'std': np.std(confidences),
                'min': np.min(confidences),
                'max': np.max(confidences),
                'median': np.median(confidences)
            },
            'similarity_stats': {
                'mean': np.mean(similarities),
                'std': np.std(similarities),
                'min': np.min(similarities),
                'max': np.max(similarities),
                'median': np.median(similarities)
            },
            'method_distribution': {}
        }
        
        # Method-specific statistics
        for result in results:
            method = result.method
            if method not in stats['method_distribution']:
                stats['method_distribution'][method] = 0
            stats['method_distribution'][method] += 1
        
        return stats
