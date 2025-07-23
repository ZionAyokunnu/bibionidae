
# =============================================================================
# Confidence Calculator (confidence_calculator.py)
# =============================================================================

"""
Advanced confidence calculation system for ortholog relationships.
Integrates multiple evidence sources to provide robust confidence scores.
"""

import numpy as np
from typing import List, Dict, Any, Optional
from ..alignment.alignment_result import AlignmentResult
from ..logger import get_logger

logger = get_logger()

class ConfidenceCalculator:
    """
    Advanced confidence calculator that integrates multiple evidence sources
    to provide robust confidence scores for ortholog relationships.
    """
    
    def __init__(self, config):
        """
        Initialize confidence calculator.
        
        Args:
            config: Configuration object with confidence parameters
        """
        self.config = config
        self.enable_confidence_scoring = config.get('enable_ortholog_confidence_scoring', True)
        
        # Confidence weighting scheme
        self.weights = config.get('confidence_weighting', {
            'sequence_similarity': 0.4,
            'alignment_quality': 0.3,
            'method_reliability': 0.2,
            'contextual_evidence': 0.1
        })
        
        # Method reliability scores
        self.method_reliability = {
            'minimap2': 0.85,      # Fast but less precise
            'biopython': 0.95,     # Slower but more precise
            'difflib_fallback': 0.6  # Least reliable
        }
        
        logger.info("Confidence calculator initialized")
        logger.info(f"  Confidence scoring: {'enabled' if self.enable_confidence_scoring else 'disabled'}")
        logger.info(f"  Weighting scheme: {self.weights}")
    
    def calculate_ortholog_confidence(self, result: AlignmentResult, 
                                    context: Optional[Dict] = None) -> float:
        """
        Calculate comprehensive confidence score for an ortholog relationship.
        
        Args:
            result: AlignmentResult object
            context: Optional contextual information
            
        Returns:
            Confidence score (0.0 to 1.0)
        """
        if not self.enable_confidence_scoring:
            return result.confidence if result.confidence > 0 else 0.5
        
        # Calculate individual components
        sequence_score = self._calculate_sequence_similarity_score(result)
        alignment_score = self._calculate_alignment_quality_score(result)
        method_score = self._calculate_method_reliability_score(result)
        context_score = self._calculate_contextual_evidence_score(result, context)
        
        # Calculate weighted confidence
        confidence = (
            self.weights['sequence_similarity'] * sequence_score +
            self.weights['alignment_quality'] * alignment_score +
            self.weights['method_reliability'] * method_score +
            self.weights['contextual_evidence'] * context_score
        )
        
        return max(0.0, min(1.0, confidence))
    
    def _calculate_sequence_similarity_score(self, result: AlignmentResult) -> float:
        """Calculate sequence similarity component of confidence."""
        # Base on identity with coverage adjustment
        identity_score = result.identity
        coverage_penalty = 1.0 - (1.0 - result.min_coverage) * 0.5  # Penalize low coverage
        
        return identity_score * coverage_penalty
    
    def _calculate_alignment_quality_score(self, result: AlignmentResult) -> float:
        """Calculate alignment quality component of confidence."""
        components = []
        
        # Coverage quality (prefer high, balanced coverage)
        coverage_balance = 1.0 - abs(result.query_coverage - result.target_coverage)
        coverage_magnitude = (result.query_coverage + result.target_coverage) / 2
        coverage_score = coverage_balance * coverage_magnitude
        components.append(coverage_score)
        
        # Alignment length quality
        if result.alignment_length > 0:
            # Prefer longer alignments, but with diminishing returns
            length_score = min(1.0, result.alignment_length / 500.0)  # Normalize to 500bp
            components.append(length_score)
        
        # Match quality
        if result.alignment_length > 0:
            match_rate = result.matches / result.alignment_length
            components.append(match_rate)
        
        return np.mean(components) if components else 0.5
    
    def _calculate_method_reliability_score(self, result: AlignmentResult) -> float:
        """Calculate method reliability component of confidence."""
        base_reliability = self.method_reliability.get(result.method, 0.5)
        
        # Adjust based on method-specific quality indicators
        if result.method == 'minimap2' and result.mapq is not None:
            # Higher MAPQ indicates better mapping quality
            mapq_bonus = min(0.15, result.mapq / 60.0 * 0.15)  # Max 15% bonus
            return min(1.0, base_reliability + mapq_bonus)
        
        elif result.method == 'biopython' and result.score is not None:
            # Higher scores indicate better alignment
            if result.score > 100:  # Good score threshold
                score_bonus = min(0.05, (result.score - 100) / 1000.0 * 0.05)  # Max 5% bonus
                return min(1.0, base_reliability + score_bonus)
        
        return base_reliability
    
    def _calculate_contextual_evidence_score(self, result: AlignmentResult, 
                                           context: Optional[Dict]) -> float:
        """Calculate contextual evidence component of confidence."""
        if not context:
            return 0.5  # Neutral score when no context available
        
        evidence_factors = []
        
        # Cross-validation evidence (if available)
        if 'alternative_confidence' in context:
            alt_conf = context['alternative_confidence']
            # High agreement between methods increases confidence
            agreement = 1.0 - abs(result.confidence - alt_conf)
            evidence_factors.append(agreement)
        
        # Synteny context (if available)
        if 'synteny_support' in context:
            synteny_score = context['synteny_support']
            evidence_factors.append(synteny_score)
        
        # Paralog context (lower confidence for high paralog counts)
        if 'paralog_count' in context:
            paralog_count = context['paralog_count']
            paralog_penalty = 1.0 / (1.0 + (paralog_count - 1) * 0.2)  # Penalty for multiple paralogs
            evidence_factors.append(paralog_penalty)
        
        return np.mean(evidence_factors) if evidence_factors else 0.5
    
    def calculate_batch_confidence(self, results: List[AlignmentResult], 
                                 contexts: Optional[List[Dict]] = None) -> List[float]:
        """
        Calculate confidence scores for a batch of results.
        
        Args:
            results: List of AlignmentResult objects
            contexts: Optional list of context dictionaries
            
        Returns:
            List of confidence scores
        """
        if not contexts:
            contexts = [None] * len(results)
        
        confidences = []
        for result, context in zip(results, contexts):
            confidence = self.calculate_ortholog_confidence(result, context)
            confidences.append(confidence)
        
        return confidences
    
    def get_confidence_statistics(self, results: List[AlignmentResult]) -> Dict[str, Any]:
        """
        Get comprehensive confidence statistics for a set of results.
        
        Args:
            results: List of AlignmentResult objects
            
        Returns:
            Dictionary with confidence statistics
        """
        if not results:
            return {}
        
        confidences = [r.confidence for r in results]
        
        stats = {
            'total_results': len(results),
            'confidence_stats': {
                'mean': np.mean(confidences),
                'std': np.std(confidences),
                'min': np.min(confidences),
                'max': np.max(confidences),
                'median': np.median(confidences),
                'q25': np.percentile(confidences, 25),
                'q75': np.percentile(confidences, 75)
            },
            'quality_distribution': {
                'high_confidence': sum(1 for c in confidences if c >= 0.8),
                'medium_confidence': sum(1 for c in confidences if 0.5 <= c < 0.8),
                'low_confidence': sum(1 for c in confidences if c < 0.5)
            },
            'method_confidence': {}
        }
        
        # Method-specific confidence statistics
        for result in results:
            method = result.method
            if method not in stats['method_confidence']:
                stats['method_confidence'][method] = []
            stats['method_confidence'][method].append(result.confidence)
        
        # Calculate method averages
        for method, confs in stats['method_confidence'].items():
            stats['method_confidence'][method] = {
                'count': len(confs),
                'mean': np.mean(confs),
                'std': np.std(confs)
            }
        
        return stats