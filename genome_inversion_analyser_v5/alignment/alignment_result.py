
# =============================================================================
# Alignment Result Data Structure (alignment_result.py)
# =============================================================================

"""
Data structures for alignment results with comprehensive metadata.
"""

from dataclasses import dataclass
from typing import Optional, Dict, Any
import pandas as pd

@dataclass
class AlignmentResult:
    """
    Represents a single alignment result with comprehensive metadata.
    Standardizes results from different alignment methods.
    """
    # Core identification
    pair_index: int
    busco_id: str
    
    # Alignment quality metrics
    identity: float
    query_coverage: float
    target_coverage: float
    min_coverage: float
    alignment_length: int
    matches: int
    
    # Method and confidence
    method: str
    confidence: float = 0.0
    similarity: float = 0.0
    
    # Method-specific data
    score: Optional[float] = None
    mapq: Optional[int] = None
    strand: Optional[str] = None
    
    # Validation metadata
    validation_method: Optional[str] = None
    alternative_confidence: Optional[float] = None
    
    def __post_init__(self):
        """Calculate derived metrics after initialization."""
        # Calculate combined similarity if not set
        if self.similarity == 0.0:
            self.similarity = self.identity * self.min_coverage
        
        # Calculate basic confidence if not set
        if self.confidence == 0.0:
            self.confidence = (self.identity + self.min_coverage) / 2
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary for DataFrame creation."""
        return {
            'pair_index': self.pair_index,
            'busco_id': self.busco_id,
            'identity': self.identity,
            'query_coverage': self.query_coverage,
            'target_coverage': self.target_coverage,
            'min_coverage': self.min_coverage,
            'alignment_length': self.alignment_length,
            'matches': self.matches,
            'method': self.method,
            'confidence': self.confidence,
            'similarity': self.similarity,
            'score': self.score,
            'mapq': self.mapq,
            'strand': self.strand,
            'validation_method': self.validation_method,
            'alternative_confidence': self.alternative_confidence
        }
    
    @classmethod
    def from_minimap2(cls, pair_index: int, busco_id: str, paf_fields: Dict) -> 'AlignmentResult':
        """
        Create AlignmentResult from Minimap2 PAF output.
        
        Args:
            pair_index: Index of sequence pair
            busco_id: BUSCO identifier
            paf_fields: Dictionary with PAF field values
            
        Returns:
            AlignmentResult instance
        """
        return cls(
            pair_index=pair_index,
            busco_id=busco_id,
            identity=paf_fields['identity'],
            query_coverage=paf_fields['query_coverage'],
            target_coverage=paf_fields['target_coverage'],
            min_coverage=min(paf_fields['query_coverage'], paf_fields['target_coverage']),
            alignment_length=paf_fields['alignment_length'],
            matches=paf_fields['matches'],
            method='minimap2',
            mapq=paf_fields.get('mapq', 0),
            strand=paf_fields.get('strand', '+')
        )
    
    @classmethod
    def from_biopython(cls, pair_index: int, busco_id: str, bio_result: Dict) -> 'AlignmentResult':
        """
        Create AlignmentResult from Biopython alignment.
        
        Args:
            pair_index: Index of sequence pair
            busco_id: BUSCO identifier
            bio_result: Dictionary with Biopython alignment data
            
        Returns:
            AlignmentResult instance
        """
        return cls(
            pair_index=pair_index,
            busco_id=busco_id,
            identity=bio_result['identity'],
            query_coverage=bio_result['query_coverage'],
            target_coverage=bio_result['target_coverage'],
            min_coverage=min(bio_result['query_coverage'], bio_result['target_coverage']),
            alignment_length=bio_result['alignment_length'],
            matches=bio_result['matches'],
            method='biopython',
            score=bio_result.get('score', 0.0)
        )
    
    @classmethod
    def from_fallback(cls, pair_index: int, busco_id: str, similarity: float) -> 'AlignmentResult':
        """
        Create AlignmentResult from fallback similarity method.
        
        Args:
            pair_index: Index of sequence pair
            busco_id: BUSCO identifier
            similarity: Similarity score from difflib or similar
            
        Returns:
            AlignmentResult instance
        """
        return cls(
            pair_index=pair_index,
            busco_id=busco_id,
            identity=similarity,
            query_coverage=1.0,
            target_coverage=1.0,
            min_coverage=1.0,
            alignment_length=100,  # Approximate
            matches=int(similarity * 100),
            method='difflib_fallback',
            score=similarity
        )