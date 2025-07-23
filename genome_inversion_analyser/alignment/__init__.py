
# =============================================================================
# Alignment Package (__init__.py)
# =============================================================================

"""
Hybrid alignment system for genome analysis.
Combines Minimap2 (fast, for long sequences) with Biopython (precise, for short sequences).
"""

from .partition import SequencePartitioner
from .biopython_runner import BioPythonAligner
from .minimap_runner import Minimap2Runner
from .hybrid_alignment import HybridAlignmentEngine
from .alignment_result import AlignmentResult

__all__ = [
    'SequencePartitioner',
    'BioPythonAligner', 
    'Minimap2Runner',
    'HybridAlignmentEngine',
    'AlignmentResult'
]