"""
Genome Inversion Analyzer
A modular toolkit for detecting and analyzing chromosomal inversions
"""

__version__ = "1.0.0"
__author__ = "Zion Ayokunnu"
__description__ = "A comprehensive toolkit for genome synteny and inversion analysis"

# Import main modules for easy access
from . import config
from . import utils
from . import core
from . import visualization

# Try to import registry if available
try:
    from . import registry
except ImportError:
    registry = None

# Try to import main functions - handle gracefully if not available
try:
    from .main import run_complete_enhanced_analysis_with_hybrid
except ImportError:
    # Fallback for compatibility
    run_complete_enhanced_analysis_with_hybrid = None

__all__ = [
    'config',
    'utils', 
    'core',
    'visualization'
]

# Add registry to exports if available
if registry is not None:
    __all__.append('registry')

# Add main function if available  
if run_complete_enhanced_analysis_with_hybrid is not None:
    __all__.append('run_complete_enhanced_analysis_with_hybrid')