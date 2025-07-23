"""
Genome Inversion Analyzer
A modular toolkit for detecting and analyzing chromosomal inversions
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__description__ = "A comprehensive toolkit for genome synteny and inversion analysis"

# Import main modules for easy access
from . import config
from . import utils
from . import core
from . import visualization

# Import main functions for direct use
from .main import run_complete_enhanced_analysis_with_hybrid

__all__ = [
    'config',
    'utils', 
    'core',
    'visualization',
    'run_complete_enhanced_analysis_with_hybrid'
]