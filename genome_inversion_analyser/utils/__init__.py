# =============================================================================
# Utility Functions (utils/__init__.py)
# =============================================================================

"""
Common utility functions used across the analysis pipeline.
"""

import pandas as pd
import numpy as np
import hashlib
from pathlib import Path

def standardize_sequence_id(id_str):
    """Standardize sequence IDs to ensure consistency across tools."""
    if pd.isna(id_str) or not id_str:
        return None
    
    # Remove common prefixes and clean up
    clean_id = str(id_str).strip()
    clean_id = clean_id.replace('>', '').split()[0]  # Take only first part
    clean_id = clean_id.replace('|', '_').replace(':', '_')  # Replace problematic chars
    
    return clean_id

def create_output_directory(config):
    """Create output directory structure."""
    base_dir = Path(config.get('base_output_dir', 'enhanced_results'))
    base_dir.mkdir(exist_ok=True)
    
    subdirs = ['plots', 'data', 'reports', 'debug', 'cache']
    for subdir in subdirs:
        (base_dir / subdir).mkdir(exist_ok=True)
    
    return base_dir

def generate_cache_key(*args):
    """Generate a cache key based on input arguments."""
    # Create hash of key arguments
    combined_str = str(sorted(str(arg) for arg in args))
    return hashlib.md5(combined_str.encode()).hexdigest()[:16]

def safe_divide(numerator, denominator, default=0.0):
    """Safely divide two numbers, returning default if denominator is zero."""
    try:
        return numerator / denominator if denominator != 0 else default
    except (TypeError, ZeroDivisionError):
        return default

def calculate_statistics(values, include_std=True):
    """Calculate basic statistics for a list of values."""
    if not values:
        return {}
    
    values = np.array(values)
    stats = {
        'count': len(values),
        'mean': np.mean(values),
        'median': np.median(values),
        'min': np.min(values),
        'max': np.max(values)
    }
    
    if include_std:
        stats['std'] = np.std(values)
    
    return stats
