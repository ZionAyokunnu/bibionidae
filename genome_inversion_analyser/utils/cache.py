"""
Caching utilities for the Genome Inversion Analyzer
Handles caching of alignment results to avoid recomputation
"""

import hashlib
import pickle
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def generate_cache_key(first_busco_df, second_busco_df, config):
    """Generate a cache key based on input data and configuration"""
    # Create hash of key configuration parameters and data
    key_config = {
        'alignment_strategy': config.get('alignment_strategy'),
        'short_sequence_threshold': config.get('short_sequence_threshold'),
        'minimap2_kmer_size': config.get('minimap2_kmer_size'),
        'base_similarity_threshold': config.get('base_similarity_threshold')
    }
    
    # Add data fingerprint
    first_hash = hashlib.md5(str(sorted(first_busco_df['busco_id'].tolist())).encode()).hexdigest()[:8]
    second_hash = hashlib.md5(str(sorted(second_busco_df['busco_id'].tolist())).encode()).hexdigest()[:8]
    config_hash = hashlib.md5(str(sorted(key_config.items())).encode()).hexdigest()[:8]
    
    return f"{first_hash}_{second_hash}_{config_hash}"


def cache_alignment_results(results, cache_key, config):
    """Cache alignment results to avoid recomputation"""
    if not config.get('alignment_cache_enabled', True):
        return
    
    cache_dir = Path(config.get('base_output_dir', 'v4/enhanced_results')) / 'cache'
    cache_dir.mkdir(exist_ok=True)
    
    cache_file = cache_dir / f"alignment_cache_{cache_key}.pkl"
    
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(results, f)
        logger.info(f"  Cached alignment results to {cache_file}")
    except Exception as e:
        logger.warning(f"  Failed to cache results: {e}")


def load_cached_alignment_results(cache_key, config):
    """Load cached alignment results if available"""
    if not config.get('alignment_cache_enabled', True):
        return None
    
    cache_dir = Path(config.get('base_output_dir', 'v4/enhanced_results')) / 'cache'
    cache_file = cache_dir / f"alignment_cache_{cache_key}.pkl"
    
    if cache_file.exists():
        try:
            with open(cache_file, 'rb') as f:
                results = pickle.load(f)
            logger.info(f"  Loaded cached alignment results from {cache_file}")
            return results
        except Exception as e:
            logger.warning(f"  Failed to load cached results: {e}")
    
    return None