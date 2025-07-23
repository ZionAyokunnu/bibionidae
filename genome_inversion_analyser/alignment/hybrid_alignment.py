
# =============================================================================
# 6. Hybrid Alignment Engine (hybrid_alignment.py)
# =============================================================================

"""
Main hybrid alignment engine that orchestrates the entire alignment workflow.
Intelligently selects alignment methods and combines results.
"""

import hashlib
import pickle
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd

from ..logger import get_logger
from .partition import SequencePartitioner
from .biopython_runner import BioPythonAligner
from .minimap_runner import Minimap2Runner
from .alignment_result import AlignmentResult
from ..orthology.scorer import AlignmentScorer
from ..orthology.rbh_filter import ReciprocalBestHitFilter
from ..orthology.converter import OrthologConverter

logger = get_logger()

class HybridAlignmentEngine:
    """
    Main hybrid alignment engine that orchestrates the complete alignment workflow.
    Combines Minimap2 and Biopython for optimal speed and accuracy.
    """
    
    def __init__(self, config):
        """
        Initialize hybrid alignment engine.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.partitioner = SequencePartitioner(config)
        self.biopython_aligner = BioPythonAligner(config)
        self.minimap2_runner = Minimap2Runner(config)
        self.scorer = AlignmentScorer(config)
        self.rbh_filter = ReciprocalBestHitFilter(config)
        self.converter = OrthologConverter(config)
        
        # Results cache
        self.cache_enabled = config.get('alignment_cache_enabled', True)
        self.cache_dir = None
        if self.cache_enabled:
            self.cache_dir = Path(config.get('base_output_dir', 'results')) / 'cache'
            self.cache_dir.mkdir(exist_ok=True)
        
        logger.info("Hybrid alignment engine initialized")
        logger.info(f"  Strategy: {config.get('alignment_strategy', 'hybrid')}")
        logger.info(f"  Caching: {'enabled' if self.cache_enabled else 'disabled'}")
    
    def run_hybrid_alignment_analysis(self, first_busco_df: pd.DataFrame, 
                                    second_busco_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run complete hybrid alignment analysis pipeline.
        
        Args:
            first_busco_df: First genome BUSCO data
            second_busco_df: Second genome BUSCO data
            
        Returns:
            Tuple of (ortholog_df, paralog_df)
        """
        logger.info("Running hybrid alignment analysis...")
        
        # Create sequence pairs
        sequence_pairs = self._create_sequence_pairs(first_busco_df, second_busco_df)
        
        # Check cache first
        cache_key = self._generate_cache_key(first_busco_df, second_busco_df)
        cached_result = self._load_cached_results(cache_key)
        
        if cached_result:
            logger.info("  Using cached alignment results")
            ortholog_df, paralog_df = cached_result
        else:
            # Run alignment pipeline
            alignment_results = self._run_alignment_pipeline(sequence_pairs)
            
            # Process results
            ortholog_df, paralog_df = self._process_alignment_results(
                alignment_results, sequence_pairs
            )
            
            # Cache results
            self._cache_results((ortholog_df, paralog_df), cache_key)
        
        # Report final statistics
        self._report_final_statistics(ortholog_df, alignment_results if 'alignment_results' in locals() else [])
        
        return ortholog_df, paralog_df
    
    def _create_sequence_pairs(self, first_busco_df: pd.DataFrame, 
                             second_busco_df: pd.DataFrame) -> List[Dict]:
        """Create sequence pairs from BUSCO dataframes."""
        # Find common BUSCO genes
        first_buscos = {row['busco_id']: row for _, row in first_busco_df.iterrows()}
        second_buscos = {row['busco_id']: row for _, row in second_busco_df.iterrows()}
        common_buscos = set(first_buscos.keys()) & set(second_buscos.keys())
        
        logger.info(f"  Species 1 BUSCOs: {len(first_buscos)}")
        logger.info(f"  Species 2 BUSCOs: {len(second_buscos)}")
        logger.info(f"  Common BUSCOs: {len(common_buscos)}")
        
        # Create sequence pairs
        sequence_pairs = []
        for busco_id in common_buscos:
            pair = {
                'busco_id': busco_id,
                'first_gene': first_buscos[busco_id],
                'second_gene': second_buscos[busco_id]
            }
            sequence_pairs.append(pair)
        
        return sequence_pairs
    
    def _run_alignment_pipeline(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """Run the complete alignment pipeline based on strategy."""
        strategy = self.config.get('alignment_strategy', 'hybrid')
        
        if strategy == 'hybrid':
            return self._run_hybrid_strategy(sequence_pairs)
        elif strategy == 'minimap2':
            return self._run_minimap2_only(sequence_pairs)
        else:  # biopython
            return self._run_biopython_only(sequence_pairs)
    
    def _run_hybrid_strategy(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """Run hybrid alignment strategy with intelligent method selection."""
        logger.info("  Using hybrid alignment strategy")
        
        # Partition sequences by length
        partitions = self.partitioner.partition_sequence_pairs(sequence_pairs)
        
        all_results = []
        
        # Process short sequences with Biopython
        if partitions['short_pairs']:
            logger.info(f"  Running Biopython on {len(partitions['short_pairs'])} short sequences...")
            short_results = self.biopython_aligner.align_sequences(partitions['short_pairs'])
            all_results.extend(short_results)
        
        # Process long sequences with Minimap2
        if partitions['long_pairs']:
            logger.info(f"  Running Minimap2 on {len(partitions['long_pairs'])} long sequences...")
            long_results = self.minimap2_runner.align_sequences(partitions['long_pairs'])
            all_results.extend(long_results)
        
        # Process buffer zone sequences
        if partitions['buffer_pairs']:
            buffer_method = self.config.get('buffer_zone_method', 'dual')
            logger.info(f"  Processing {len(partitions['buffer_pairs'])} buffer zone sequences with {buffer_method}...")
            
            if buffer_method == 'dual' and self.config.get('cross_validate_buffer_zone', True):
                buffer_results = self._process_buffer_zone_dual(partitions['buffer_pairs'])
            elif buffer_method == 'minimap2':
                buffer_results = self.minimap2_runner.align_sequences(partitions['buffer_pairs'])
            else:  # biopython
                buffer_results = self.biopython_aligner.align_sequences(partitions['buffer_pairs'])
            
            all_results.extend(buffer_results)
        
        # Process mixed length sequences with Biopython (safer)
        if partitions['mixed_pairs']:
            logger.info(f"  Running Biopython on {len(partitions['mixed_pairs'])} mixed-length sequences...")
            mixed_results = self.biopython_aligner.align_sequences(partitions['mixed_pairs'])
            all_results.extend(mixed_results)
        
        return all_results
    
    def _run_minimap2_only(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """Run Minimap2-only alignment strategy."""
        logger.info("  Using Minimap2-only alignment strategy")
        return self.minimap2_runner.align_sequences(sequence_pairs)
    
    def _run_biopython_only(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """Run Biopython-only alignment strategy."""
        logger.info("  Using Biopython-only alignment strategy")
        return self.biopython_aligner.align_sequences(sequence_pairs)
    
    def _process_buffer_zone_dual(self, buffer_pairs: List[Dict]) -> List[AlignmentResult]:
        """Process buffer zone sequences with both methods and select best results."""
        logger.info("    Running dual validation on buffer zone sequences...")
        
        # Run both methods
        bio_results = self.biopython_aligner.align_sequences(buffer_pairs)
        mm2_results = self.minimap2_runner.align_sequences(buffer_pairs)
        
        # Select best results
        return self._select_best_buffer_results(bio_results, mm2_results)
    
    def _select_best_buffer_results(self, bio_results: List[AlignmentResult], 
                                  mm2_results: List[AlignmentResult]) -> List[AlignmentResult]:
        """Select best results when both methods are run on buffer zone."""
        # Create mappings by BUSCO ID
        bio_map = {r.busco_id: r for r in bio_results}
        mm2_map = {r.busco_id: r for r in mm2_results}
        
        best_results = []
        all_buscos = set(bio_map.keys()) | set(mm2_map.keys())
        
        for busco_id in all_buscos:
            bio_result = bio_map.get(busco_id)
            mm2_result = mm2_map.get(busco_id)
            
            if bio_result and mm2_result:
                # Both methods succeeded - choose best based on confidence
                if bio_result.confidence >= mm2_result.confidence:
                    bio_result.validation_method = 'dual_validated'
                    bio_result.alternative_confidence = mm2_result.confidence
                    best_results.append(bio_result)
                else:
                    mm2_result.validation_method = 'dual_validated'
                    mm2_result.alternative_confidence = bio_result.confidence
                    best_results.append(mm2_result)
            elif bio_result:
                bio_result.validation_method = 'biopython_only'
                best_results.append(bio_result)
            elif mm2_result:
                mm2_result.validation_method = 'minimap2_only'
                best_results.append(mm2_result)
        
        return best_results
    
    def _process_alignment_results(self, alignment_results: List[AlignmentResult], 
                                 sequence_pairs: List[Dict]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Process alignment results into ortholog and paralog dataframes."""
        logger.info("  Processing alignment results...")
        
        # Normalize scores and calculate confidence
        logger.info("    Normalizing scores and calculating confidence...")
        normalized_results = self.scorer.normalize_alignment_scores(alignment_results)
        
        # Apply reciprocal best hit filtering
        if self.config.get('use_reciprocal_best_hits', True):
            logger.info("    Applying reciprocal best hit filtering...")
            filtered_results = self.rbh_filter.apply_rbh_filtering(normalized_results)
        else:
            filtered_results = normalized_results
        
        # Convert to ortholog pairs format
        logger.info("    Converting to ortholog pairs format...")
        ortholog_df = self.converter.convert_to_ortholog_pairs(filtered_results, sequence_pairs)
        
        # Create empty paralog dataframe for now (could be expanded)
        paralog_df = pd.DataFrame()
        
        return ortholog_df, paralog_df
    
    def _generate_cache_key(self, first_busco_df: pd.DataFrame, second_busco_df: pd.DataFrame) -> str:
        """Generate cache key based on input data and configuration."""
        key_config = {
            'alignment_strategy': self.config.get('alignment_strategy'),
            'short_sequence_threshold': self.config.get('short_sequence_threshold'),
            'minimap2_kmer_size': self.config.get('minimap2_kmer_size'),
            'base_similarity_threshold': self.config.get('base_similarity_threshold')
        }
        
        # Create data fingerprint
        first_hash = hashlib.md5(str(sorted(first_busco_df['busco_id'].tolist())).encode()).hexdigest()[:8]
        second_hash = hashlib.md5(str(sorted(second_busco_df['busco_id'].tolist())).encode()).hexdigest()[:8]
        config_hash = hashlib.md5(str(sorted(key_config.items())).encode()).hexdigest()[:8]
        
        return f"{first_hash}_{second_hash}_{config_hash}"
    
    def _load_cached_results(self, cache_key: str) -> Optional[Tuple[pd.DataFrame, pd.DataFrame]]:
        """Load cached alignment results if available."""
        if not self.cache_enabled or not self.cache_dir:
            return None
        
        cache_file = self.cache_dir / f"alignment_cache_{cache_key}.pkl"
        
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    results = pickle.load(f)
                logger.info(f"    Loaded cached results from {cache_file}")
                return results
            except Exception as e:
                logger.warning(f"    Failed to load cached results: {e}")
        
        return None
    
    def _cache_results(self, results: Tuple[pd.DataFrame, pd.DataFrame], cache_key: str):
        """Cache alignment results for future use."""
        if not self.cache_enabled or not self.cache_dir:
            return
        
        cache_file = self.cache_dir / f"alignment_cache_{cache_key}.pkl"
        
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(results, f)
            logger.info(f"    Cached results to {cache_file}")
        except Exception as e:
            logger.warning(f"    Failed to cache results: {e}")
    
    def _report_final_statistics(self, ortholog_df: pd.DataFrame, alignment_results: List[AlignmentResult]):
        """Report comprehensive final statistics."""
        logger.info(f"  Alignment results:")
        logger.info(f"    Final ortholog pairs: {len(ortholog_df)}")
        
        if alignment_results:
            logger.info(f"    Total alignments attempted: {len(alignment_results)}")
            
            # Method usage statistics
            method_counts = {}
            for result in alignment_results:
                method = result.method
                method_counts[method] = method_counts.get(method, 0) + 1
            
            logger.info(f"    Methods used: {method_counts}")
            
            # Quality statistics
            avg_identity = sum(r.identity for r in alignment_results) / len(alignment_results)
            avg_confidence = sum(r.confidence for r in alignment_results) / len(alignment_results)
            logger.info(f"    Average identity: {avg_identity:.3f}")
            logger.info(f"    Average confidence: {avg_confidence:.3f}")
        
        if len(ortholog_df) > 0:
            # Ortholog-specific statistics
            if 'similarity' in ortholog_df.columns:
                logger.info(f"    Final average similarity: {ortholog_df['similarity'].mean():.3f}")
            if 'alignment_method' in ortholog_df.columns:
                final_methods = ortholog_df['alignment_method'].value_counts()
                logger.info(f"    Final method distribution: {final_methods.to_dict()}")
