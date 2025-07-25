"""
Alignment module for the Genome Inversion Analyzer
Handles hybrid Minimap2 + Biopython alignment system
"""

import pandas as pd
import numpy as np
import logging
import subprocess
import tempfile
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from difflib import SequenceMatcher
from Bio.Align import PairwiseAligner

# Import from utils module
from ..utils import standardize_sequence_id

logger = logging.getLogger(__name__)


def partition_sequences_by_length(ortholog_pairs, config):
    """Partition sequence pairs by length for optimal alignment method selection"""
    short_threshold = config.get('short_sequence_threshold', 500)
    long_threshold = config.get('long_sequence_threshold', 1500)
    
    partitions = {
        'short_pairs': [],      # Use Biopython
        'long_pairs': [],       # Use Minimap2  
        'buffer_pairs': [],     # Use both or fallback method
        'mixed_pairs': []       # One short, one long - needs special handling
    }
    
    for pair in ortholog_pairs:
        len1 = pair['first_gene']['gene_length']
        len2 = pair['second_gene']['gene_length']
        
        max_len = max(len1, len2)
        min_len = min(len1, len2)
        
        if max_len <= short_threshold:
            partitions['short_pairs'].append(pair)
        elif min_len >= long_threshold:
            partitions['long_pairs'].append(pair)
        elif short_threshold <= max_len <= long_threshold:
            partitions['buffer_pairs'].append(pair)
        else:
            # One sequence much longer than the other
            partitions['mixed_pairs'].append(pair)
    
    # Log partition statistics
    if config.get('detailed_alignment_logging', False):
        logger.info(f"  Sequence partitioning:")
        logger.info(f"    Short pairs (≤{short_threshold}bp): {len(partitions['short_pairs'])}")
        logger.info(f"    Long pairs (≥{long_threshold}bp): {len(partitions['long_pairs'])}")
        logger.info(f"    Buffer zone pairs: {len(partitions['buffer_pairs'])}")
        logger.info(f"    Mixed length pairs: {len(partitions['mixed_pairs'])}")
    
    return partitions


def create_minimap2_fasta(sequence_pairs, temp_dir):
    """Create FASTA files for minimap2 alignment"""
    query_file = temp_dir / "queries.fasta"
    target_file = temp_dir / "targets.fasta"
    
    # Create sequence mappings
    query_map = {}
    target_map = {}
    
    with open(query_file, 'w') as qf, open(target_file, 'w') as tf:
        for i, pair in enumerate(sequence_pairs):
            # Standardize IDs
            query_id = f"query_{i}_{pair['first_gene']['busco_id']}"
            target_id = f"target_{i}_{pair['second_gene']['busco_id']}"
            
            # Write sequences
            qf.write(f">{query_id}\n{pair['first_gene']['gene_sequence']}\n")
            tf.write(f">{target_id}\n{pair['second_gene']['gene_sequence']}\n")
            
            # Store mappings
            query_map[query_id] = i
            target_map[target_id] = i
    
    return query_file, target_file, query_map, target_map


def run_minimap2_alignment(sequence_pairs, config):
    """Run minimap2 alignment on sequence pairs - FIXED VERSION"""
    if not sequence_pairs:
        return []
    
    results = []
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create FASTA files
            query_file, target_file, query_map, target_map = create_minimap2_fasta(sequence_pairs, temp_path)
            output_file = temp_path / "alignments.paf"
            
            # Build minimap2 command - FIXED
            cmd = ['minimap2']
            
            # Add preset (fix: use -x instead of --preset)
            preset = config.get('minimap2_preset', '--sr')
            if preset.startswith('--'):
                preset = preset[2:]  # Remove -- prefix
            cmd.extend(['-x', preset])
            
            # Add other parameters
            cmd.extend(['-k', str(config.get('minimap2_kmer_size', 13))])
            cmd.extend(['-t', str(config.get('minimap2_threads', 4))])
            
            # Fix: Use -s for minimum score, not --score-N
            min_score = config.get('minimap2_min_score', 100)
            cmd.extend(['-s', str(min_score)])
            
            # Add extra flags (but parse them properly)
            extra_flags = config.get('minimap2_extra_flags', '-c --cs')
            if extra_flags:
                # Split and add each flag separately
                for flag in extra_flags.split():
                    if flag.strip():
                        cmd.append(flag.strip())
            
            # Add input files
            cmd.append(str(target_file))
            cmd.append(str(query_file))
            
            # Log the actual command for debugging
            logger.info(f"    Running minimap2 command: {' '.join(cmd)}")
            
            # Run minimap2
            with open(output_file, 'w') as out_f:
                process = subprocess.run(
                    cmd, 
                    stdout=out_f, 
                    stderr=subprocess.PIPE, 
                    timeout=config.get('timeout_per_alignment', 30) * len(sequence_pairs),
                    text=True
                )
            
            if process.returncode != 0:
                logger.warning(f"Minimap2 failed: {process.stderr}")
                return []
            
            # Parse PAF output
            results = parse_minimap2_paf(output_file, sequence_pairs, query_map, target_map, config)
            
    except subprocess.TimeoutExpired:
        logger.error("Minimap2 alignment timed out")
        return []
    except FileNotFoundError:
        logger.error("Minimap2 not found in PATH. Install minimap2 or use 'alignment_strategy': 'biopython'")
        return []
    except Exception as e:
        logger.error(f"Minimap2 alignment failed: {e}")
        return []
    
    return results


def parse_minimap2_paf(paf_file, sequence_pairs, query_map, target_map, config):
    """Parse minimap2 PAF output and extract alignment statistics"""
    results = []
    min_identity = config.get('minimap2_min_identity', 0.7)
    
    try:
        with open(paf_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:
                        # Parse PAF fields
                        query_name = fields[0]
                        query_len = int(fields[1])
                        query_start = int(fields[2])
                        query_end = int(fields[3])
                        strand = fields[4]
                        target_name = fields[5]
                        target_len = int(fields[6])
                        target_start = int(fields[7])
                        target_end = int(fields[8])
                        matches = int(fields[9])
                        alignment_len = int(fields[10])
                        mapq = int(fields[11])
                        
                        # Calculate metrics
                        query_coverage = (query_end - query_start) / query_len
                        target_coverage = (target_end - target_start) / target_len
                        identity = matches / alignment_len if alignment_len > 0 else 0
                        
                        # Extract pair index from sequence name
                        try:
                            pair_idx = int(query_name.split('_')[1])
                            if pair_idx < len(sequence_pairs):
                                pair = sequence_pairs[pair_idx]
                                
                                # Only keep high-quality alignments
                                if identity >= min_identity and query_coverage >= 0.5 and target_coverage >= 0.5:
                                    result = {
                                        'pair_index': pair_idx,
                                        'busco_id': pair['first_gene']['busco_id'],
                                        'identity': identity,
                                        'query_coverage': query_coverage,
                                        'target_coverage': target_coverage,
                                        'min_coverage': min(query_coverage, target_coverage),
                                        'alignment_length': alignment_len,
                                        'matches': matches,
                                        'mapq': mapq,
                                        'strand': strand,
                                        'method': 'minimap2'
                                    }
                                    results.append(result)
                        except (IndexError, ValueError):
                            continue
    except Exception as e:
        logger.error(f"Error parsing PAF file: {e}")
    
    return results


def run_biopython_alignment_batch(sequence_pairs_batch, config):
    """Run Biopython alignment on a batch of sequence pairs"""
    results = []
    
    # Handle both dict and simple config
    if isinstance(config, dict):
        match_score = config.get('biopython_match_score', 2)
        mismatch_score = config.get('biopython_mismatch_score', -1)
        gap_open_score = config.get('biopython_gap_open_score', -2)
        gap_extend_score = config.get('biopython_gap_extend_score', -0.5)
        mode = config.get('biopython_mode', 'local')
        fallback_enabled = config.get('fallback_to_simple_similarity', True)
    else:
        # Default values if config is not available
        match_score = 2
        mismatch_score = -1
        gap_open_score = -2
        gap_extend_score = -0.5
        mode = 'local'
        fallback_enabled = True
    
    # Setup aligner
    try:
        aligner = PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        aligner.open_gap_score = gap_open_score
        aligner.extend_gap_score = gap_extend_score
        aligner.mode = mode
        use_aligner = True
    except Exception as e:
        logger.warning(f"Failed to setup Biopython aligner: {e}")
        use_aligner = False
    
    for pair_idx, pair in enumerate(sequence_pairs_batch):
        try:
            seq1 = pair['first_gene']['gene_sequence']
            seq2 = pair['second_gene']['gene_sequence']
            busco_id = pair['first_gene']['busco_id']
            
            if use_aligner:
                # Try Biopython alignment
                try:
                    alignments = aligner.align(seq1, seq2)
                    if alignments:
                        best_alignment = alignments[0]
                        
                        # Calculate detailed statistics
                        try:
                            alignment_str = str(best_alignment)
                            alignment_lines = alignment_str.split('\n')
                            if len(alignment_lines) >= 3:
                                seq1_aligned = alignment_lines[0]
                                seq2_aligned = alignment_lines[2]
                                
                                matches = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) 
                                            if a == b and a != '-')
                                alignment_length = len(seq1_aligned)
                                
                                # Coverage calculations
                                query_coverage = (alignment_length - seq1_aligned.count('-')) / len(seq1)
                                target_coverage = (alignment_length - seq2_aligned.count('-')) / len(seq2)
                                identity = matches / alignment_length if alignment_length > 0 else 0
                            else:
                                # Fallback calculation
                                matches = int(best_alignment.score / 2)  # Rough estimate
                                alignment_length = max(len(seq1), len(seq2))
                                query_coverage = 0.8  # Conservative estimate
                                target_coverage = 0.8
                                identity = matches / alignment_length if alignment_length > 0 else 0
                        except Exception as e:
                            # Simple fallback
                            matches = int(best_alignment.score / 2)
                            alignment_length = max(len(seq1), len(seq2))
                            query_coverage = 0.7
                            target_coverage = 0.7
                            identity = matches / alignment_length if alignment_length > 0 else 0
                        
                        result = {
                            'pair_index': pair_idx,
                            'busco_id': busco_id,
                            'identity': identity,
                            'query_coverage': query_coverage,
                            'target_coverage': target_coverage,
                            'min_coverage': min(query_coverage, target_coverage),
                            'alignment_length': alignment_length,
                            'matches': matches,
                            'score': best_alignment.score,
                            'method': 'biopython'
                        }
                        results.append(result)
                        continue
                except Exception as e:
                    # Biopython failed, fall back to difflib
                    pass
            
            # Fallback to difflib
            if fallback_enabled:
                similarity = SequenceMatcher(None, seq1, seq2).ratio()
                len_ratio = min(len(seq1), len(seq2)) / max(len(seq1), len(seq2))
                
                result = {
                    'pair_index': pair_idx,
                    'busco_id': busco_id,
                    'identity': similarity,
                    'query_coverage': 1.0,
                    'target_coverage': 1.0,
                    'min_coverage': 1.0,
                    'alignment_length': min(len(seq1), len(seq2)),
                    'matches': int(similarity * min(len(seq1), len(seq2))),
                    'score': similarity * len_ratio,
                    'method': 'difflib_fallback'
                }
                results.append(result)
                
        except Exception as e:
            # Log error but continue with next pair
            logger.warning(f"Failed to process sequence pair {pair_idx}: {e}")
            continue
    
    return results


def run_simple_biopython_alignment(sequence_pairs, config):
    """Simple single-threaded Biopython alignment as fallback"""
    logger.info(f"  Using simple sequential alignment for {len(sequence_pairs)} sequences...")
    
    results = []
    batch_size = 100
    total_batches = (len(sequence_pairs) + batch_size - 1) // batch_size
    
    for i in range(0, len(sequence_pairs), batch_size):
        batch = sequence_pairs[i:i + batch_size]
        batch_results = run_biopython_alignment_batch(batch, config)
        results.extend(batch_results)
        
        # Progress reporting
        current_batch = (i // batch_size) + 1
        processed = min(i + batch_size, len(sequence_pairs))
        logger.info(f"    Progress: {processed}/{len(sequence_pairs)} sequences ({processed/len(sequence_pairs)*100:.1f}%) - Batch {current_batch}/{total_batches}")
    
    return results


def run_parallel_biopython_alignment(sequence_pairs, config):
    """Run Biopython alignment with multiprocessing"""
    if not sequence_pairs:
        return []
    
    batch_size = config.get('biopython_batch_size', 100)
    max_workers = min(multiprocessing.cpu_count(), config.get('minimap2_threads', 4))
    
    # Split into batches
    batches = [sequence_pairs[i:i + batch_size] for i in range(0, len(sequence_pairs), batch_size)]
    
    all_results = []
    
    # Extract only serializable config parameters
    simple_config = {
        'biopython_match_score': config.get('biopython_match_score', 2),
        'biopython_mismatch_score': config.get('biopython_mismatch_score', -1),
        'biopython_gap_open_score': config.get('biopython_gap_open_score', -2),
        'biopython_gap_extend_score': config.get('biopython_gap_extend_score', -0.5),
        'biopython_mode': config.get('biopython_mode', 'local'),
        'fallback_to_simple_similarity': config.get('fallback_to_simple_similarity', True),
        'timeout_per_alignment': config.get('timeout_per_alignment', 30),
        'detailed_alignment_logging': config.get('detailed_alignment_logging', False)
    }
    
    if config.get('enable_parallel_alignment', True) and len(batches) > 1 and len(sequence_pairs) > 200:
        # Parallel processing for large datasets only
        logger.info(f"    Using parallel processing with {max_workers} workers")
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = []
                for batch in batches:
                    future = executor.submit(run_biopython_alignment_batch, batch, simple_config)
                    futures.append(future)
                
                for i, future in enumerate(futures):
                    try:
                        batch_results = future.result(timeout=simple_config['timeout_per_alignment'])
                        all_results.extend(batch_results)
                        
                        if simple_config['detailed_alignment_logging']:
                            logger.info(f"    Completed batch {i+1}/{len(batches)}")
                        elif i % 5 == 0:  # Progress every 5 batches
                            logger.info(f"    Progress: {i+1}/{len(batches)} batches completed")
                    except Exception as e:
                        logger.warning(f"Batch {i+1} failed: {e}")
                        # Try sequential processing for failed batch
                        try:
                            batch_results = run_biopython_alignment_batch(batches[i], simple_config)
                            all_results.extend(batch_results)
                            logger.info(f"    Batch {i+1} recovered with sequential processing")
                        except Exception as e2:
                            logger.error(f"    Batch {i+1} failed completely: {e2}")
        except Exception as e:
            logger.warning(f"Parallel processing failed: {e}. Falling back to sequential processing")
            # Fall back to sequential processing
            for i, batch in enumerate(batches):
                batch_results = run_biopython_alignment_batch(batch, simple_config)
                all_results.extend(batch_results)
                
                if i % 10 == 0:
                    logger.info(f"    Processed {i * batch_size}/{len(sequence_pairs)} sequence pairs")
    else:
        # Sequential processing
        logger.info(f"    Using sequential processing for {len(sequence_pairs)} sequence pairs")
        for i, batch in enumerate(batches):
            batch_results = run_biopython_alignment_batch(batch, simple_config)
            all_results.extend(batch_results)
            
            # Progress reporting
            if i % 10 == 0 or i == len(batches) - 1:
                processed = min((i + 1) * batch_size, len(sequence_pairs))
                logger.info(f"    Processed {processed}/{len(sequence_pairs)} sequence pairs ({processed/len(sequence_pairs)*100:.1f}%)")
    
    return all_results


def normalize_alignment_scores(alignment_results, config):
    """Normalize alignment scores from different methods to a unified scale"""
    if not config.get('normalize_scores', True):
        return alignment_results
    
    weights = config.get('confidence_weighting', {
        'identity': 0.4,
        'coverage': 0.3,
        'length_ratio': 0.2,
        'score_quality': 0.1
    })
    
    normalized_results = []
    
    for result in alignment_results:
        # Calculate unified confidence score
        identity_score = result['identity']
        coverage_score = result['min_coverage']
        
        # Length ratio - calculated from original pair data if available
        length_ratio = 1.0  # Default if not available
        
        # Score quality - method-dependent normalization
        if result['method'] == 'minimap2':
            score_quality = min(1.0, result.get('mapq', 0) / 60.0)  # MAPQ is typically 0-60
        elif result['method'] == 'biopython':
            score_quality = min(1.0, max(0.0, (result.get('score', 0) + 100) / 200.0))  # Rough normalization
        else:
            score_quality = result.get('identity', 0.5)  # Fallback methods
        
        # Calculate weighted confidence
        confidence = (
            weights['identity'] * identity_score +
            weights['coverage'] * coverage_score +
            weights['length_ratio'] * length_ratio +
            weights['score_quality'] * score_quality
        )
        
        # Add normalized scores to result
        result['confidence'] = confidence
        result['similarity'] = identity_score * coverage_score  # Combined similarity metric
        
        normalized_results.append(result)
    
    return normalized_results


def apply_reciprocal_best_hit_filtering(alignment_results, config):
    """Apply reciprocal best hit filtering to remove many-to-many matches"""
    if not config.get('use_reciprocal_best_hits', True):
        return alignment_results
    
    # Group results by BUSCO ID
    busco_results = {}
    for result in alignment_results:
        busco_id = result['busco_id']
        if busco_id not in busco_results:
            busco_results[busco_id] = []
        busco_results[busco_id].append(result)
    
    # Apply RBH filtering per BUSCO
    filtered_results = []
    identity_gap_threshold = config.get('identity_gap_threshold', 0.02)
    
    for busco_id, results in busco_results.items():
        if len(results) <= 1:
            # Only one result, keep it
            filtered_results.extend(results)
        else:
            # Multiple results - apply tie-breaking
            # Sort by confidence, then identity, then coverage
            results.sort(key=lambda x: (x['confidence'], x['identity'], x['min_coverage']), reverse=True)
            
            best_result = results[0]
            
            # Check for ties within identity gap threshold
            ties = [r for r in results if abs(r['identity'] - best_result['identity']) <= identity_gap_threshold]
            
            if len(ties) == 1:
                filtered_results.append(best_result)
            else:
                # Multiple ties - use additional tie-breakers
                ties.sort(key=lambda x: (x['matches'], -x.get('mapq', 0), x['alignment_length']), reverse=True)
                filtered_results.append(ties[0])
    
    return filtered_results


def convert_alignment_results_to_ortholog_pairs(alignment_results, sequence_pairs, config):
    """Convert alignment results back to ortholog pairs format"""
    ortholog_pairs = []
    
    # Create mapping from pair index to sequence pair
    pair_map = {i: pair for i, pair in enumerate(sequence_pairs)}
    
    for result in alignment_results:
        pair_idx = result.get('pair_index')
        if pair_idx is not None and pair_idx in pair_map:
            pair = pair_map[pair_idx]
            first_gene = pair['first_gene']
            second_gene = pair['second_gene']
            
            ortholog_pair = {
                'busco_id': result['busco_id'],
                'first_chr': standardize_sequence_id(first_gene['sequence_id']),
                'first_start': first_gene['gene_start'],
                'first_end': first_gene['gene_end'],
                'first_strand': first_gene['strand'],
                'second_chr': standardize_sequence_id(second_gene['sequence_id']),
                'second_start': second_gene['gene_start'],
                'second_end': second_gene['gene_end'],
                'second_strand': second_gene['strand'],
                'similarity': result['similarity'],
                'confidence': result['confidence'],
                'first_length': first_gene['gene_length'],
                'second_length': second_gene['gene_length'],
                'identity': result['identity'],
                'coverage': result['min_coverage'],
                'alignment_method': result['method'],
                'alignment_length': result['alignment_length'],
                'matches': result['matches'],
                'first_paralog_rank': first_gene.get('paralog_rank', 1),
                'second_paralog_rank': second_gene.get('paralog_rank', 1),
                'first_paralog_count': first_gene.get('paralog_count', 1),
                'second_paralog_count': second_gene.get('paralog_count', 1),
                'length_ratio': min(first_gene['gene_length'], second_gene['gene_length']) / max(first_gene['gene_length'], second_gene['gene_length']),
                'mapping_type': 'ortholog'
            }
            
            # Add method-specific metadata
            if result['method'] == 'minimap2':
                ortholog_pair['mapq'] = result.get('mapq', 0)
            elif result['method'] == 'biopython':
                ortholog_pair['alignment_score'] = result.get('score', 0)
            
            # Add validation information if available
            if 'validation_method' in result:
                ortholog_pair['validation_method'] = result['validation_method']
                if 'alternative_confidence' in result:
                    ortholog_pair['alternative_confidence'] = result['alternative_confidence']
            
            ortholog_pairs.append(ortholog_pair)
    
    return ortholog_pairs


def run_hybrid_alignment_analysis(first_busco_df, second_busco_df, config):
    """Main orchestrator for hybrid alignment analysis"""
    logger.info("Running hybrid alignment analysis...")
    
    # Find common BUSCO genes and create pairs
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
    
    # Check if minimap2 is available
    minimap2_available = check_minimap2_available()
    if not minimap2_available:
        logger.warning("  Minimap2 not available, forcing Biopython-only mode")
        config['alignment_strategy'] = 'biopython'
    
    # Get strategy from config, default to hybrid if minimap2 is available
    strategy = config.get('alignment_strategy', 'hybrid' if minimap2_available else 'biopython')
    logger.info(f"  Using alignment strategy: {strategy}")
    
    # Partition sequences by length for optimal alignment method
    partitions = partition_sequences_by_length(sequence_pairs, config)
    
    # Run alignments using appropriate methods
    all_results = []
    
    if strategy == 'hybrid' and minimap2_available:
        # HYBRID MODE: Use different methods for different sequence lengths
        
        # Short sequences -> Biopython
        if partitions['short_pairs']:
            logger.info(f"  Running Biopython alignment on {len(partitions['short_pairs'])} short sequences...")
            try:
                short_results = run_parallel_biopython_alignment(partitions['short_pairs'], config)
            except Exception as e:
                logger.warning(f"Parallel alignment failed for short sequences: {e}")
                short_results = run_simple_biopython_alignment(partitions['short_pairs'], config)
            all_results.extend(short_results)
        
        # Long sequences -> Minimap2
        if partitions['long_pairs']:
            logger.info(f"  Running Minimap2 alignment on {len(partitions['long_pairs'])} long sequences...")
            long_results = run_minimap2_alignment(partitions['long_pairs'], config)
            if not long_results:
                logger.warning("  Minimap2 failed, falling back to Biopython for long sequences")
                try:
                    long_results = run_parallel_biopython_alignment(partitions['long_pairs'], config)
                except Exception as e:
                    logger.warning(f"Parallel alignment failed for long sequences: {e}")
                    long_results = run_simple_biopython_alignment(partitions['long_pairs'], config)
            all_results.extend(long_results)
        
        # Buffer zone -> Both methods or fallback
        if partitions['buffer_pairs']:
            buffer_method = config.get('buffer_zone_method', 'dual')
            logger.info(f"  Processing {len(partitions['buffer_pairs'])} buffer zone sequences with {buffer_method} method...")
            
            if buffer_method == 'dual' and config.get('cross_validate_buffer_zone', True):
                # Run both methods and compare
                try:
                    bio_results = run_parallel_biopython_alignment(partitions['buffer_pairs'], config)
                except Exception as e:
                    logger.warning(f"Parallel alignment failed for buffer zone: {e}")
                    bio_results = run_simple_biopython_alignment(partitions['buffer_pairs'], config)
                mm2_results = run_minimap2_alignment(partitions['buffer_pairs'], config)
                if mm2_results:
                    buffer_results = select_best_buffer_results(bio_results, mm2_results, config)
                else:
                    logger.warning("  Minimap2 failed for buffer zone, using Biopython results")
                    buffer_results = bio_results
            elif buffer_method == 'minimap2':
                buffer_results = run_minimap2_alignment(partitions['buffer_pairs'], config)
                if not buffer_results:
                    logger.warning("  Minimap2 failed for buffer zone, falling back to Biopython")
                    try:
                        buffer_results = run_parallel_biopython_alignment(partitions['buffer_pairs'], config)
                    except Exception as e:
                        logger.warning(f"Parallel alignment failed for buffer zone: {e}")
                        buffer_results = run_simple_biopython_alignment(partitions['buffer_pairs'], config)
            else:  # Default to Biopython
                try:
                    buffer_results = run_parallel_biopython_alignment(partitions['buffer_pairs'], config)
                except Exception as e:
                    logger.warning(f"Parallel alignment failed for buffer zone: {e}")
                    buffer_results = run_simple_biopython_alignment(partitions['buffer_pairs'], config)
            
            all_results.extend(buffer_results)
        
        # Mixed length pairs -> Biopython (safer for heterogeneous lengths)
        if partitions['mixed_pairs']:
            logger.info(f"  Running Biopython alignment on {len(partitions['mixed_pairs'])} mixed-length sequences...")
            try:
                mixed_results = run_parallel_biopython_alignment(partitions['mixed_pairs'], config)
            except Exception as e:
                logger.warning(f"Parallel alignment failed for mixed sequences: {e}")
                mixed_results = run_simple_biopython_alignment(partitions['mixed_pairs'], config)
            all_results.extend(mixed_results)
    
    elif strategy == 'minimap2' and minimap2_available:
        # MINIMAP2 ONLY MODE
        logger.info(f"  Running Minimap2 alignment on all {len(sequence_pairs)} sequences...")
        all_results = run_minimap2_alignment(sequence_pairs, config)
        if not all_results:
            logger.warning("  Minimap2 failed completely, falling back to Biopython")
            try:
                all_results = run_parallel_biopython_alignment(sequence_pairs, config)
            except Exception as e:
                logger.warning(f"Parallel alignment failed: {e}")
                logger.info("  Falling back to simple sequential alignment...")
                all_results = run_simple_biopython_alignment(sequence_pairs, config)
    
    else:  # BIOPYTHON ONLY MODE
        logger.info(f"  Running Biopython alignment on all {len(sequence_pairs)} sequences...")
        try:
            all_results = run_parallel_biopython_alignment(sequence_pairs, config)
        except Exception as e:
            logger.warning(f"Parallel alignment failed: {e}")
            logger.info("  Falling back to simple sequential alignment...")
            all_results = run_simple_biopython_alignment(sequence_pairs, config)
    
    # Normalize scores and calculate confidence
    logger.info("  Normalizing alignment scores and calculating confidence...")
    normalized_results = normalize_alignment_scores(all_results, config)
    
    # Apply reciprocal best hit filtering
    if config.get('use_reciprocal_best_hits', True):
        logger.info("  Applying reciprocal best hit filtering...")
        filtered_results = apply_reciprocal_best_hit_filtering(normalized_results, config)
    else:
        filtered_results = normalized_results
    
    # Convert results to ortholog pairs format
    ortholog_pairs = convert_alignment_results_to_ortholog_pairs(
        filtered_results, sequence_pairs, config
    )
    
    # Report statistics
    logger.info(f"  Alignment results:")
    logger.info(f"    Total alignments attempted: {len(sequence_pairs)}")
    logger.info(f"    Successful alignments: {len(all_results)}")
    logger.info(f"    After RBH filtering: {len(filtered_results)}")
    logger.info(f"    Final ortholog pairs: {len(ortholog_pairs)}")
    
    if filtered_results:
        method_counts = {}
        for result in filtered_results:
            method = result['method']
            method_counts[method] = method_counts.get(method, 0) + 1
        
        logger.info(f"    Methods used: {method_counts}")
        
        avg_identity = np.mean([r['identity'] for r in filtered_results])
        avg_confidence = np.mean([r['confidence'] for r in filtered_results])
        logger.info(f"    Average identity: {avg_identity:.3f}")
        logger.info(f"    Average confidence: {avg_confidence:.3f}")
    
    return pd.DataFrame(ortholog_pairs), pd.DataFrame()  # Return empty paralog_df for compatibility


def check_minimap2_available():
    """Check if minimap2 is available in the system PATH"""
    try:
        result = subprocess.run(['minimap2', '--version'], 
                              capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            logger.info(f"  Minimap2 detected: {result.stdout.strip()}")
            return True
        else:
            logger.warning("  Minimap2 found but returned non-zero exit code")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        logger.warning("  Minimap2 not found in PATH")
        return False


def select_best_buffer_results(bio_results, mm2_results, config):
    """Select best results when both Biopython and Minimap2 are run on buffer zone"""
    # Create mapping of results by BUSCO ID
    bio_map = {r['busco_id']: r for r in bio_results}
    mm2_map = {r['busco_id']: r for r in mm2_results}
    
    best_results = []
    
    all_buscos = set(bio_map.keys()) | set(mm2_map.keys())
    
    for busco_id in all_buscos:
        bio_result = bio_map.get(busco_id)
        mm2_result = mm2_map.get(busco_id)
        
        if bio_result and mm2_result:
            # Both methods succeeded - choose best based on confidence
            bio_confidence = bio_result('confidence', 0.0)
            mm2_confidence = mm2_result.get('confidence', 0.0)
            if bio_confidence > mm2_confidence:
                bio_result['validation_method'] = 'dual_validated'
                bio_result['alternative_confidence'] = mm2_result['confidence']
                best_results.append(bio_result)
            else:
                mm2_result['validation_method'] = 'dual_validated'
                mm2_result['alternative_confidence'] = bio_result['confidence']
                best_results.append(mm2_result)
        elif bio_result:
            bio_result['validation_method'] = 'biopython_only'
            best_results.append(bio_result)
        elif mm2_result:
            mm2_result['validation_method'] = 'minimap2_only'
            best_results.append(mm2_result)
    
    return best_results


def setup_hybrid_sequence_aligner(config):
    """Setup hybrid aligner system based on configuration"""
    strategy = config.get('alignment_strategy', 'hybrid')
    
    if strategy == 'biopython':
        # Traditional Biopython aligner
        aligner = PairwiseAligner()
        aligner.match_score = config.get('biopython_match_score', 2)
        aligner.mismatch_score = config.get('biopython_mismatch_score', -1)
        aligner.open_gap_score = config.get('biopython_gap_open_score', -2)
        aligner.extend_gap_score = config.get('biopython_gap_extend_score', -0.5)
        aligner.mode = config.get('biopython_mode', 'local')
        return aligner
    
    elif strategy in ['minimap2', 'hybrid']:
        # Check if minimap2 is available
        try:
            subprocess.run(['minimap2', '--version'], capture_output=True, check=True)
            logger.info(f"  Using {strategy} alignment strategy with minimap2")
            return None  # No Biopython aligner needed for pure minimap2
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("  Minimap2 not found, falling back to Biopython")
            config['alignment_strategy'] = 'biopython'
            return setup_hybrid_sequence_aligner(config)
    
    return None