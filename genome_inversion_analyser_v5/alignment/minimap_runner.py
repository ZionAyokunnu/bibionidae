
# =============================================================================
# Minimap2 Runner (minimap_runner.py)
# =============================================================================

"""
Minimap2 subprocess runner with comprehensive PAF parsing.
Optimized for long sequences with fast alignment and detailed output parsing.
"""

import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from ..logger import get_logger
from .alignment_result import AlignmentResult

logger = get_logger()

class Minimap2Runner:
    """
    Minimap2 subprocess runner with comprehensive PAF output parsing.
    Handles temporary file management and detailed alignment statistics.
    """
    
    def __init__(self, config):
        """
        Initialize Minimap2 runner.
        
        Args:
            config: Configuration object with Minimap2 parameters
        """
        self.config = config
        self.minimap2_available = self._check_minimap2_availability()
        
        if self.minimap2_available:
            logger.info("Minimap2 runner initialized")
            logger.info(f"  Preset: {config.get('minimap2_preset', '--sr')}")
            logger.info(f"  K-mer size: {config.get('minimap2_kmer_size', 13)}")
            logger.info(f"  Threads: {config.get('minimap2_threads', 4)}")
            logger.info(f"  Min identity: {config.get('minimap2_min_identity', 0.7)}")
        else:
            logger.warning("Minimap2 not available - will fallback to Biopython")
    
    def _check_minimap2_availability(self) -> bool:
        """Check if Minimap2 is available in the system PATH."""
        try:
            result = subprocess.run(
                ['minimap2', '--version'], 
                capture_output=True, 
                timeout=10,
                text=True
            )
            if result.returncode == 0:
                version = result.stdout.strip()
                logger.info(f"  Minimap2 version: {version}")
                return True
        except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
            pass
        
        return False
    
    def align_sequences(self, sequence_pairs: List[Dict]) -> List[AlignmentResult]:
        """
        Align sequence pairs using Minimap2.
        
        Args:
            sequence_pairs: List of sequence pair dictionaries
            
        Returns:
            List of AlignmentResult objects
        """
        if not sequence_pairs:
            return []
        
        if not self.minimap2_available:
            logger.error("Minimap2 not available for alignment")
            return []
        
        logger.info(f"Running Minimap2 alignment on {len(sequence_pairs)} sequence pairs...")
        
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                # Create FASTA files
                query_file, target_file, pair_mapping = self._create_minimap2_fasta(
                    sequence_pairs, temp_path
                )
                
                # Run Minimap2
                paf_results = self._run_minimap2(query_file, target_file, temp_path)
                
                # Parse results
                alignment_results = self._parse_paf_results(
                    paf_results, sequence_pairs, pair_mapping
                )
                
                # Cleanup if requested
                if self.config.get('temp_file_cleanup', True):
                    # Temporary directory will be cleaned up automatically
                    pass
                
                return alignment_results
                
        except Exception as e:
            logger.error(f"Minimap2 alignment failed: {e}")
            return []
    
    def _create_minimap2_fasta(self, sequence_pairs: List[Dict], temp_dir: Path) -> Tuple[Path, Path, Dict]:
        """
        Create FASTA files for Minimap2 alignment.
        
        Args:
            sequence_pairs: List of sequence pairs
            temp_dir: Temporary directory path
            
        Returns:
            Tuple of (query_file, target_file, pair_mapping)
        """
        query_file = temp_dir / "queries.fasta"
        target_file = temp_dir / "targets.fasta"
        
        # Create sequence mappings for result parsing
        pair_mapping = {}
        
        with open(query_file, 'w') as qf, open(target_file, 'w') as tf:
            for i, pair in enumerate(sequence_pairs):
                # Create standardized sequence IDs
                query_id = f"query_{i}_{pair['first_gene']['busco_id']}"
                target_id = f"target_{i}_{pair['second_gene']['busco_id']}"
                
                # Write sequences
                qf.write(f">{query_id}\n{pair['first_gene']['gene_sequence']}\n")
                tf.write(f">{target_id}\n{pair['second_gene']['gene_sequence']}\n")
                
                # Store mapping for result parsing
                pair_mapping[query_id] = i
                pair_mapping[target_id] = i
        
        return query_file, target_file, pair_mapping
    
    def _run_minimap2(self, query_file: Path, target_file: Path, temp_dir: Path) -> Path:
        """
        Execute Minimap2 alignment.
        
        Args:
            query_file: Path to query FASTA file
            target_file: Path to target FASTA file
            temp_dir: Temporary directory for output
            
        Returns:
            Path to PAF output file
        """
        output_file = temp_dir / "alignments.paf"
        
        # Build Minimap2 command
        cmd = [
            'minimap2',
            self.config.get('minimap2_preset', '--sr'),
            '-k', str(self.config.get('minimap2_kmer_size', 13)),
            '-t', str(self.config.get('minimap2_threads', 4)),
            '--score-N', str(self.config.get('minimap2_min_score', 100))
        ]
        
        # Add extra flags
        extra_flags = self.config.get('minimap2_extra_flags', '-c --cs')
        if extra_flags:
            cmd.extend(extra_flags.split())
        
        # Add input files
        cmd.extend([str(target_file), str(query_file)])
        
        # Remove empty arguments
        cmd = [c for c in cmd if c.strip()]
        
        if self.config.get('detailed_alignment_logging', False):
            logger.info(f"    Running: {' '.join(cmd[:6])}...")
        
        # Execute Minimap2
        timeout = self.config.get('timeout_per_alignment', 30) * len(open(query_file).readlines()) // 2
        
        with open(output_file, 'w') as out_f:
            process = subprocess.run(
                cmd,
                stdout=out_f,
                stderr=subprocess.PIPE,
                timeout=timeout,
                text=True
            )
        
        if process.returncode != 0:
            error_msg = process.stderr.strip() if process.stderr else "Unknown error"
            raise subprocess.CalledProcessError(process.returncode, cmd, error_msg)
        
        return output_file
    
    def _parse_paf_results(self, paf_file: Path, sequence_pairs: List[Dict], 
                          pair_mapping: Dict) -> List[AlignmentResult]:
        """
        Parse Minimap2 PAF output into AlignmentResult objects.
        
        Args:
            paf_file: Path to PAF output file
            sequence_pairs: Original sequence pairs
            pair_mapping: Mapping from sequence IDs to pair indices
            
        Returns:
            List of AlignmentResult objects
        """
        results = []
        min_identity = self.config.get('minimap2_min_identity', 0.7)
        
        try:
            with open(paf_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line.strip():
                        try:
                            result = self._parse_paf_line(line, sequence_pairs, pair_mapping, min_identity)
                            if result:
                                results.append(result)
                        except Exception as e:
                            if self.config.get('detailed_alignment_logging', False):
                                logger.warning(f"Failed to parse PAF line {line_num}: {e}")
                            continue
        
        except Exception as e:
            logger.error(f"Error reading PAF file: {e}")
        
        logger.info(f"  Parsed {len(results)} high-quality alignments from Minimap2")
        return results
    
    def _parse_paf_line(self, line: str, sequence_pairs: List[Dict], 
                       pair_mapping: Dict, min_identity: float) -> Optional[AlignmentResult]:
        """
        Parse a single PAF line into an AlignmentResult.
        
        Args:
            line: PAF format line
            sequence_pairs: Original sequence pairs
            pair_mapping: Sequence ID to pair index mapping
            min_identity: Minimum identity threshold
            
        Returns:
            AlignmentResult object or None if below threshold
        """
        fields = line.strip().split('\t')
        if len(fields) < 12:
            return None
        
        # Parse PAF fields (format: query_name, query_len, query_start, query_end, strand, 
        # target_name, target_len, target_start, target_end, matches, alignment_len, mapq)
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
        
        # Calculate coverage and identity
        query_coverage = (query_end - query_start) / query_len if query_len > 0 else 0
        target_coverage = (target_end - target_start) / target_len if target_len > 0 else 0
        identity = matches / alignment_len if alignment_len > 0 else 0
        
        # Quality filters
        if identity < min_identity or query_coverage < 0.5 or target_coverage < 0.5:
            return None
        
        # Extract pair index from query name
        try:
            pair_idx = int(query_name.split('_')[1])
            if pair_idx >= len(sequence_pairs):
                return None
                
            busco_id = sequence_pairs[pair_idx]['first_gene']['busco_id']
            
            # Create PAF result dictionary
            paf_data = {
                'identity': identity,
                'query_coverage': query_coverage,
                'target_coverage': target_coverage,
                'alignment_length': alignment_len,
                'matches': matches,
                'mapq': mapq,
                'strand': strand
            }
            
            return AlignmentResult.from_minimap2(pair_idx, busco_id, paf_data)
            
        except (IndexError, ValueError):
            return None