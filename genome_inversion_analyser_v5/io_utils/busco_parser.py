# =============================================================================
# BUSCO Parser (busco_parser.py)
# =============================================================================

"""
Enhanced BUSCO table parser with comprehensive validation and error handling.
Handles complex BUSCO formats and provides detailed parsing statistics.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple, NamedTuple
from dataclasses import dataclass
from collections import defaultdict, Counter

from ..logger import get_logger
from ..utils import standardize_sequence_id

logger = get_logger()

@dataclass
class BuscoRecord:
    """Represents a single BUSCO gene record with validation."""
    busco_id: str
    status: str
    sequence: str
    gene_start: int
    gene_end: int
    strand: str
    score: Optional[float] = None
    length: Optional[int] = None
    line_number: Optional[int] = None
    original_start: Optional[int] = None
    original_end: Optional[int] = None
    is_paralog: bool = False
    paralog_rank: int = 1
    paralog_count: int = 1
    paralog_cluster: str = ""
    
    def __post_init__(self):
        """Validate record after initialization."""
        # Ensure coordinates are properly ordered
        if self.gene_start > self.gene_end:
            self.gene_start, self.gene_end = self.gene_end, self.gene_start
        
        # Calculate length if not provided
        if self.length is None:
            self.length = self.gene_end - self.gene_start
        
        # Set paralog cluster
        if not self.paralog_cluster:
            self.paralog_cluster = f"{self.busco_id}_singleton"
    
    def to_dict(self) -> Dict:
        """Convert record to dictionary."""
        return {
            'busco_id': self.busco_id,
            'status': self.status,
            'sequence': standardize_sequence_id(self.sequence),  # Standardize here
            'gene_start': self.gene_start,
            'gene_end': self.gene_end,
            'strand': self.strand,
            'score': self.score,
            'length': self.length,
            'line_number': self.line_number,
            'original_start': self.original_start,
            'original_end': self.original_end,
            'is_paralog': self.is_paralog,
            'paralog_rank': self.paralog_rank,
            'paralog_count': self.paralog_count,
            'paralog_cluster': self.paralog_cluster
        }

class BuscoParser:
    """
    Enhanced BUSCO table parser with comprehensive validation and error handling.
    Correctly handles negative strand genes and provides detailed parsing statistics.
    """
    
    def __init__(self, busco_path: str):
        """
        Initialize BUSCO parser.
        
        Args:
            busco_path: Path to BUSCO full_table.tsv file
        """
        self.busco_path = Path(busco_path)
        self.records = []
        self.parsing_errors = []
        self.duplicate_entries = []
        self.statistics = {}
        
        if not self.busco_path.exists():
            raise FileNotFoundError(f"BUSCO file not found: {self.busco_path}")
        
        logger.info(f"Initializing BUSCO parser for: {self.busco_path}")
    
    def parse_busco_table(self, config) -> pd.DataFrame:
        """
        Parse BUSCO table with enhanced validation.
        
        Args:
            config: Configuration object with parsing settings
            
        Returns:
            DataFrame with parsed BUSCO records
        """
        logger.info(f"Parsing BUSCO table: {self.busco_path}")
        
        # Read file and filter comments
        with open(self.busco_path, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]
        
        # Parse each line
        busco_records = []
        seen_entries = set()
        
        for line_num, line in enumerate(lines, 1):
            if line.strip():
                record, error = self._parse_busco_line(line, line_num, config)
                
                if record:
                    # Check for duplicates if enabled
                    if config.get('enable_duplicate_handling', False):
                        entry_key = (record.busco_id, record.sequence, 
                                   record.gene_start, record.gene_end)
                        
                        if entry_key in seen_entries:
                            self.duplicate_entries.append(line_num)
                            continue
                        seen_entries.add(entry_key)
                    
                    busco_records.append(record)
                
                if error:
                    self.parsing_errors.append(error)
        
        self.records = busco_records
        
        # Convert to DataFrame
        busco_df = pd.DataFrame([record.to_dict() for record in busco_records])
        
        # Detect and annotate paralogs
        if config.get('enable_paralog_detection', False):
            busco_df = self._detect_and_annotate_paralogs(busco_df, config)
        
        # Calculate and report statistics
        self._calculate_parsing_statistics(lines, busco_df, config)
        
        return busco_df
    
    def _parse_busco_line(self, line: str, line_num: int, config) -> Tuple[Optional[BuscoRecord], Optional[str]]:
        """
        Parse a single BUSCO table line.
        
        Args:
            line: Line to parse
            line_num: Line number for error reporting
            config: Configuration object
            
        Returns:
            Tuple of (BuscoRecord or None, error_message or None)
        """
        try:
            parts = line.strip().split('\t')
            
            if len(parts) < 6:
                return None, f"Line {line_num}: Too few columns ({len(parts)})"
            
            # Extract basic fields
            busco_id = parts[0]
            status = parts[1]
            sequence = standardize_sequence_id(parts[2])  # Standardize sequence ID immediately
            start_str = parts[3]
            end_str = parts[4]
            strand = parts[5] if len(parts) > 5 else '+'
            
            # Validate sequence ID
            if sequence is None:
                return None, f"Line {line_num}: Invalid sequence ID"
            
            # Validate coordinates
            if start_str == 'N/A' or end_str == 'N/A':
                return None, f"Line {line_num}: N/A coordinates"
            
            try:
                coord1 = int(start_str)
                coord2 = int(end_str)
            except ValueError as e:
                return None, f"Line {line_num}: Invalid coordinates - {e}"
            
            # Determine proper start/end and strand
            gene_start = min(coord1, coord2)
            gene_end = max(coord1, coord2)
            
            # Correct strand based on coordinate order
            if coord1 > coord2 and strand == '+':
                strand = '-'
            elif coord1 < coord2 and strand == '-':
                strand = '+'
            
            # Extract optional fields
            score = None
            if len(parts) > 6 and parts[6] != 'N/A':
                try:
                    score = float(parts[6])
                except ValueError:
                    pass
            
            length = None
            if len(parts) > 7 and parts[7] != 'N/A':
                try:
                    length = int(parts[7])
                except ValueError:
                    length = gene_end - gene_start
            else:
                length = gene_end - gene_start
            
            # Create record
            record = BuscoRecord(
                busco_id=busco_id,
                status=status,
                sequence=sequence,
                gene_start=gene_start,
                gene_end=gene_end,
                strand=strand,
                score=score,
                length=length,
                line_number=line_num,
                original_start=coord1,
                original_end=coord2
            )
            
            return record, None
            
        except Exception as e:
            return None, f"Line {line_num}: Parsing error - {e}"
    
    def _detect_and_annotate_paralogs(self, busco_df: pd.DataFrame, config) -> pd.DataFrame:
        """
        Enhanced paralog detection with detailed analysis.
        
        Args:
            busco_df: DataFrame with BUSCO records
            config: Configuration object
            
        Returns:
            DataFrame with paralog annotations
        """
        logger.info("  Detecting and analyzing paralogs...")
        
        # Count occurrences of each BUSCO ID
        busco_counts = busco_df['busco_id'].value_counts()
        paralogs = busco_counts[busco_counts > 1]
        
        # Add paralog annotations
        busco_df['is_paralog'] = busco_df['busco_id'].isin(paralogs.index)
        busco_df['paralog_count'] = busco_df['busco_id'].map(busco_counts)
        
        # Initialize rank and cluster columns
        busco_df['paralog_rank'] = 1
        busco_df['paralog_cluster'] = busco_df['busco_id'] + '_singleton'
        
        # Process each BUSCO ID with paralogs
        for busco_id in paralogs.index:
            paralog_group = busco_df[busco_df['busco_id'] == busco_id].copy()
            
            # Multi-criteria ranking
            ranks = self._rank_paralogs(paralog_group)
            
            # Update main dataframe
            busco_df.loc[paralog_group.index, 'paralog_rank'] = ranks
            busco_df.loc[paralog_group.index, 'paralog_cluster'] = f"{busco_id}_cluster"
        
        logger.info(f"    Found {len(paralogs)} BUSCO IDs with paralogs")
        logger.info(f"    Total paralogous genes: {len(busco_df[busco_df['is_paralog']])}")
        logger.info(f"    Max paralogs for single BUSCO: {busco_counts.max()}")
        
        return busco_df
    
    def _rank_paralogs(self, paralog_group: pd.DataFrame) -> pd.Series:
        """
        Rank paralogs using multiple criteria.
        
        Args:
            paralog_group: DataFrame with paralogs for a single BUSCO ID
            
        Returns:
            Series with paralog ranks
        """
        ranking_criteria = []
        weights = []
        
        # Primary: score (if available)
        if 'score' in paralog_group.columns and paralog_group['score'].notna().any():
            ranking_criteria.append(paralog_group['score'].fillna(0))
            weights.append(0.5)
        
        # Secondary: length (if available)
        if 'length' in paralog_group.columns and paralog_group['length'].notna().any():
            ranking_criteria.append(paralog_group['length'].fillna(0))
            weights.append(0.3)
        
        # Tertiary: genomic position (earlier positions get higher rank)
        if paralog_group['gene_start'].notna().any():
            ranking_criteria.append(-paralog_group['gene_start'].fillna(float('inf')))
            weights.append(0.2)
        
        if ranking_criteria:
            # Normalize weights
            weights = np.array(weights) / sum(weights)
            # Combine criteria with weights
            combined_score = sum(w * criteria for w, criteria in zip(weights, ranking_criteria))
            ranks = combined_score.rank(ascending=False, method='first')
        else:
            # If no ranking criteria available, rank by dataframe order
            ranks = pd.Series(range(1, len(paralog_group) + 1), index=paralog_group.index)
        
        return ranks
    
    def _calculate_parsing_statistics(self, lines: List[str], busco_df: pd.DataFrame, config):
        """Calculate comprehensive parsing statistics."""
        total_lines = len(lines)
        successful = len(busco_df)
        failed = len(self.parsing_errors)
        
        success_rate = successful / total_lines * 100 if total_lines > 0 else 0
        
        self.statistics = {
            'total_lines': total_lines,
            'successful_parses': successful,
            'failed_parses': failed,
            'success_rate': success_rate,
            'duplicate_entries': len(self.duplicate_entries),
            'unique_buscos': busco_df['busco_id'].nunique() if len(busco_df) > 0 else 0,
            'status_counts': busco_df['status'].value_counts().to_dict() if len(busco_df) > 0 else {},
            'paralogs_detected': len(busco_df[busco_df['is_paralog']]) if 'is_paralog' in busco_df.columns else 0
        }
        
        # Report results
        logger.info(f"  Parsing results:")
        logger.info(f"    Total data lines: {total_lines}")
        logger.info(f"    Successfully parsed: {successful} ({success_rate:.1f}%)")
        logger.info(f"    Parsing errors: {failed}")
        
        if failed > 0 and config.get('enable_debug_output', False):
            logger.info(f"    First few errors:")
            for error in self.parsing_errors[:3]:
                logger.info(f"      {error}")
        
        if self.duplicate_entries and config.get('enable_debug_output', False):
            logger.warning(f"  {len(self.duplicate_entries)} duplicate entries removed")
    
    def filter_busco_genes(self, busco_df: pd.DataFrame, config, quality_info=None) -> pd.DataFrame:
        """
        Enhanced BUSCO filtering with adaptive parameters.
        
        Args:
            busco_df: DataFrame with BUSCO records
            config: Configuration object
            quality_info: Assembly quality information for adaptive filtering
            
        Returns:
            Filtered DataFrame
        """
        logger.info("Enhanced BUSCO gene filtering...")
        
        initial_count = len(busco_df)
        filtering_stats = {'initial': initial_count}
        
        # Get adaptive parameters
        if quality_info and config.get('enable_adaptive_thresholds', False):
            filter_params = quality_info['adjustments']
        else:
            filter_params = {
                'min_busco_length': config.get('base_min_busco_length', 150),
                'similarity_threshold': config.get('base_similarity_threshold', 0.5)
            }
        
        # Filter by status
        status_filter = config.get('busco_status_filter', ['Complete'])
        filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
        filtering_stats['status_filter'] = len(filtered_df)
        
        # Filter by length
        if 'length' in filtered_df.columns and filter_params.get('min_busco_length'):
            length_before = len(filtered_df)
            min_length = filter_params['min_busco_length']
            filtered_df = filtered_df[filtered_df['length'] >= min_length]
            filtering_stats['length_filter'] = len(filtered_df)
            logger.info(f"  Length filter (>={min_length}bp): {length_before} -> {len(filtered_df)}")
        
        # Remove genes with missing coordinates
        coord_before = len(filtered_df)
        filtered_df = filtered_df.dropna(subset=['gene_start', 'gene_end'])
        filtering_stats['coordinate_filter'] = len(filtered_df)
        
        # Paralog handling
        if config.get('enable_paralog_detection', False):
            if config.get('enable_paralog_ortholog_mapping', False):
                # Keep all paralogs for downstream analysis
                pass
            else:
                # Keep only best paralog per BUSCO ID
                paralog_before = len(filtered_df)
                filtered_df = filtered_df.sort_values('paralog_rank').groupby('busco_id').first().reset_index()
                filtering_stats['paralog_filter'] = len(filtered_df)
                logger.info(f"  Paralog filter (best per BUSCO): {paralog_before} -> {len(filtered_df)}")
        
        final_count = len(filtered_df)
        exclusion_rate = (initial_count - final_count) / initial_count if initial_count > 0 else 0
        
        # Warning about excessive exclusions
        if config.get('enable_exclusion_warnings', False):
            if exclusion_rate > 0.7:
                logger.warning(f"Very high exclusion rate: {exclusion_rate:.1%} of BUSCOs excluded")
                logger.warning("Consider relaxing filtering parameters or checking assembly quality")
            elif exclusion_rate > 0.5:
                logger.warning(f"High exclusion rate: {exclusion_rate:.1%} of BUSCOs excluded")
        
        logger.info(f"  Filtering summary: {initial_count} -> {final_count} ({exclusion_rate:.1%} excluded)")
        
        # Add filtering statistics
        if config.get('enable_debug_output', False):
            filtered_df.attrs['filtering_stats'] = filtering_stats
        
        return filtered_df
    
    def get_parsing_statistics(self) -> Dict:
        """Get comprehensive parsing statistics."""
        return self.statistics.copy()
    
    def get_parsing_errors(self) -> List[str]:
        """Get list of parsing errors."""
        return self.parsing_errors.copy()