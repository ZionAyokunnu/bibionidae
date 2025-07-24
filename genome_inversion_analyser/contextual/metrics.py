"""
Contextual Metrics Module for Genome Inversion Analyzer
Adds biological context to inversion analysis with functional annotations
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from collections import defaultdict
import tempfile
import subprocess

logger = logging.getLogger(__name__)


class ContextualMetrics:
    """
    Analyzes inversions in biological context with functional annotations
    """
    
    def __init__(self, registry, config: Dict = None):
        self.registry = registry
        self.config = config or {}
        self.functional_data = {}
        self.repeat_data = {}
        self.gc_content = {}
        
    def add_functional_annotations(self, gff_file: Union[str, Path], 
                                 annotation_type: str = "genes") -> bool:
        """
        Load functional annotations from GFF file
        
        Args:
            gff_file: Path to GFF3 file with gene/feature annotations
            annotation_type: Type of annotations (genes, regulatory, etc.)
            
        Returns:
            Success status
        """
        try:
            logger.info(f"Loading {annotation_type} annotations from {gff_file}")
            
            # Parse GFF file
            annotations = self._parse_gff_file(gff_file)
            
            if not annotations.empty:  # Fixed: use .empty instead of bool(annotations)
                self.functional_data[annotation_type] = annotations
                
                # Register in registry
                self.registry.register_file(
                    f'{annotation_type}_annotations',
                    annotations,
                    'csv',
                    f'{annotation_type.title()} functional annotations',
                    metadata={
                        'source': str(gff_file),
                        'annotation_type': annotation_type,
                        'total_features': len(annotations)
                    }
                )
                
                logger.info(f"Loaded {len(annotations)} {annotation_type} annotations")
                return True
            else:
                logger.warning(f"No annotations found in {gff_file}")
                return False
                
        except Exception as e:
            logger.error(f"Failed to load annotations from {gff_file}: {e}")
            return False
    
    def add_repeat_data(self, repeat_file: Union[str, Path], 
                       file_format: str = "bed") -> bool:
        """
        Load repeat element data
        
        Args:
            repeat_file: Path to repeat element file (BED or GFF format)
            file_format: Format of the file ('bed' or 'gff')
            
        Returns:
            Success status
        """
        try:
            logger.info(f"Loading repeat data from {repeat_file}")
            
            if file_format.lower() == 'bed':
                repeats = self._parse_bed_file(repeat_file)
            elif file_format.lower() == 'gff':
                repeats = self._parse_gff_file(repeat_file)
            else:
                raise ValueError(f"Unsupported format: {file_format}")
            
            if not repeats.empty:  # Fixed: use .empty instead of bool(repeats)
                self.repeat_data['repeats'] = repeats
                
                # Register in registry
                self.registry.register_file(
                    'repeat_elements',
                    repeats,
                    'csv',
                    'Repeat element annotations',
                    metadata={
                        'source': str(repeat_file),
                        'format': file_format,
                        'total_repeats': len(repeats)
                    }
                )
                
                logger.info(f"Loaded {len(repeats)} repeat elements")
                return True
            else:
                logger.warning(f"No repeat data found in {repeat_file}")
                return False
                
        except Exception as e:
            logger.error(f"Failed to load repeat data from {repeat_file}: {e}")
            return False
    
    def calculate_gc_content(self, fasta_file: Union[str, Path], 
                           window_size: int = 1000) -> bool:
        """
        Calculate GC content across genome in windows
        
        Args:
            fasta_file: Path to genome FASTA file
            window_size: Window size for GC calculation
            
        Returns:
            Success status
        """
        try:
            logger.info(f"Calculating GC content from {fasta_file} (window: {window_size})")
            
            from Bio import SeqIO
            
            gc_windows = []
            
            for record in SeqIO.parse(fasta_file, 'fasta'):
                seq_str = str(record.seq).upper()
                chrom = record.id
                
                # Calculate GC in sliding windows
                for i in range(0, len(seq_str) - window_size + 1, window_size):
                    window_seq = seq_str[i:i + window_size]
                    
                    if len(window_seq) >= window_size // 2:  # Only count substantial windows
                        gc_count = window_seq.count('G') + window_seq.count('C')
                        total_bases = len(window_seq) - window_seq.count('N')
                        
                        if total_bases > 0:
                            gc_content = gc_count / total_bases
                            
                            gc_windows.append({
                                'chrom': chrom,
                                'start': i,
                                'end': i + window_size,
                                'gc_content': gc_content,
                                'window_size': len(window_seq),
                                'n_count': window_seq.count('N')
                            })
            
            if gc_windows:
                gc_df = pd.DataFrame(gc_windows)
                self.gc_content['windows'] = gc_df
                
                # Register in registry
                self.registry.register_file(
                    'gc_content_windows',
                    gc_df,
                    'csv',
                    f'GC content in {window_size}bp windows',
                    metadata={
                        'source': str(fasta_file),
                        'window_size': window_size,
                        'total_windows': len(gc_df)
                    }
                )
                
                logger.info(f"Calculated GC content for {len(gc_windows)} windows")
                return True
            else:
                logger.warning("No GC content windows calculated")
                return False
                
        except ImportError:
            logger.error("BioPython required for GC content calculation")
            return False
        except Exception as e:
            logger.error(f"Failed to calculate GC content: {e}")
            return False
    
    def compute_inversion_rate_per_mb(self, inversion_df: pd.DataFrame, 
                                    assembly_stats: Dict) -> Dict:
        """
        Calculate inversion rate normalized by genome size
        
        Args:
            inversion_df: DataFrame with inversion data
            assembly_stats: Dictionary with genome statistics
            
        Returns:
            Dictionary with rate metrics
        """
        logger.info("Computing inversion rate per Mb")
        
        total_inversions = len(inversion_df)
        genome_size_mb = assembly_stats.get('total_length', 0) / 1_000_000
        
        if genome_size_mb == 0:
            logger.warning("Genome size not available, cannot calculate rate per Mb")
            return {}
        
        rate_per_mb = total_inversions / genome_size_mb
        
        # Calculate rates by inversion type if available
        rates_by_type = {}
        if 'inversion_type' in inversion_df.columns:
            for inv_type in inversion_df['inversion_type'].unique():
                type_count = len(inversion_df[inversion_df['inversion_type'] == inv_type])
                rates_by_type[inv_type] = type_count / genome_size_mb
        
        # Calculate rates by chromosome
        rates_by_chr = {}
        if 'first_chr' in inversion_df.columns:
            for chrom in inversion_df['first_chr'].unique():
                chr_count = len(inversion_df[inversion_df['first_chr'] == chrom])
                rates_by_chr[chrom] = chr_count / genome_size_mb
        
        rate_metrics = {
            'total_inversions': total_inversions,
            'genome_size_mb': genome_size_mb,
            'inversions_per_mb': rate_per_mb,
            'rates_by_type': rates_by_type,
            'rates_by_chromosome': rates_by_chr
        }
        
        # Register rate metrics
        self.registry.register_file(
            'inversion_rate_metrics',
            rate_metrics,
            'json',
            'Inversion rate calculations per Mb',
            dependencies=['inversion_events'],
            metadata={'rate_per_mb': rate_per_mb}
        )
        
        logger.info(f"Inversion rate: {rate_per_mb:.3f} inversions per Mb")
        return rate_metrics
    
    def analyze_functional_overlap(self, inversion_df: pd.DataFrame, 
                                 overlap_type: str = "genes") -> pd.DataFrame:
        """
        Analyze overlap between inversions and functional elements
        
        Args:
            inversion_df: DataFrame with inversion coordinates
            overlap_type: Type of functional elements to check
            
        Returns:
            DataFrame with overlap analysis
        """
        logger.info(f"Analyzing inversion overlap with {overlap_type}")
        
        if overlap_type not in self.functional_data:
            logger.warning(f"No {overlap_type} data loaded")
            return pd.DataFrame()
        
        annotations = self.functional_data[overlap_type]
        overlaps = []
        
        for _, inversion in inversion_df.iterrows():
            # For now, use chromosome-level overlap
            # In real implementation, would use actual breakpoint coordinates
            inv_chr = inversion.get('first_chr', '')
            
            # Find overlapping annotations on same chromosome
            chr_annotations = annotations[annotations['chrom'] == inv_chr]
            
            # Convert list to string for hashable storage
            feature_types_list = chr_annotations['feature_type'].unique().tolist() if 'feature_type' in chr_annotations.columns else []
            feature_types_str = ';'.join(feature_types_list)  # Convert list to semicolon-separated string
            
            overlap_info = {
                'inversion_id': inversion.get('start_gene', f"inv_{len(overlaps)}"),
                'chromosome': inv_chr,
                'inversion_type': inversion.get('inversion_type', 'unknown'),
                'overlapping_features': len(chr_annotations),
                'feature_types': feature_types_str,  # Store as string instead of list
                'has_overlap': len(chr_annotations) > 0
            }
            
            # Add specific overlap details if coordinates available
            if all(col in inversion for col in ['first_start', 'first_end']):
                inv_start = inversion['first_start']
                inv_end = inversion['first_end']
                
                # Find features that actually overlap with inversion region
                overlapping = chr_annotations[
                    (chr_annotations['start'] < inv_end) & 
                    (chr_annotations['end'] > inv_start)
                ]
                
                overlap_info.update({
                    'precise_overlaps': len(overlapping),
                    'overlap_fraction': len(overlapping) / len(chr_annotations) if len(chr_annotations) > 0 else 0
                })
            
            overlaps.append(overlap_info)
        
        overlap_df = pd.DataFrame(overlaps)
        
        if not overlap_df.empty:
            # Register overlap analysis
            self.registry.register_file(
                f'inversion_{overlap_type}_overlap',
                overlap_df,
                'csv',
                f'Analysis of inversion overlap with {overlap_type}',
                dependencies=['inversion_events', f'{overlap_type}_annotations'],
                metadata={
                    'overlap_type': overlap_type,
                    'total_inversions': len(overlap_df),
                    'inversions_with_overlap': len(overlap_df[overlap_df['has_overlap']])
                }
            )
            
            logger.info(f"Found overlaps for {len(overlap_df[overlap_df['has_overlap']])} inversions")
        
        return overlap_df
    
    def analyze_repeat_correlation(self, inversion_df: pd.DataFrame) -> Dict:
        """
        Analyze correlation between inversions and repeat elements
        
        Args:
            inversion_df: DataFrame with inversion data
            
        Returns:
            Dictionary with correlation analysis
        """
        logger.info("Analyzing inversion-repeat correlation")
        
        if 'repeats' not in self.repeat_data:
            logger.warning("No repeat data loaded")
            return {}
        
        repeats = self.repeat_data['repeats']
        
        # Calculate repeat density per chromosome
        repeat_density = {}
        for chrom in repeats['chrom'].unique():
            chr_repeats = repeats[repeats['chrom'] == chrom]
            total_repeat_length = chr_repeats['length'].sum() if 'length' in chr_repeats.columns else len(chr_repeats)
            repeat_density[chrom] = total_repeat_length
        
        # Calculate inversion density per chromosome
        inversion_density = {}
        for chrom in inversion_df['first_chr'].unique():
            chr_inversions = inversion_df[inversion_df['first_chr'] == chrom]
            inversion_density[chrom] = len(chr_inversions)
        
        # Calculate correlation
        common_chroms = set(repeat_density.keys()) & set(inversion_density.keys())
        
        if len(common_chroms) > 2:
            repeat_values = [repeat_density[chrom] for chrom in common_chroms]
            inversion_values = [inversion_density[chrom] for chrom in common_chroms]
            
            correlation = np.corrcoef(repeat_values, inversion_values)[0, 1]
            
            correlation_data = {
                'repeat_inversion_correlation': correlation,
                'chromosomes_analyzed': list(common_chroms),
                'repeat_density_by_chr': repeat_density,
                'inversion_density_by_chr': inversion_density
            }
            
            # Register correlation analysis
            self.registry.register_file(
                'repeat_inversion_correlation',
                correlation_data,
                'json',
                'Correlation analysis between repeats and inversions',
                dependencies=['inversion_events', 'repeat_elements'],
                metadata={'correlation': correlation}
            )
            
            logger.info(f"Repeat-inversion correlation: {correlation:.3f}")
            return correlation_data
        else:
            logger.warning("Insufficient data for correlation analysis")
            return {}
    
    def analyze_gc_correlation(self, inversion_df: pd.DataFrame) -> Dict:
        """
        Analyze correlation between inversions and GC content
        
        Args:
            inversion_df: DataFrame with inversion data
            
        Returns:
            Dictionary with GC correlation analysis
        """
        logger.info("Analyzing inversion-GC content correlation")
        
        if 'windows' not in self.gc_content:
            logger.warning("No GC content data available")
            return {}
        
        gc_windows = self.gc_content['windows']
        
        # For each inversion, find overlapping GC windows and calculate average GC
        inversion_gc = []
        
        for _, inversion in inversion_df.iterrows():
            inv_chr = inversion.get('first_chr', '')
            
            # Get GC windows for this chromosome
            chr_gc = gc_windows[gc_windows['chrom'] == inv_chr]
            
            if len(chr_gc) > 0:
                avg_gc = chr_gc['gc_content'].mean()
                inversion_gc.append({
                    'inversion_id': inversion.get('start_gene', f"inv_{len(inversion_gc)}"),
                    'chromosome': inv_chr,
                    'avg_gc_content': avg_gc,
                    'gc_windows_count': len(chr_gc)
                })
        
        if inversion_gc:
            gc_analysis_df = pd.DataFrame(inversion_gc)
            
            # Calculate overall statistics
            overall_stats = {
                'mean_gc_at_inversions': gc_analysis_df['avg_gc_content'].mean(),
                'std_gc_at_inversions': gc_analysis_df['avg_gc_content'].std(),
                'genome_wide_mean_gc': gc_windows['gc_content'].mean(),
                'genome_wide_std_gc': gc_windows['gc_content'].std()
            }
            
            # Register GC analysis
            self.registry.register_file(
                'inversion_gc_analysis',
                gc_analysis_df,
                'csv',
                'GC content analysis at inversion sites',
                dependencies=['inversion_events', 'gc_content_windows'],
                metadata=overall_stats
            )
            
            logger.info(f"GC content at inversions: {overall_stats['mean_gc_at_inversions']:.3f}")
            return overall_stats
        else:
            logger.warning("No GC data available for inversion analysis")
            return {}
    
    def detect_positional_patterns(self, inversion_df: pd.DataFrame, 
                                 chromosome_info: Dict) -> Dict:
        """
        Detect positional patterns in inversion distribution
        
        Args:
            inversion_df: DataFrame with inversion data
            chromosome_info: Dictionary with chromosome lengths and centromere info
            
        Returns:
            Dictionary with positional pattern analysis
        """
        logger.info("Detecting positional patterns in inversions")
        
        patterns = {
            'telomeric_enrichment': {},
            'centromeric_enrichment': {},
            'positional_distribution': {}
        }
        
        for chrom in inversion_df['first_chr'].unique():
            chr_inversions = inversion_df[inversion_df['first_chr'] == chrom]
            
            if chrom in chromosome_info:
                chr_length = chromosome_info[chrom].get('length', 0)
                centromere_pos = chromosome_info[chrom].get('centromere', chr_length // 2)
                
                if chr_length > 0:
                    # Calculate distances to telomeres
                    telomere_distances = []
                    for _, inv in chr_inversions.iterrows():
                        pos = inv.get('first_start', chr_length // 2)
                        dist_to_start = pos
                        dist_to_end = chr_length - pos
                        min_telomere_dist = min(dist_to_start, dist_to_end)
                        telomere_distances.append(min_telomere_dist / chr_length)  # Normalized
                    
                    # Calculate distances to centromere
                    centromere_distances = []
                    for _, inv in chr_inversions.iterrows():
                        pos = inv.get('first_start', chr_length // 2)
                        cent_dist = abs(pos - centromere_pos) / chr_length  # Normalized
                        centromere_distances.append(cent_dist)
                    
                    patterns['telomeric_enrichment'][chrom] = {
                        'mean_distance_to_telomere': np.mean(telomere_distances) if telomere_distances else 0,
                        'telomeric_inversions': sum(1 for d in telomere_distances if d < 0.1)  # Within 10% of telomere
                    }
                    
                    patterns['centromeric_enrichment'][chrom] = {
                        'mean_distance_to_centromere': np.mean(centromere_distances) if centromere_distances else 0,
                        'centromeric_inversions': sum(1 for d in centromere_distances if d < 0.1)  # Within 10% of centromere
                    }
                    
                    # Divide chromosome into bins and count inversions
                    n_bins = 10
                    bin_size = chr_length / n_bins
                    bin_counts = [0] * n_bins
                    
                    for _, inv in chr_inversions.iterrows():
                        pos = inv.get('first_start', chr_length // 2)
                        bin_idx = min(int(pos / bin_size), n_bins - 1)
                        bin_counts[bin_idx] += 1
                    
                    patterns['positional_distribution'][chrom] = bin_counts
        
        # Register positional analysis
        self.registry.register_file(
            'positional_patterns',
            patterns,
            'json',
            'Positional pattern analysis of inversions',
            dependencies=['inversion_events'],
            metadata={
                'chromosomes_analyzed': list(patterns['telomeric_enrichment'].keys()),
                'pattern_types': ['telomeric', 'centromeric', 'positional_distribution']
            }
        )
        
        logger.info(f"Analyzed positional patterns for {len(patterns['telomeric_enrichment'])} chromosomes")
        return patterns
    
    def generate_contextual_summary(self, inversion_df: pd.DataFrame, 
                                  assembly_stats: Dict,
                                  chromosome_info: Dict = None) -> Dict:
        """
        Generate comprehensive contextual metrics summary
        
        Args:
            inversion_df: DataFrame with inversion data
            assembly_stats: Assembly statistics
            chromosome_info: Chromosome information
            
        Returns:
            Complete contextual analysis summary
        """
        logger.info("Generating comprehensive contextual metrics summary")
        
        summary = {
            'rate_metrics': {},
            'functional_overlaps': {},
            'repeat_correlation': {},
            'gc_correlation': {},
            'positional_patterns': {}
        }
        
        # 1. Inversion rate metrics
        summary['rate_metrics'] = self.compute_inversion_rate_per_mb(inversion_df, assembly_stats)
        
        # 2. Functional overlaps (if data available)
        for annotation_type in self.functional_data.keys():
            overlap_df = self.analyze_functional_overlap(inversion_df, annotation_type)
            if not overlap_df.empty:
                summary['functional_overlaps'][annotation_type] = {
                    'total_inversions_analyzed': len(overlap_df),
                    'inversions_with_overlap': len(overlap_df[overlap_df['has_overlap']]),
                    'overlap_percentage': len(overlap_df[overlap_df['has_overlap']]) / len(overlap_df) * 100
                }
        
        # 3. Repeat correlation (if data available)
        if self.repeat_data:
            summary['repeat_correlation'] = self.analyze_repeat_correlation(inversion_df)
        
        # 4. GC correlation (if data available)
        if self.gc_content:
            summary['gc_correlation'] = self.analyze_gc_correlation(inversion_df)
        
        # 5. Positional patterns (if chromosome info available)
        if chromosome_info:
            summary['positional_patterns'] = self.detect_positional_patterns(inversion_df, chromosome_info)
        
        # Register complete summary
        self.registry.register_file(
            'contextual_metrics_summary',
            summary,
            'json',
            'Comprehensive contextual metrics analysis',
            dependencies=['inversion_events'],
            metadata={
                'metrics_computed': list(summary.keys()),
                'has_functional_data': bool(self.functional_data),
                'has_repeat_data': bool(self.repeat_data),
                'has_gc_data': bool(self.gc_content)
            }
        )
        
        logger.info("Contextual metrics summary generated")
        return summary
    
    def _parse_gff_file(self, gff_file: Path) -> pd.DataFrame:
        """Parse GFF3 file and return DataFrame"""
        annotations = []
        
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    annotation = {
                        'chrom': fields[0],
                        'source': fields[1],
                        'feature_type': fields[2],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'score': fields[5] if fields[5] != '.' else None,
                        'strand': fields[6],
                        'phase': fields[7] if fields[7] != '.' else None,
                        'attributes': fields[8]
                    }
                    annotations.append(annotation)
        
        return pd.DataFrame(annotations)
    
    def _parse_bed_file(self, bed_file: Path) -> pd.DataFrame:
        """Parse BED file and return DataFrame"""
        bed_data = []
        
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    bed_entry = {
                        'chrom': fields[0],
                        'start': int(fields[1]),
                        'end': int(fields[2]),
                        'length': int(fields[2]) - int(fields[1])
                    }
                    
                    # Add optional fields if present
                    if len(fields) > 3:
                        bed_entry['name'] = fields[3]
                    if len(fields) > 4:
                        bed_entry['score'] = fields[4]
                    if len(fields) > 5:
                        bed_entry['strand'] = fields[5]
                    
                    bed_data.append(bed_entry)
        
        return pd.DataFrame(bed_data)