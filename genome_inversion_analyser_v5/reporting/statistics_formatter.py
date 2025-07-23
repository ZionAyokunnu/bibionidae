# =============================================================================
# Statistics Formatter (statistics_formatter.py)
# =============================================================================

"""
Statistical formatting utilities for genome analysis results.
Provides formatted statistical summaries and data aggregation.
"""

from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
import numpy as np
from datetime import datetime

from ..logger import get_logger

logger = get_logger()

class StatisticsFormatter:
    """
    Formats and aggregates statistical data from genome analysis results.
    Creates summary statistics and formatted output for reports.
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        Initialize statistics formatter.
        
        Args:
            config: Optional configuration dictionary
        """
        self.config = config or {}
        self.precision = self.config.get('statistics_precision', 3)
        self.include_percentiles = self.config.get('include_percentiles', True)
        
        logger.info("Statistics formatter initialized")
    
    def format_comprehensive_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create comprehensive statistical summary from analysis results.
        
        Args:
            results: Complete analysis results dictionary
            
        Returns:
            Formatted statistics dictionary
        """
        logger.info("Formatting comprehensive statistics...")
        
        formatted_stats = {
            'summary': self._create_summary_statistics(results),
            'synteny': self._format_synteny_statistics(results),
            'inversions': self._format_inversion_statistics(results),
            'rearrangements': self._format_rearrangement_statistics(results),
            'quality': self._format_quality_statistics(results),
            'timing': self._format_timing_statistics(results)
        }
        
        # Add validation statistics if available
        if 'synteny_validation' in results or 'inversion_validation' in results:
            formatted_stats['validation'] = self._format_validation_statistics(results)
        
        logger.info("  Statistics formatting completed")
        return formatted_stats
    
    def _create_summary_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Create high-level summary statistics."""
        summary = {}
        
        # Count main analysis results
        if 'synteny_df' in results:
            summary['total_synteny_blocks'] = len(results['synteny_df'])
        
        if 'inversion_df' in results:
            summary['total_inversions'] = len(results['inversion_df'])
        
        if 'rearrangement_df' in results:
            summary['total_rearrangements'] = len(results['rearrangement_df'])
        
        if 'ortholog_df' in results:
            summary['total_orthologs'] = len(results['ortholog_df'])
        
        # Quality summary
        if 'first_quality' in results and 'second_quality' in results:
            q1_score = results['first_quality']['quality_score']
            q2_score = results['second_quality']['quality_score']
            summary['average_quality_score'] = round((q1_score + q2_score) / 2, self.precision)
            summary['quality_difference'] = round(abs(q1_score - q2_score), self.precision)
        
        # Analysis coverage
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            ortholog_df = results['ortholog_df']
            summary['chromosomes_analyzed'] = {
                'first_genome': len(ortholog_df['first_chr'].unique()),
                'second_genome': len(ortholog_df['second_chr'].unique())
            }
        
        return summary
    
    def _format_synteny_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format synteny analysis statistics."""
        stats = {}
        
        if 'synteny_df' not in results or len(results['synteny_df']) == 0:
            return {'total_blocks': 0, 'message': 'No synteny blocks detected'}
        
        synteny_df = results['synteny_df']
        stats['total_blocks'] = len(synteny_df)
        
        # Block size statistics
        if 'block_size' in synteny_df.columns:
            block_sizes = synteny_df['block_size']
            stats['block_size_stats'] = self._calculate_descriptive_stats(block_sizes, 'genes')
        
        # Synteny type distribution
        if 'synteny_type' in synteny_df.columns:
            type_counts = synteny_df['synteny_type'].value_counts()
            stats['type_distribution'] = type_counts.to_dict()
            stats['type_percentages'] = (type_counts / len(synteny_df) * 100).round(1).to_dict()
        
        # Confidence statistics
        if 'confidence' in synteny_df.columns:
            confidence_scores = synteny_df['confidence']
            stats['confidence_stats'] = self._calculate_descriptive_stats(confidence_scores, 'score')
            stats['high_confidence_blocks'] = int((confidence_scores >= 0.8).sum())
            stats['high_confidence_percentage'] = round(
                (confidence_scores >= 0.8).sum() / len(synteny_df) * 100, 1
            )
        
        # Chromosome distribution
        if 'first_chr' in synteny_df.columns:
            chr_counts = synteny_df['first_chr'].value_counts()
            stats['chromosome_distribution'] = {
                'most_blocks': {'chromosome': chr_counts.index[0], 'count': int(chr_counts.iloc[0])},
                'least_blocks': {'chromosome': chr_counts.index[-1], 'count': int(chr_counts.iloc[-1])},
                'chromosomes_with_blocks': len(chr_counts)
            }
        
        return stats
    
    def _format_inversion_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format inversion analysis statistics."""
        stats = {}
        
        if 'inversion_df' not in results or len(results['inversion_df']) == 0:
            return {'total_inversions': 0, 'message': 'No inversions detected'}
        
        inversion_df = results['inversion_df']
        stats['total_inversions'] = len(inversion_df)
        
        # Size statistics
        if 'size_genes' in inversion_df.columns:
            sizes = inversion_df['size_genes']
            stats['size_stats'] = self._calculate_descriptive_stats(sizes, 'genes')
            stats['total_genes_inverted'] = int(sizes.sum())
        
        # Type distribution
        if 'inversion_type' in inversion_df.columns:
            type_counts = inversion_df['inversion_type'].value_counts()
            stats['type_distribution'] = type_counts.to_dict()
            stats['type_percentages'] = (type_counts / len(inversion_df) * 100).round(1).to_dict()
        
        # Detection method statistics
        if 'detection_method' in inversion_df.columns:
            method_counts = inversion_df['detection_method'].value_counts()
            stats['detection_methods'] = method_counts.to_dict()
        
        # Confidence statistics
        if 'confidence' in inversion_df.columns:
            confidence_scores = inversion_df['confidence']
            stats['confidence_stats'] = self._calculate_descriptive_stats(confidence_scores, 'score')
            stats['high_confidence_inversions'] = int((confidence_scores >= 0.8).sum())
            stats['high_confidence_percentage'] = round(
                (confidence_scores >= 0.8).sum() / len(inversion_df) * 100, 1
            )
        
        # Positional analysis
        if 'first_start' in inversion_df.columns and 'first_end' in inversion_df.columns:
            lengths = inversion_df['first_end'] - inversion_df['first_start']
            stats['positional_stats'] = self._calculate_descriptive_stats(lengths, 'bp')
        
        return stats
    
    def _format_rearrangement_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format rearrangement analysis statistics."""
        stats = {}
        
        if 'rearrangement_df' not in results or len(results['rearrangement_df']) == 0:
            return {'total_rearrangements': 0, 'message': 'No rearrangements detected'}
        
        rearr_df = results['rearrangement_df']
        stats['total_rearrangements'] = len(rearr_df)
        
        # Type distribution
        if 'type' in rearr_df.columns:
            type_counts = rearr_df['type'].value_counts()
            stats['type_distribution'] = type_counts.to_dict()
            stats['type_percentages'] = (type_counts / len(rearr_df) * 100).round(1).to_dict()
            
            # Specific statistics for each type
            for rearr_type in type_counts.index:
                type_data = rearr_df[rearr_df['type'] == rearr_type]
                type_stats = {}
                
                if rearr_type == 'chromosome_split' and 'split_ratio' in type_data.columns:
                    ratios = type_data['split_ratio']
                    type_stats['split_ratio_stats'] = self._calculate_descriptive_stats(ratios, 'ratio')
                
                elif rearr_type == 'chromosome_fusion' and 'fusion_ratio' in type_data.columns:
                    ratios = type_data['fusion_ratio']
                    type_stats['fusion_ratio_stats'] = self._calculate_descriptive_stats(ratios, 'ratio')
                
                elif rearr_type == 'reciprocal_translocation' and 'reciprocity_balance' in type_data.columns:
                    balances = type_data['reciprocity_balance']
                    type_stats['balance_stats'] = self._calculate_descriptive_stats(balances, 'balance')
                
                if type_stats:
                    stats[f'{rearr_type}_details'] = type_stats
        
        # Confidence statistics
        if 'confidence' in rearr_df.columns:
            confidence_scores = rearr_df['confidence']
            stats['confidence_stats'] = self._calculate_descriptive_stats(confidence_scores, 'score')
            stats['high_confidence_rearrangements'] = int((confidence_scores >= 0.8).sum())
            stats['high_confidence_percentage'] = round(
                (confidence_scores >= 0.8).sum() / len(rearr_df) * 100, 1
            )
        
        return stats
    
    def _format_quality_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format assembly quality statistics."""
        stats = {}
        
        # First genome quality
        if 'first_quality' in results:
            stats['first_genome_quality'] = self._format_individual_quality(
                results['first_quality'], 'first_genome'
            )
        
        # Second genome quality
        if 'second_quality' in results:
            stats['second_genome_quality'] = self._format_individual_quality(
                results['second_quality'], 'second_genome'
            )
        
        # Quality comparison
        if 'first_quality' in results and 'second_quality' in results:
            q1 = results['first_quality']
            q2 = results['second_quality']
            
            stats['quality_comparison'] = {
                'score_difference': round(abs(q1['quality_score'] - q2['quality_score']), self.precision),
                'class_match': q1['quality_class'] == q2['quality_class'],
                'better_quality': 'first' if q1['quality_score'] > q2['quality_score'] else 'second',
                'relative_improvement': round(
                    abs(q1['quality_score'] - q2['quality_score']) / min(q1['quality_score'], q2['quality_score']) * 100, 1
                )
            }
        
        return stats
    
    def _format_individual_quality(self, quality_data: Dict[str, Any], genome_name: str) -> Dict[str, Any]:
        """Format individual genome quality statistics."""
        stats = {
            'overall_quality': {
                'score': round(quality_data['quality_score'], self.precision),
                'class': quality_data['quality_class']
            }
        }
        
        # Detailed metrics
        if 'metrics' in quality_data:
            metrics = quality_data['metrics']
            formatted_metrics = {}
            
            # Size metrics
            if 'total_length' in metrics:
                formatted_metrics['total_length_mbp'] = round(metrics['total_length'] / 1e6, 1)
            
            if 'n_contigs' in metrics:
                formatted_metrics['contig_count'] = metrics['n_contigs']
            
            if 'n50' in metrics:
                formatted_metrics['n50_mbp'] = round(metrics['n50'] / 1e6, 2)
            
            # BUSCO metrics
            if 'busco_completeness' in metrics:
                formatted_metrics['busco_completeness_percent'] = round(metrics['busco_completeness'] * 100, 1)
            
            if 'busco_duplication' in metrics:
                formatted_metrics['busco_duplication_percent'] = round(metrics['busco_duplication'] * 100, 1)
            
            stats['detailed_metrics'] = formatted_metrics
        
        return stats
    
    def _format_validation_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format statistical validation results."""
        stats = {}
        
        # Synteny validation
        if 'synteny_validation' in results:
            synteny_val = results['synteny_validation']
            stats['synteny_validation'] = {
                'overall_passed': synteny_val.get('overall_validation', {}).get('validated', False),
                'confidence_level': synteny_val.get('overall_validation', {}).get('confidence_level', 'unknown')
            }
            
            if 'correlation_validation' in synteny_val:
                corr_val = synteny_val['correlation_validation']
                stats['synteny_validation']['correlation'] = {
                    'mean_correlation': round(corr_val.get('mean_correlation', 0), self.precision),
                    'p_value': corr_val.get('p_value', 1.0),
                    'significance': 'significant' if corr_val.get('p_value', 1.0) < 0.05 else 'not_significant'
                }
        
        # Inversion validation
        if 'inversion_validation' in results:
            inversion_val = results['inversion_validation']
            stats['inversion_validation'] = {
                'overall_passed': inversion_val.get('overall_validation', {}).get('validated', False),
                'confidence_level': inversion_val.get('overall_validation', {}).get('confidence_level', 'unknown')
            }
        
        return stats
    
    def _format_timing_statistics(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Format analysis timing statistics."""
        stats = {}
        
        if 'timing' in results:
            timing_data = results['timing']
            
            # Total analysis time
            if 'total_time' in timing_data:
                stats['total_analysis_time_seconds'] = round(timing_data['total_time'], 2)
                stats['total_analysis_time_formatted'] = self._format_duration(timing_data['total_time'])
            
            # Individual step timings
            step_timings = {}
            for step, duration in timing_data.items():
                if step != 'total_time' and isinstance(duration, (int, float)):
                    step_timings[step] = {
                        'seconds': round(duration, 2),
                        'formatted': self._format_duration(duration),
                        'percentage': round(duration / timing_data.get('total_time', duration) * 100, 1)
                    }
            
            if step_timings:
                stats['step_timings'] = step_timings
        
        return stats
    
    def _calculate_descriptive_stats(self, data: pd.Series, unit: str) -> Dict[str, Any]:
        """Calculate descriptive statistics for a data series."""
        if len(data) == 0:
            return {'count': 0, 'unit': unit}
        
        stats = {
            'count': len(data),
            'mean': round(data.mean(), self.precision),
            'median': round(data.median(), self.precision),
            'std': round(data.std(), self.precision),
            'min': round(data.min(), self.precision),
            'max': round(data.max(), self.precision),
            'unit': unit
        }
        
        # Add percentiles if enabled
        if self.include_percentiles:
            percentiles = [25, 75, 90, 95]
            for p in percentiles:
                stats[f'p{p}'] = round(data.quantile(p/100), self.precision)
        
        return stats
    
    def _format_duration(self, seconds: float) -> str:
        """Format duration in seconds to human-readable string."""
        if seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = seconds / 60
            return f"{minutes:.1f}m"
        else:
            hours = seconds / 3600
            return f"{hours:.1f}h"
    
    def create_summary_table(self, results: Dict[str, Any]) -> pd.DataFrame:
        """Create summary table of key statistics."""
        summary_data = []
        
        # Basic counts
        categories = [
            ('Synteny Blocks', 'synteny_df'),
            ('Inversions', 'inversion_df'),
            ('Rearrangements', 'rearrangement_df'),
            ('Orthologs', 'ortholog_df')
        ]
        
        for category, key in categories:
            if key in results:
                count = len(results[key])
                summary_data.append({
                    'Category': category,
                    'Count': count,
                    'Status': 'Detected' if count > 0 else 'None found'
                })
        
        # Quality scores
        if 'first_quality' in results and 'second_quality' in results:
            q1 = results['first_quality']
            q2 = results['second_quality']
            
            summary_data.extend([
                {
                    'Category': 'First Genome Quality',
                    'Count': f"{q1['quality_score']:.3f}",
                    'Status': q1['quality_class'].title()
                },
                {
                    'Category': 'Second Genome Quality',
                    'Count': f"{q2['quality_score']:.3f}",
                    'Status': q2['quality_class'].title()
                }
            ])
        
        return pd.DataFrame(summary_data)
    
    def export_formatted_statistics(self, results: Dict[str, Any], output_path: str) -> bool:
        """Export formatted statistics to JSON file."""
        try:
            formatted_stats = self.format_comprehensive_statistics(results)
            
            # Add metadata
            export_data = {
                'metadata': {
                    'generated_at': datetime.now().isoformat(),
                    'formatter_version': '1.0',
                    'precision': self.precision,
                    'include_percentiles': self.include_percentiles
                },
                'statistics': formatted_stats
            }
            
            import json
            with open(output_path, 'w') as f:
                json.dump(export_data, f, indent=2, default=str)
            
            logger.info(f"Formatted statistics exported to: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to export formatted statistics: {e}")
            return False