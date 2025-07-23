# =============================================================================
# Quality Plots (quality_plots.py)
# =============================================================================

"""
Assembly quality visualization system.
Creates comprehensive quality assessment plots and assembly statistics.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional
import seaborn as sns

from ..logger import get_logger

logger = get_logger()

class QualityPlotter:
    """
    Assembly quality visualization system that creates comprehensive
    quality assessment plots and assembly comparison visualizations.
    """
    
    def __init__(self, config):
        """
        Initialize quality plotter.
        
        Args:
            config: Configuration object with plotting parameters
        """
        self.config = config
        self.figure_size = config.get('quality_plot_size', (14, 10))
        self.dpi = config.get('dpi', 300)
        
        # Set style
        plt.style.use('default')
        plt.rcParams['figure.dpi'] = self.dpi
        plt.rcParams['font.size'] = config.get('font_size', 12)
        
        logger.info("Quality plotter initialized")
    
    def create_assembly_quality_comparison(self, first_quality: Dict[str, Any],
                                         second_quality: Dict[str, Any],
                                         output_path: str = None) -> str:
        """
        Create comprehensive assembly quality comparison plots.
        
        Args:
            first_quality: Quality assessment for first assembly
            second_quality: Quality assessment for second assembly
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        logger.info("Creating assembly quality comparison plots...")
        
        fig, axes = plt.subplots(2, 3, figsize=self.figure_size)
        fig.suptitle('Assembly Quality Comparison', fontsize=16, fontweight='bold')
        
        # 1. Overall quality scores
        self._plot_quality_scores(axes[0, 0], first_quality, second_quality)
        
        # 2. Assembly statistics
        self._plot_assembly_statistics(axes[0, 1], first_quality, second_quality)
        
        # 3. BUSCO completeness
        self._plot_busco_comparison(axes[0, 2], first_quality, second_quality)
        
        # 4. N-statistics
        self._plot_n_statistics(axes[1, 0], first_quality, second_quality)
        
        # 5. Contiguity metrics
        self._plot_contiguity_metrics(axes[1, 1], first_quality, second_quality)
        
        # 6. Quality radar chart
        self._plot_quality_radar(axes[1, 2], first_quality, second_quality)
        
        # Save plot
        if output_path is None:
            output_path = "assembly_quality_comparison.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Quality comparison plots saved: {output_path}")
        return output_path
    
    def create_busco_assessment_plot(self, first_busco_df: pd.DataFrame,
                                   second_busco_df: pd.DataFrame,
                                   output_path: str = None) -> str:
        """
        Create detailed BUSCO assessment visualization.
        
        Args:
            first_busco_df: First genome BUSCO data
            second_busco_df: Second genome BUSCO data
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        logger.info("Creating BUSCO assessment plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('BUSCO Completeness Assessment', fontsize=14, fontweight='bold')
        
        # 1. BUSCO status distribution
        self._plot_busco_status_distribution(axes[0, 0], first_busco_df, second_busco_df)
        
        # 2. BUSCO score distributions
        self._plot_busco_score_distributions(axes[0, 1], first_busco_df, second_busco_df)
        
        # 3. Missing BUSCOs comparison
        self._plot_missing_buscos_comparison(axes[1, 0], first_busco_df, second_busco_df)
        
        # 4. Duplication analysis
        self._plot_duplication_analysis(axes[1, 1], first_busco_df, second_busco_df)
        
        # Save plot
        if output_path is None:
            output_path = "busco_assessment.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  BUSCO assessment plots saved: {output_path}")
        return output_path
    
    def _plot_quality_scores(self, ax, first_quality: Dict[str, Any], 
                           second_quality: Dict[str, Any]):
        """Plot overall quality scores comparison."""
        scores = [first_quality['quality_score'], second_quality['quality_score']]
        classes = [first_quality['quality_class'], second_quality['quality_class']]
        labels = ['First Genome', 'Second Genome']
        colors = ['#1f77b4', '#ff7f0e']
        
        bars = ax.bar(labels, scores, color=colors, alpha=0.7, edgecolor='black')
        ax.set_ylabel('Quality Score')
        ax.set_title('Overall Quality Scores', fontweight='bold')
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels and quality classes
        for bar, score, quality_class in zip(bars, scores, classes):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                   f'{score:.3f}', ha='center', va='bottom', fontweight='bold')
            ax.text(bar.get_x() + bar.get_width()/2., height - 0.1,
                   quality_class, ha='center', va='center', 
                   style='italic', fontsize=10)
    
    def _plot_assembly_statistics(self, ax, first_quality: Dict[str, Any],
                                second_quality: Dict[str, Any]):
        """Plot basic assembly statistics."""
        metrics1 = first_quality.get('metrics', {})
        metrics2 = second_quality.get('metrics', {})
        
        stats_to_plot = ['total_length', 'n_contigs', 'max_contig']
        labels = ['Total Length\n(Mbp)', 'Number of\nContigs', 'Longest Contig\n(Mbp)']
        
        values1 = []
        values2 = []
        
        for stat in stats_to_plot:
            val1 = metrics1.get(stat, 0)
            val2 = metrics2.get(stat, 0)
            
            if stat in ['total_length', 'max_contig']:
                val1 = val1 / 1e6  # Convert to Mbp
                val2 = val2 / 1e6
            
            values1.append(val1)
            values2.append(val2)
        
        x = np.arange(len(labels))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, values1, width, label='First Genome', 
                      color='#1f77b4', alpha=0.7)
        bars2 = ax.bar(x + width/2, values2, width, label='Second Genome', 
                      color='#ff7f0e', alpha=0.7)
        
        ax.set_xlabel('Metrics')
        ax.set_ylabel('Values')
        ax.set_title('Assembly Statistics', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                           f'{height:.1f}', ha='center', va='bottom', fontsize=8)
    
    def _plot_busco_comparison(self, ax, first_quality: Dict[str, Any],
                             second_quality: Dict[str, Any]):
        """Plot BUSCO completeness comparison."""
        metrics1 = first_quality.get('metrics', {})
        metrics2 = second_quality.get('metrics', {})
        
        busco_metrics = ['busco_completeness', 'busco_fragmentation', 'busco_duplication', 'busco_missing_rate']
        labels = ['Complete', 'Fragmented', 'Duplicated', 'Missing']
        
        values1 = [metrics1.get(metric, 0) * 100 for metric in busco_metrics]
        values2 = [metrics2.get(metric, 0) * 100 for metric in busco_metrics]
        
        x = np.arange(len(labels))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, values1, width, label='First Genome', 
                      color='#2ca02c', alpha=0.7)
        bars2 = ax.bar(x + width/2, values2, width, label='Second Genome', 
                      color='#d62728', alpha=0.7)
        
        ax.set_xlabel('BUSCO Categories')
        ax.set_ylabel('Percentage (%)')
        ax.set_title('BUSCO Completeness', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                           f'{height:.1f}%', ha='center', va='bottom', fontsize=8)
    
    def _plot_n_statistics(self, ax, first_quality: Dict[str, Any],
                         second_quality: Dict[str, Any]):
        """Plot N-statistics comparison."""
        metrics1 = first_quality.get('metrics', {})
        metrics2 = second_quality.get('metrics', {})
        
        n_stats = ['n50', 'n90']
        labels = ['N50 (Mbp)', 'N90 (Mbp)']
        
        values1 = [metrics1.get(stat, 0) / 1e6 for stat in n_stats]
        values2 = [metrics2.get(stat, 0) / 1e6 for stat in n_stats]
        
        x = np.arange(len(labels))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, values1, width, label='First Genome', 
                      color='#9467bd', alpha=0.7)
        bars2 = ax.bar(x + width/2, values2, width, label='Second Genome', 
                      color='#8c564b', alpha=0.7)
        
        ax.set_xlabel('N-Statistics')
        ax.set_ylabel('Size (Mbp)')
        ax.set_title('Contiguity Metrics', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height + height*0.02,
                           f'{height:.1f}', ha='center', va='bottom', fontsize=8)
    
    def _plot_contiguity_metrics(self, ax, first_quality: Dict[str, Any],
                               second_quality: Dict[str, Any]):
        """Plot contiguity and fragmentation metrics."""
        metrics1 = first_quality.get('metrics', {})
        metrics2 = second_quality.get('metrics', {})
        
        # Calculate log-scale metrics for better visualization
        frag1 = metrics1.get('assembly_fragmentation', 1)
        frag2 = metrics2.get('assembly_fragmentation', 1)
        
        contig_cv1 = metrics1.get('contig_length_cv', 1)
        contig_cv2 = metrics2.get('contig_length_cv', 1)
        
        categories = ['Assembly\nFragmentation', 'Contig Length\nCV']
        values1 = [frag1, contig_cv1]
        values2 = [frag2, contig_cv2]
        
        x = np.arange(len(categories))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, values1, width, label='First Genome', 
                      color='#17becf', alpha=0.7)
        bars2 = ax.bar(x + width/2, values2, width, label='Second Genome', 
                      color='#bcbd22', alpha=0.7)
        
        ax.set_xlabel('Fragmentation Metrics')
        ax.set_ylabel('Values (log scale)')
        ax.set_yscale('log')
        ax.set_title('Assembly Fragmentation', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_quality_radar(self, ax, first_quality: Dict[str, Any],
                          second_quality: Dict[str, Any]):
        """Plot quality radar chart."""
        # Define quality dimensions
        metrics1 = first_quality.get('metrics', {})
        metrics2 = second_quality.get('metrics', {})
        
        # Normalize metrics to 0-1 scale
        dimensions = ['BUSCO\nCompleteness', 'N50\nScore', 'Contiguity\nScore', 'Overall\nQuality']
        
        # Calculate normalized scores
        busco1 = metrics1.get('busco_completeness', 0)
        busco2 = metrics2.get('busco_completeness', 0)
        
        n50_1 = min(1.0, metrics1.get('n50', 0) / 10e6)  # Normalize by 10Mbp
        n50_2 = min(1.0, metrics2.get('n50', 0) / 10e6)
        
        contiguity1 = max(0, 1 - metrics1.get('assembly_fragmentation', 10) / 10)
        contiguity2 = max(0, 1 - metrics2.get('assembly_fragmentation', 10) / 10)
        
        overall1 = first_quality['quality_score']
        overall2 = second_quality['quality_score']
        
        values1 = [busco1, n50_1, contiguity1, overall1]
        values2 = [busco2, n50_2, contiguity2, overall2]
        
        # Create radar chart
        angles = np.linspace(0, 2 * np.pi, len(dimensions), endpoint=False).tolist()
        values1 += values1[:1]  # Complete the circle
        values2 += values2[:1]
        angles += angles[:1]
        
        ax.plot(angles, values1, 'o-', linewidth=2, label='First Genome', color='#1f77b4')
        ax.fill(angles, values1, alpha=0.25, color='#1f77b4')
        ax.plot(angles, values2, 'o-', linewidth=2, label='Second Genome', color='#ff7f0e')
        ax.fill(angles, values2, alpha=0.25, color='#ff7f0e')
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(dimensions)
        ax.set_ylim(0, 1)
        ax.set_title('Quality Profile Comparison', fontweight='bold')
        ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1.0))
        ax.grid(True)
    
    def _plot_busco_status_distribution(self, ax, first_busco_df: pd.DataFrame,
                                      second_busco_df: pd.DataFrame):
        """Plot BUSCO status distribution."""
        status_counts1 = first_busco_df['status'].value_counts()
        status_counts2 = second_busco_df['status'].value_counts()
        
        # Combine all possible statuses
        all_statuses = set(status_counts1.index) | set(status_counts2.index)
        
        percentages1 = [(status_counts1.get(status, 0) / len(first_busco_df)) * 100 
                       for status in all_statuses]
        percentages2 = [(status_counts2.get(status, 0) / len(second_busco_df)) * 100 
                       for status in all_statuses]
        
        x = np.arange(len(all_statuses))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, percentages1, width, label='First Genome', alpha=0.7)
        bars2 = ax.bar(x + width/2, percentages2, width, label='Second Genome', alpha=0.7)
        
        ax.set_xlabel('BUSCO Status')
        ax.set_ylabel('Percentage (%)')
        ax.set_title('BUSCO Status Distribution', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(list(all_statuses), rotation=45)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    def _plot_busco_score_distributions(self, ax, first_busco_df: pd.DataFrame,
                                      second_busco_df: pd.DataFrame):
        """Plot BUSCO score distributions."""
        if 'score' in first_busco_df.columns and 'score' in second_busco_df.columns:
            complete1 = first_busco_df[first_busco_df['status'] == 'Complete']['score']
            complete2 = second_busco_df[second_busco_df['status'] == 'Complete']['score']
            
            ax.hist(complete1, bins=20, alpha=0.7, label='First Genome', density=True)
            ax.hist(complete2, bins=20, alpha=0.7, label='Second Genome', density=True)
            
            ax.set_xlabel('BUSCO Score')
            ax.set_ylabel('Density')
            ax.set_title('Complete BUSCO Score Distribution', fontweight='bold')
            ax.legend()
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'Score data not available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('BUSCO Scores')
    
    def _plot_missing_buscos_comparison(self, ax, first_busco_df: pd.DataFrame,
                                      second_busco_df: pd.DataFrame):
        """Plot missing BUSCOs comparison."""
        missing1 = set(first_busco_df[first_busco_df['status'] == 'Missing']['busco_id'])
        missing2 = set(second_busco_df[second_busco_df['status'] == 'Missing']['busco_id'])
        
        # Calculate overlaps
        common_missing = len(missing1 & missing2)
        unique_missing1 = len(missing1 - missing2)
        unique_missing2 = len(missing2 - missing1)
        
        labels = ['Common\nMissing', 'Unique to\nFirst', 'Unique to\nSecond']
        values = [common_missing, unique_missing1, unique_missing2]
        colors = ['#d62728', '#1f77b4', '#ff7f0e']
        
        bars = ax.bar(labels, values, color=colors, alpha=0.7)
        ax.set_ylabel('Number of BUSCOs')
        ax.set_title('Missing BUSCO Analysis', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, value in zip(bars, values):
            if value > 0:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                       f'{value}', ha='center', va='bottom', fontweight='bold')
    
    def _plot_duplication_analysis(self, ax, first_busco_df: pd.DataFrame,
                                 second_busco_df: pd.DataFrame):
        """Plot duplication analysis."""
        dup1 = len(first_busco_df[first_busco_df['status'] == 'Duplicated'])
        dup2 = len(second_busco_df[second_busco_df['status'] == 'Duplicated'])
        
        total1 = len(first_busco_df)
        total2 = len(second_busco_df)
        
        dup_rate1 = (dup1 / total1) * 100 if total1 > 0 else 0
        dup_rate2 = (dup2 / total2) * 100 if total2 > 0 else 0
        
        labels = ['First Genome', 'Second Genome']
        values = [dup_rate1, dup_rate2]
        colors = ['#9467bd', '#8c564b']
        
        bars = ax.bar(labels, values, color=colors, alpha=0.7)
        ax.set_ylabel('Duplication Rate (%)')
        ax.set_title('BUSCO Duplication Rates', fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.2,
                   f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')

