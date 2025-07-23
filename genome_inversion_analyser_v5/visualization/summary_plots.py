# =============================================================================
# Summary Plotter (summary_plots.py)
# =============================================================================

"""
Comprehensive summary plotting system for analysis results.
Creates publication-quality summary plots and statistical visualizations.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional
from pathlib import Path

from ..logger import get_logger

logger = get_logger()

class SummaryPlotter:
    """
    Comprehensive summary plotter for genomic analysis results.
    Creates publication-quality statistical plots and summaries.
    """
    
    def __init__(self, config):
        """
        Initialize summary plotter.
        
        Args:
            config: Configuration object with plotting parameters
        """
        self.config = config
        self.figure_size = (15, 12)
        self.dpi = config.get('dpi', 300)
        self.color_palette = config.get('color_palette', 'viridis')
        
        # Set style
        plt.style.use('default')
        sns.set_palette(self.color_palette)
        plt.rcParams['figure.dpi'] = self.dpi
        plt.rcParams['font.size'] = config.get('font_size', 12)
        
        logger.info("Summary plotter initialized")
    
    def create_analysis_dashboard(self, results: Dict[str, Any], 
                                output_path: str = None) -> str:
        """
        Create comprehensive analysis dashboard with key metrics.
        
        Args:
            results: Dictionary with analysis results
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        logger.info("Creating comprehensive analysis dashboard...")
        
        fig = plt.figure(figsize=self.figure_size)
        gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
        fig.suptitle('Genome Comparison Analysis Dashboard', fontsize=16, fontweight='bold')
        
        # Panel 1: Overall statistics
        ax1 = fig.add_subplot(gs[0, :2])
        self._plot_overall_statistics(ax1, results)
        
        # Panel 2: Quality scores
        ax2 = fig.add_subplot(gs[0, 2:])
        self._plot_quality_summary(ax2, results)
        
        # Panel 3: Chromosome mapping overview
        ax3 = fig.add_subplot(gs[1, :])
        self._plot_chromosome_mapping_overview(ax3, results)
        
        # Panel 4: Individual component summaries
        ax4 = fig.add_subplot(gs[2, 0])
        self._plot_synteny_summary(ax4, results)
        
        ax5 = fig.add_subplot(gs[2, 1])
        self._plot_inversion_summary(ax5, results)
        
        ax6 = fig.add_subplot(gs[2, 2])
        self._plot_rearrangement_summary(ax6, results)
        
        ax7 = fig.add_subplot(gs[2, 3])
        self._plot_method_summary(ax7, results)
        
        # Save plot
        if output_path is None:
            output_path = "analysis_dashboard.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Analysis dashboard saved: {output_path}")
        return output_path
    
    def create_ortholog_quality_plots(self, ortholog_df: pd.DataFrame,
                                    output_path: str = None) -> str:
        """
        Create comprehensive ortholog quality assessment plots.
        
        Args:
            ortholog_df: DataFrame with ortholog data
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        if len(ortholog_df) == 0:
            logger.warning("No ortholog data for quality plots")
            return None
        
        logger.info("Creating ortholog quality assessment plots...")
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('Ortholog Quality Assessment', fontsize=14, fontweight='bold')
        
        # 1. Similarity distribution
        if 'similarity' in ortholog_df.columns:
            axes[0, 0].hist(ortholog_df['similarity'], bins=30, alpha=0.7, 
                          edgecolor='black', color='skyblue')
            axes[0, 0].axvline(ortholog_df['similarity'].mean(), color='red', 
                             linestyle='--', label=f'Mean: {ortholog_df["similarity"].mean():.3f}')
            axes[0, 0].set_xlabel('Similarity Score')
            axes[0, 0].set_ylabel('Number of Orthologs')
            axes[0, 0].set_title('Similarity Score Distribution')
            axes[0, 0].legend()
            axes[0, 0].grid(True, alpha=0.3)
        
        # 2. Confidence distribution
        if 'confidence' in ortholog_df.columns:
            axes[0, 1].hist(ortholog_df['confidence'], bins=30, alpha=0.7, 
                          edgecolor='black', color='lightgreen')
            axes[0, 1].axvline(ortholog_df['confidence'].mean(), color='red', 
                             linestyle='--', label=f'Mean: {ortholog_df["confidence"].mean():.3f}')
            axes[0, 1].set_xlabel('Confidence Score')
            axes[0, 1].set_ylabel('Number of Orthologs')
            axes[0, 1].set_title('Confidence Score Distribution')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # 3. Identity vs Coverage scatter
        if 'identity' in ortholog_df.columns and 'coverage' in ortholog_df.columns:
            scatter = axes[0, 2].scatter(ortholog_df['identity'], ortholog_df['coverage'],
                                      alpha=0.6, s=10, c=ortholog_df.get('confidence', 'blue'))
            axes[0, 2].set_xlabel('Sequence Identity')
            axes[0, 2].set_ylabel('Coverage')
            axes[0, 2].set_title('Identity vs Coverage')
            axes[0, 2].grid(True, alpha=0.3)
            
            if 'confidence' in ortholog_df.columns:
                cbar = plt.colorbar(scatter, ax=axes[0, 2])
                cbar.set_label('Confidence')
        
        # 4. Gene length distribution
        if 'first_length' in ortholog_df.columns:
            axes[1, 0].hist(ortholog_df['first_length'], bins=30, alpha=0.7,
                          edgecolor='black', color='orange')
            axes[1, 0].set_xlabel('Gene Length (bp)')
            axes[1, 0].set_ylabel('Number of Genes')
            axes[1, 0].set_title('Gene Length Distribution')
            axes[1, 0].grid(True, alpha=0.3)
        
        # 5. Alignment method usage
        if 'alignment_method' in ortholog_df.columns:
            method_counts = ortholog_df['alignment_method'].value_counts()
            colors = plt.cm.Set3(range(len(method_counts)))
            axes[1, 1].pie(method_counts.values, labels=method_counts.index,
                         autopct='%1.1f%%', colors=colors)
            axes[1, 1].set_title('Alignment Methods Used')
        
        # 6. Length ratio vs similarity
        if 'length_ratio' in ortholog_df.columns and 'similarity' in ortholog_df.columns:
            axes[1, 2].scatter(ortholog_df['length_ratio'], ortholog_df['similarity'],
                             alpha=0.6, s=10)
            axes[1, 2].set_xlabel('Length Ratio')
            axes[1, 2].set_ylabel('Similarity')
            axes[1, 2].set_title('Length Ratio vs Similarity')
            axes[1, 2].grid(True, alpha=0.3)
        
        # Save plot
        if output_path is None:
            output_path = "ortholog_quality_assessment.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Ortholog quality plots saved: {output_path}")
        return output_path
    
    def create_synteny_analysis_plots(self, synteny_df: pd.DataFrame,
                                    output_path: str = None) -> str:
        """
        Create comprehensive synteny analysis plots.
        
        Args:
            synteny_df: DataFrame with synteny blocks
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        if len(synteny_df) == 0:
            logger.warning("No synteny data for analysis plots")
            return None
        
        logger.info("Creating synteny analysis plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Synteny Block Analysis', fontsize=14, fontweight='bold')
        
        # 1. Block size distribution
        axes[0, 0].hist(synteny_df['block_size'], bins=20, alpha=0.7, 
                       edgecolor='black', color='orange')
        axes[0, 0].set_xlabel('Block Size (genes)')
        axes[0, 0].set_ylabel('Number of Blocks')
        axes[0, 0].set_title('Synteny Block Size Distribution')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Add statistics
        mean_size = synteny_df['block_size'].mean()
        axes[0, 0].axvline(mean_size, color='red', linestyle='--', 
                          label=f'Mean: {mean_size:.1f}')
        axes[0, 0].legend()
        
        # 2. Synteny types
        if 'synteny_type' in synteny_df.columns:
            type_counts = synteny_df['synteny_type'].value_counts()
            colors = plt.cm.Pastel1(range(len(type_counts)))
            axes[0, 1].pie(type_counts.values, labels=type_counts.index, 
                          autopct='%1.1f%%', colors=colors)
            axes[0, 1].set_title('Synteny Types')
        
        # 3. Position correlation distribution
        if 'position_correlation' in synteny_df.columns:
            axes[1, 0].hist(synteny_df['position_correlation'], bins=20, alpha=0.7,
                          edgecolor='black', color='lightcoral')
            axes[1, 0].set_xlabel('Position Correlation')
            axes[1, 0].set_ylabel('Number of Blocks')
            axes[1, 0].set_title('Position Correlation Distribution')
            axes[1, 0].grid(True, alpha=0.3)
            
            # Add vertical lines for significance thresholds
            axes[1, 0].axvline(0.5, color='green', linestyle='--', alpha=0.7, label='Threshold')
            axes[1, 0].axvline(-0.5, color='green', linestyle='--', alpha=0.7)
            axes[1, 0].legend()
        
        # 4. Confidence vs block size scatter
        if 'confidence' in synteny_df.columns:
            scatter = axes[1, 1].scatter(synteny_df['block_size'], synteny_df['confidence'],
                                       alpha=0.7, c=synteny_df['block_size'], cmap='viridis')
            axes[1, 1].set_xlabel('Block Size (genes)')
            axes[1, 1].set_ylabel('Confidence Score')
            axes[1, 1].set_title('Block Size vs Confidence')
            axes[1, 1].grid(True, alpha=0.3)
            
            plt.colorbar(scatter, ax=axes[1, 1], label='Block Size')
        
        # Save plot
        if output_path is None:
            output_path = "synteny_analysis_plots.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Synteny analysis plots saved: {output_path}")
        return output_path
    
    def create_inversion_analysis_plots(self, inversion_df: pd.DataFrame,
                                      output_path: str = None) -> str:
        """
        Create comprehensive inversion analysis plots.
        
        Args:
            inversion_df: DataFrame with inversion data
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        if len(inversion_df) == 0:
            logger.warning("No inversion data for analysis plots")
            return None
        
        logger.info("Creating inversion analysis plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Inversion Analysis', fontsize=14, fontweight='bold')
        
        # 1. Inversion size distribution
        if 'size_genes' in inversion_df.columns:
            axes[0, 0].hist(inversion_df['size_genes'], bins=20, alpha=0.7,
                          edgecolor='black', color='crimson')
            axes[0, 0].set_xlabel('Inversion Size (genes)')
            axes[0, 0].set_ylabel('Number of Inversions')
            axes[0, 0].set_title('Inversion Size Distribution')
            axes[0, 0].grid(True, alpha=0.3)
            
            # Add statistics
            mean_size = inversion_df['size_genes'].mean()
            axes[0, 0].axvline(mean_size, color='blue', linestyle='--',
                              label=f'Mean: {mean_size:.1f}')
            axes[0, 0].legend()
        
        # 2. Inversion types
        if 'inversion_type' in inversion_df.columns:
            type_counts = inversion_df['inversion_type'].value_counts()
            colors = plt.cm.Reds(np.linspace(0.4, 0.8, len(type_counts)))
            bars = axes[0, 1].bar(range(len(type_counts)), type_counts.values, color=colors)
            axes[0, 1].set_xticks(range(len(type_counts)))
            axes[0, 1].set_xticklabels([label.replace('_', '\n') for label in type_counts.index],
                                      rotation=45, ha='right')
            axes[0, 1].set_ylabel('Number of Inversions')
            axes[0, 1].set_title('Inversion Types')
            axes[0, 1].grid(True, alpha=0.3)
            
            # Add value labels
            for bar, value in zip(bars, type_counts.values):
                height = bar.get_height()
                axes[0, 1].text(bar.get_x() + bar.get_width()/2., height + 0.1,
                               f'{value}', ha='center', va='bottom')
        
        # 3. Chromosome distribution of inversions
        if 'first_chr' in inversion_df.columns:
            chr_counts = inversion_df['first_chr'].value_counts()
            colors = plt.cm.Set3(range(len(chr_counts)))
            axes[1, 0].pie(chr_counts.values, labels=chr_counts.index,
                          autopct='%1.1f%%', colors=colors)
            axes[1, 0].set_title('Inversions by Chromosome')
        
        # 4. Confidence distribution
        if 'confidence' in inversion_df.columns:
            axes[1, 1].hist(inversion_df['confidence'], bins=15, alpha=0.7,
                          edgecolor='black', color='purple')
            axes[1, 1].set_xlabel('Confidence Score')
            axes[1, 1].set_ylabel('Number of Inversions')
            axes[1, 1].set_title('Inversion Confidence Distribution')
            axes[1, 1].grid(True, alpha=0.3)
            
            # Add mean line
            mean_conf = inversion_df['confidence'].mean()
            axes[1, 1].axvline(mean_conf, color='red', linestyle='--',
                              label=f'Mean: {mean_conf:.3f}')
            axes[1, 1].legend()
        
        # Save plot
        if output_path is None:
            output_path = "inversion_analysis_plots.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Inversion analysis plots saved: {output_path}")
        return output_path
    
    def _plot_overall_statistics(self, ax, results: Dict[str, Any]):
        """Plot overall analysis statistics."""
        categories = []
        values = []
        colors = []
        
        if 'ortholog_df' in results:
            categories.append('Orthologs')
            values.append(len(results['ortholog_df']))
            colors.append('#1f77b4')
        
        if 'synteny_df' in results:
            categories.append('Synteny\nBlocks')
            values.append(len(results['synteny_df']))
            colors.append('#ff7f0e')
        
        if 'inversion_df' in results:
            categories.append('Inversions')
            values.append(len(results['inversion_df']))
            colors.append('#d62728')
        
        if 'rearrangement_df' in results:
            categories.append('Rearrangements')
            values.append(len(results['rearrangement_df']))
            colors.append('#2ca02c')
        
        if categories:
            bars = ax.bar(categories, values, color=colors, alpha=0.7, edgecolor='black')
            ax.set_ylabel('Count')
            ax.set_title('Overall Analysis Results', fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add value labels on bars
            for bar, value in zip(bars, values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01*max(values),
                       f'{value}', ha='center', va='bottom', fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No statistics available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('Overall Statistics')
    
    def _plot_quality_summary(self, ax, results: Dict[str, Any]):
        """Plot assembly quality summary."""
        quality_data = []
        labels = []
        
        if 'first_quality' in results:
            quality_data.append(results['first_quality']['quality_score'])
            labels.append('First Genome')
        
        if 'second_quality' in results:
            quality_data.append(results['second_quality']['quality_score'])
            labels.append('Second Genome')
        
        if quality_data:
            colors = ['#2E8B57', '#4682B4'][:len(quality_data)]
            bars = ax.bar(labels, quality_data, color=colors, alpha=0.7, edgecolor='black')
            ax.set_ylabel('Quality Score')
            ax.set_title('Assembly Quality Comparison', fontweight='bold')
            ax.set_ylim(0, 1)
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add value labels
            for bar, value in zip(bars, quality_data):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                       f'{value:.3f}', ha='center', va='bottom', fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No quality data available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('Assembly Quality')
    
    def _plot_chromosome_mapping_overview(self, ax, results: Dict[str, Any]):
        """Plot chromosome mapping overview."""
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            ortholog_df = results['ortholog_df']
            
            # Create chromosome mapping matrix
            first_chrs = sorted(ortholog_df['first_chr'].unique())
            second_chrs = sorted(ortholog_df['second_chr'].unique())
            
            # Create matrix of gene counts
            matrix = np.zeros((len(second_chrs), len(first_chrs)))
            
            for _, row in ortholog_df.iterrows():
                first_idx = first_chrs.index(row['first_chr'])
                second_idx = second_chrs.index(row['second_chr'])
                matrix[second_idx, first_idx] += 1
            
            # Create heatmap
            im = ax.imshow(matrix, cmap='Blues', aspect='auto')
            
            # Set ticks and labels
            ax.set_xticks(range(len(first_chrs)))
            ax.set_yticks(range(len(second_chrs)))
            ax.set_xticklabels([chr_name[:10] for chr_name in first_chrs], rotation=45, ha='right')
            ax.set_yticklabels([chr_name[:10] for chr_name in second_chrs])
            
            ax.set_xlabel('First Genome Chromosomes')
            ax.set_ylabel('Second Genome Chromosomes')
            ax.set_title('Chromosome Mapping Matrix', fontweight='bold')
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Number of Genes')
        else:
            ax.text(0.5, 0.5, 'No mapping data available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
            ax.set_title('Chromosome Mapping')
    
    def _plot_synteny_summary(self, ax, results: Dict[str, Any]):
        """Plot synteny analysis summary."""
        if 'synteny_df' in results and len(results['synteny_df']) > 0:
            synteny_df = results['synteny_df']
            
            if 'synteny_type' in synteny_df.columns:
                type_counts = synteny_df['synteny_type'].value_counts()
                colors = plt.cm.Set2(range(len(type_counts)))
                
                wedges, texts, autotexts = ax.pie(type_counts.values, labels=type_counts.index,
                                                 autopct='%1.0f%%', colors=colors, startangle=90)
                
                # Improve text visibility
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
            else:
                ax.text(0.5, 0.5, f'{len(synteny_df)} synteny blocks', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No synteny data', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
        
        ax.set_title('Synteny Analysis', fontweight='bold')
    
    def _plot_inversion_summary(self, ax, results: Dict[str, Any]):
        """Plot inversion analysis summary."""
        if 'inversion_df' in results and len(results['inversion_df']) > 0:
            inversion_df = results['inversion_df']
            
            if 'inversion_type' in inversion_df.columns:
                type_counts = inversion_df['inversion_type'].value_counts()
                colors = plt.cm.Set1(range(len(type_counts)))
                
                bars = ax.bar(range(len(type_counts)), type_counts.values, color=colors, alpha=0.7)
                ax.set_xticks(range(len(type_counts)))
                ax.set_xticklabels([label.replace('_', '\n') for label in type_counts.index], 
                                  rotation=45, ha='right', fontsize=8)
                ax.set_ylabel('Count')
                ax.grid(True, alpha=0.3, axis='y')
                
                # Add value labels
                for bar, value in zip(bars, type_counts.values):
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + 0.01*max(type_counts.values),
                           f'{value}', ha='center', va='bottom', fontsize=8, fontweight='bold')
            else:
                ax.text(0.5, 0.5, f'{len(inversion_df)} inversions', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No inversions detected', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
        
        ax.set_title('Inversion Analysis', fontweight='bold')
    
    def _plot_rearrangement_summary(self, ax, results: Dict[str, Any]):
        """Plot rearrangement analysis summary."""
        if 'rearrangement_df' in results and len(results['rearrangement_df']) > 0:
            rearr_df = results['rearrangement_df']
            
            if 'type' in rearr_df.columns:
                type_counts = rearr_df['type'].value_counts()
                colors = plt.cm.Pastel1(range(len(type_counts)))
                
                wedges, texts, autotexts = ax.pie(type_counts.values, 
                                                 labels=[label.replace('_', '\n') for label in type_counts.index],
                                                 autopct='%1.0f%%', colors=colors, startangle=90)
                
                # Improve text visibility
                for autotext in autotexts:
                    autotext.set_fontweight('bold')
                    autotext.set_fontsize(8)
            else:
                ax.text(0.5, 0.5, f'{len(rearr_df)} rearrangements', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'No rearrangements detected', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
        
        ax.set_title('Rearrangement Analysis', fontweight='bold')
    
    def _plot_method_summary(self, ax, results: Dict[str, Any]):
        """Plot alignment method usage summary."""
        method_data = None
        
        # Try to get method data from ortholog results
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            ortholog_df = results['ortholog_df']
            if 'alignment_method' in ortholog_df.columns:
                method_data = ortholog_df['alignment_method'].value_counts()
        
        if method_data is not None and len(method_data) > 0:
            colors = plt.cm.Set3(range(len(method_data)))
            bars = ax.bar(range(len(method_data)), method_data.values, color=colors, alpha=0.7)
            ax.set_xticks(range(len(method_data)))
            ax.set_xticklabels(method_data.index, rotation=45, ha='right')
            ax.set_ylabel('Number of Alignments')
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add value labels
            for bar, value in zip(bars, method_data.values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01*max(method_data.values),
                       f'{value}', ha='center', va='bottom', fontsize=8, fontweight='bold')
        else:
            ax.text(0.5, 0.5, 'Method data\nnot available', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12)
        
        ax.set_title('Alignment Methods', fontweight='bold')