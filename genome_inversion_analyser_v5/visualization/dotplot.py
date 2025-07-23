# =============================================================================
# Synteny Dotplot Visualizer (dotplot.py)
# =============================================================================

"""
Advanced synteny dotplot visualization with confidence coloring and interactive features.
Creates publication-quality synteny plots with comprehensive annotation.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import seaborn as sns

from ..logger import get_logger

logger = get_logger()

class SyntenyDotplot:
    """
    Advanced synteny dotplot visualizer with confidence coloring and comprehensive annotation.
    Supports multiple plot types and customization options.
    """
    
    def __init__(self, config):
        """
        Initialize synteny dotplot visualizer.
        
        Args:
            config: Configuration object with plotting parameters
        """
        self.config = config
        self.figure_size = config.get('dotplot_size', (12, 10))
        self.dpi = config.get('dpi', 300)
        self.synteny_color = config.get('synteny_color', '#1f77b4')
        self.inversion_color = config.get('inversion_color', '#d62728')
        self.confidence_alpha = config.get('confidence_alpha', True)
        self.show_synteny_blocks = config.get('show_synteny_blocks', True)
        self.show_breakpoints = config.get('show_breakpoints', True)
        self.show_labels = config.get('dotplot_show_labels', False)
        
        # Set matplotlib style
        plt.style.use('default')
        plt.rcParams['figure.dpi'] = self.dpi
        plt.rcParams['font.size'] = config.get('font_size', 12)
        
        logger.info("Synteny dotplot visualizer initialized")
        logger.info(f"  Figure size: {self.figure_size}")
        logger.info(f"  Confidence alpha: {'enabled' if self.confidence_alpha else 'disabled'}")
    
    def create_busco_synteny_dotplot(self, ortholog_df: pd.DataFrame, 
                                   output_path: str = None,
                                   title: str = "BUSCO Synteny Analysis") -> str:
        """
        Create BUSCO-based synteny dotplot showing gene order relationships.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            output_path: Optional output file path
            title: Plot title
            
        Returns:
            Path to saved plot file
        """
        if len(ortholog_df) == 0:
            logger.warning("No ortholog data available for dotplot")
            return None
        
        logger.info(f"Creating BUSCO synteny dotplot with {len(ortholog_df)} orthologs...")
        
        # Prepare data
        plot_data = self._prepare_dotplot_data(ortholog_df)
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        # Plot syntenic and inverted points
        self._plot_ortholog_points(ax, plot_data)
        
        # Add synteny block lines if enabled
        if self.show_synteny_blocks:
            self._add_synteny_block_lines(ax, plot_data)
        
        # Add breakpoints if enabled
        if self.show_breakpoints:
            self._add_breakpoints(ax, plot_data)
        
        # Customize plot
        self._customize_dotplot(ax, plot_data, title)
        
        # Save plot
        if output_path is None:
            output_path = "busco_synteny_dotplot.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Synteny dotplot saved: {output_path}")
        return output_path
    
    def create_chromosome_dotplot(self, ortholog_df: pd.DataFrame,
                                first_chr: str, second_chr: str,
                                output_path: str = None) -> str:
        """
        Create detailed dotplot for a specific chromosome pair.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            first_chr: First chromosome identifier
            second_chr: Second chromosome identifier
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        # Filter data for specific chromosome pair
        chr_data = ortholog_df[
            (ortholog_df['first_chr'] == first_chr) & 
            (ortholog_df['second_chr'] == second_chr)
        ]
        
        if len(chr_data) == 0:
            logger.warning(f"No data for chromosome pair {first_chr} vs {second_chr}")
            return None
        
        logger.info(f"Creating chromosome dotplot: {first_chr} vs {second_chr}")
        
        # Create figure with additional detail
        fig, (ax_main, ax_hist_x, ax_hist_y) = plt.subplots(
            2, 2, figsize=(self.figure_size[0] + 2, self.figure_size[1] + 2),
            gridspec_kw={'height_ratios': [1, 4], 'width_ratios': [4, 1]}
        )
        
        # Remove unused subplot
        ax_hist_x.remove()
        
        # Main dotplot
        plot_data = self._prepare_dotplot_data(chr_data, use_genomic_coords=True)
        self._plot_ortholog_points(ax_main, plot_data)
        
        if self.show_synteny_blocks:
            self._add_synteny_block_lines(ax_main, plot_data)
        
        # Histograms showing gene density
        self._add_gene_density_histograms(ax_hist_y, ax_main, chr_data)
        
        # Customize
        title = f"Chromosome Synteny: {first_chr} vs {second_chr}"
        self._customize_dotplot(ax_main, plot_data, title)
        ax_main.set_xlabel(f'{first_chr} Position (bp)')
        ax_main.set_ylabel(f'{second_chr} Position (bp)')
        
        # Save plot
        if output_path is None:
            output_path = f"chromosome_dotplot_{first_chr}_vs_{second_chr}.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Chromosome dotplot saved: {output_path}")
        return output_path
    
    def create_multi_chromosome_dotplot(self, ortholog_df: pd.DataFrame,
                                      output_path: str = None) -> str:
        """
        Create multi-panel dotplot showing all chromosome pairs.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        # Get unique chromosome combinations
        chr_pairs = ortholog_df.groupby(['first_chr', 'second_chr']).size().reset_index(name='count')
        chr_pairs = chr_pairs.sort_values('count', ascending=False)
        
        n_pairs = len(chr_pairs)
        if n_pairs == 0:
            logger.warning("No chromosome pairs found for multi-panel dotplot")
            return None
        
        logger.info(f"Creating multi-chromosome dotplot for {n_pairs} chromosome pairs...")
        
        # Calculate subplot layout
        n_cols = min(4, n_pairs)
        n_rows = (n_pairs + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4))
        if n_pairs == 1:
            axes = [axes]
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        
        # Create individual dotplots
        for idx, (_, row) in enumerate(chr_pairs.iterrows()):
            if idx >= n_rows * n_cols:
                break
                
            ax_row = idx // n_cols
            ax_col = idx % n_cols
            ax = axes[ax_row, ax_col] if n_rows > 1 else axes[ax_col]
            
            # Filter data for this chromosome pair
            pair_data = ortholog_df[
                (ortholog_df['first_chr'] == row['first_chr']) & 
                (ortholog_df['second_chr'] == row['second_chr'])
            ]
            
            # Create mini dotplot
            plot_data = self._prepare_dotplot_data(pair_data, use_genomic_coords=True)
            self._plot_ortholog_points(ax, plot_data, point_size=2)
            
            # Customize mini plot
            ax.set_title(f"{row['first_chr']} vs {row['second_chr']}\n({row['count']} genes)")
            ax.set_xlabel(f"{row['first_chr']}")
            ax.set_ylabel(f"{row['second_chr']}")
            
            # Format tick labels
            ax.ticklabel_format(style='scientific', axis='both', scilimits=(0,0))
        
        # Remove empty subplots
        for idx in range(n_pairs, n_rows * n_cols):
            ax_row = idx // n_cols
            ax_col = idx % n_cols
            if n_rows > 1:
                axes[ax_row, ax_col].remove()
            else:
                axes[ax_col].remove()
        
        # Save plot
        if output_path is None:
            output_path = "multi_chromosome_dotplot.png"
        
        plt.suptitle("Multi-Chromosome Synteny Overview", fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Multi-chromosome dotplot saved: {output_path}")
        return output_path
    
    def _prepare_dotplot_data(self, ortholog_df: pd.DataFrame, 
                            use_genomic_coords: bool = False) -> Dict[str, Any]:
        """Prepare data for dotplot visualization."""
        df = ortholog_df.copy()
        
        if use_genomic_coords:
            # Use actual genomic coordinates
            df['x_coord'] = df['first_start']
            df['y_coord'] = df['second_start']
        else:
            # Create gene order indices for each chromosome
            first_order = {}
            second_order = {}
            
            for chr_name in df['first_chr'].unique():
                chr_genes = df[df['first_chr'] == chr_name].sort_values('first_start')
                first_order.update({row['busco_id']: i for i, (_, row) in enumerate(chr_genes.iterrows())})
            
            for chr_name in df['second_chr'].unique():
                chr_genes = df[df['second_chr'] == chr_name].sort_values('second_start')
                second_order.update({row['busco_id']: i for i, (_, row) in enumerate(chr_genes.iterrows())})
            
            # Add order indices to dataframe
            df['x_coord'] = df['busco_id'].map(first_order)
            df['y_coord'] = df['busco_id'].map(second_order)
        
        # Remove genes without coordinates
        df = df.dropna(subset=['x_coord', 'y_coord'])
        
        # Classify points by strand orientation
        df['is_syntenic'] = df['first_strand'] == df['second_strand']
        
        return {
            'data': df,
            'syntenic_data': df[df['is_syntenic']],
            'inverted_data': df[~df['is_syntenic']],
            'use_genomic_coords': use_genomic_coords
        }
    
    def _plot_ortholog_points(self, ax, plot_data: Dict[str, Any], point_size: int = 20):
        """Plot ortholog points with appropriate coloring."""
        syntenic_data = plot_data['syntenic_data']
        inverted_data = plot_data['inverted_data']
        
        # Plot syntenic points
        if len(syntenic_data) > 0:
            alpha_vals = self._get_alpha_values(syntenic_data) if self.confidence_alpha else 0.7
            scatter = ax.scatter(
                syntenic_data['x_coord'], syntenic_data['y_coord'],
                c=self.synteny_color, alpha=alpha_vals, s=point_size,
                label='Syntenic', edgecolors='none'
            )
        
        # Plot inverted points
        if len(inverted_data) > 0:
            alpha_vals = self._get_alpha_values(inverted_data) if self.confidence_alpha else 0.7
            scatter = ax.scatter(
                inverted_data['x_coord'], inverted_data['y_coord'],
                c=self.inversion_color, alpha=alpha_vals, s=point_size,
                label='Inverted', edgecolors='none'
            )
    
    def _get_alpha_values(self, data: pd.DataFrame) -> np.ndarray:
        """Calculate alpha values based on confidence scores."""
        if 'confidence' in data.columns:
            # Map confidence to alpha range [0.3, 1.0]
            confidences = data['confidence'].values
            return 0.3 + 0.7 * confidences
        else:
            return np.full(len(data), 0.7)
    
    def _add_synteny_block_lines(self, ax, plot_data: Dict[str, Any]):
        """Add lines connecting synteny blocks."""
        df = plot_data['data']
        
        # Group by chromosome pairs and draw connecting lines
        for (first_chr, second_chr), group in df.groupby(['first_chr', 'second_chr']):
            if len(group) >= 3:  # Only draw lines for substantial blocks
                syntenic_group = group[group['is_syntenic']]
                if len(syntenic_group) >= 2:
                    # Sort by first genome coordinate
                    sorted_group = syntenic_group.sort_values('x_coord')
                    ax.plot(
                        sorted_group['x_coord'], sorted_group['y_coord'],
                        color='gray', alpha=0.4, linewidth=0.8, zorder=0
                    )
    
    def _add_breakpoints(self, ax, plot_data: Dict[str, Any]):
        """Add markers for potential breakpoints."""
        df = plot_data['data']
        
        # Identify potential breakpoints (large gaps in synteny)
        for (first_chr, second_chr), group in df.groupby(['first_chr', 'second_chr']):
            if len(group) >= 5:
                sorted_group = group.sort_values('x_coord')
                
                # Calculate gaps
                x_gaps = np.diff(sorted_group['x_coord'])
                gap_threshold = np.percentile(x_gaps, 90)  # Large gaps
                
                breakpoint_indices = np.where(x_gaps > gap_threshold)[0]
                
                for bp_idx in breakpoint_indices:
                    x_pos = sorted_group['x_coord'].iloc[bp_idx + 1]
                    y_pos = sorted_group['y_coord'].iloc[bp_idx + 1]
                    ax.scatter(x_pos, y_pos, marker='x', s=50, c='red', 
                             alpha=0.8, linewidths=2, zorder=10)
    
    def _add_gene_density_histograms(self, ax_hist, ax_main, chr_data: pd.DataFrame):
        """Add gene density histograms to chromosome dotplot."""
        # Y-axis histogram (second chromosome)
        ax_hist.hist(chr_data['second_start'], bins=20, orientation='horizontal',
                    alpha=0.6, color='lightblue', edgecolor='black')
        ax_hist.set_ylabel('Gene Density')
        ax_hist.set_xlabel('')
        
        # Match y-axis with main plot
        ax_hist.set_ylim(ax_main.get_ylim())
    
    def _customize_dotplot(self, ax, plot_data: Dict[str, Any], title: str):
        """Apply consistent customization to dotplots."""
        # Set title and labels
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        if plot_data['use_genomic_coords']:
            ax.set_xlabel('First Genome Position (bp)')
            ax.set_ylabel('Second Genome Position (bp)')
            # Format axis labels
            ax.ticklabel_format(style='scientific', axis='both', scilimits=(0,0))
        else:
            ax.set_xlabel('Gene Order in First Genome')
            ax.set_ylabel('Gene Order in Second Genome')
        
        # Add legend
        ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Equal aspect ratio for better visualization
        ax.set_aspect('auto')
        
        # Add statistics text box
        self._add_statistics_box(ax, plot_data)
    
    def _add_statistics_box(self, ax, plot_data: Dict[str, Any]):
        """Add statistics text box to the plot."""
        df = plot_data['data']
        
        total_genes = len(df)
        syntenic_genes = len(plot_data['syntenic_data'])
        inverted_genes = len(plot_data['inverted_data'])
        
        syntenic_rate = syntenic_genes / total_genes if total_genes > 0 else 0
        inversion_rate = inverted_genes / total_genes if total_genes > 0 else 0
        
        # Calculate average confidence if available
        if 'confidence' in df.columns:
            avg_confidence = df['confidence'].mean()
            confidence_text = f"Avg Confidence: {avg_confidence:.3f}\n"
        else:
            confidence_text = ""
        
        stats_text = (
            f"Total Genes: {total_genes}\n"
            f"Syntenic: {syntenic_genes} ({syntenic_rate:.1%})\n"
            f"Inverted: {inverted_genes} ({inversion_rate:.1%})\n"
            f"{confidence_text}"
        )
        
        # Add text box
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)