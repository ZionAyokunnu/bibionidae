# =============================================================================
# Analysis Dashboard (dashboard.py)
# =============================================================================

"""
Comprehensive analysis dashboard that integrates all visualization components.
Creates a unified interface for viewing all analysis results.
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional
from pathlib import Path

from .dotplot import SyntenyDotplot
from .summary_plots import SummaryPlotter
from .network_plots import RearrangementNetworkPlotter
from .quality_plots import QualityPlotter
from ..logger import get_logger

logger = get_logger()

class AnalysisDashboard:
    """
    Comprehensive analysis dashboard that creates a unified visualization
    of all genome comparison results.
    """
    
    def __init__(self, config):
        """
        Initialize analysis dashboard.
        
        Args:
            config: Configuration object with dashboard parameters
        """
        self.config = config
        self.output_dir = Path(config.get('output_directory', 'output'))
        self.create_pdf_report = config.get('create_pdf_report', True)
        self.figure_size = config.get('dashboard_size', (16, 12))
        self.dpi = config.get('dpi', 300)
        
        # Initialize component plotters
        self.dotplot_maker = SyntenyDotplot(config)
        self.summary_plotter = SummaryPlotter(config)
        self.network_plotter = RearrangementNetworkPlotter(config)
        self.quality_plotter = QualityPlotter(config)
        
        # Set matplotlib parameters
        plt.rcParams['figure.dpi'] = self.dpi
        plt.rcParams['font.size'] = config.get('font_size', 12)
        
        logger.info("Analysis dashboard initialized")
        logger.info(f"  PDF report: {'enabled' if self.create_pdf_report else 'disabled'}")
    
    def create_comprehensive_dashboard(self, results: Dict[str, Any],
                                     output_path: str = None) -> Dict[str, str]:
        """
        Create comprehensive analysis dashboard with all visualizations.
        
        Args:
            results: Dictionary with all analysis results
            output_path: Optional output directory path
            
        Returns:
            Dictionary with paths to generated plots
        """
        logger.info("Creating comprehensive analysis dashboard...")
        
        if output_path is None:
            output_path = self.output_dir / 'plots'
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        
        plot_paths = {}
        
        # 1. Main dashboard overview
        dashboard_path = self._create_main_dashboard(results, output_path / 'main_dashboard.png')
        if dashboard_path:
            plot_paths['main_dashboard'] = dashboard_path
        
        # 2. Synteny dotplots
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            dotplot_path = self.dotplot_maker.create_busco_synteny_dotplot(
                results['ortholog_df'], output_path / 'synteny_dotplot.png'
            )
            if dotplot_path:
                plot_paths['synteny_dotplot'] = dotplot_path
        
        # 3. Quality assessment plots
        if 'first_quality' in results and 'second_quality' in results:
            quality_path = self.quality_plotter.create_assembly_quality_comparison(
                results['first_quality'], results['second_quality'], 
                output_path / 'quality_comparison.png'
            )
            if quality_path:
                plot_paths['quality_comparison'] = quality_path
        
        # 4. Ortholog quality plots
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            ortholog_quality_path = self.summary_plotter.create_ortholog_quality_plots(
                results['ortholog_df'], output_path / 'ortholog_quality.png'
            )
            if ortholog_quality_path:
                plot_paths['ortholog_quality'] = ortholog_quality_path
        
        # 5. Network plots
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            network_path = self.network_plotter.create_chromosome_relationship_network(
                results['ortholog_df'], results.get('rearrangement_df'),
                output_path / 'chromosome_network.png'
            )
            if network_path:
                plot_paths['chromosome_network'] = network_path
        
        # 6. PDF Report (if enabled)
        if self.create_pdf_report:
            pdf_path = self._create_pdf_report(results, plot_paths, output_path / 'analysis_report.pdf')
            if pdf_path:
                plot_paths['pdf_report'] = pdf_path
        
        logger.info(f"  Dashboard creation completed: {len(plot_paths)} files generated")
        return plot_paths
    
    def _create_main_dashboard(self, results: Dict[str, Any], output_path: str) -> str:
        """Create main analysis dashboard with key metrics."""
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
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Main dashboard saved: {output_path}")
        return str(output_path)
    
    def _plot_overall_statistics(self, ax, results: Dict[str, Any]):
        """Plot overall analysis statistics."""
        # Extract statistics from results
        categories = []
        values = []
        colors = []
        
        # Ortholog statistics
        if 'ortholog_df' in results and len(results['ortholog_df']) > 0:
            categories.append('Orthologs')
            values.append(len(results['ortholog_df']))
            colors.append('#1f77b4')
        
        # Synteny blocks
        if 'synteny_df' in results:
            categories.append('Synteny\nBlocks')
            values.append(len(results['synteny_df']))
            colors.append('#ff7f0e')
        
        # Inversions
        if 'inversion_df' in results:
            categories.append('Inversions')
            values.append(len(results['inversion_df']))
            colors.append('#d62728')
        
        # Rearrangements
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
            
            # Add quality class annotations
            if 'first_quality' in results:
                ax.text(0, quality_data[0] - 0.1, results['first_quality']['quality_class'],
                       ha='center', va='top', fontsize=10, style='italic')
            if 'second_quality' in results and len(quality_data) > 1:
                ax.text(1, quality_data[1] - 0.1, results['second_quality']['quality_class'],
                       ha='center', va='top', fontsize=10, style='italic')
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
            ax.set_title('Chromosome Mapping Matrix (Gene Counts)', fontweight='bold')
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Number of Genes')
            
            # Add text annotations for significant mappings
            for i in range(len(second_chrs)):
                for j in range(len(first_chrs)):
                    if matrix[i, j] > 0:
                        text = ax.text(j, i, int(matrix[i, j]), ha='center', va='center',
                                     color='white' if matrix[i, j] > matrix.max()/2 else 'black',
                                     fontsize=8, fontweight='bold')
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
                ax.set_title('Synteny Types', fontweight='bold')
                
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
                                  rotation=45, ha='right')
                ax.set_ylabel('Count')
                ax.set_title('Inversion Types', fontweight='bold')
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
            ax.set_title('Alignment Methods', fontweight='bold')
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
    
    def _create_pdf_report(self, results: Dict[str, Any], plot_paths: Dict[str, str],
                          output_path: str) -> str:
        """Create comprehensive PDF report with all plots."""
        try:
            with PdfPages(output_path) as pdf:
                # Page 1: Main dashboard
                if 'main_dashboard' in plot_paths:
                    fig = plt.figure(figsize=(11, 8.5))
                    img = plt.imread(plot_paths['main_dashboard'])
                    plt.imshow(img)
                    plt.axis('off')
                    plt.title('Genome Comparison Analysis Dashboard', fontsize=16, fontweight='bold', pad=20)
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close()
                
                # Page 2: Synteny dotplot
                if 'synteny_dotplot' in plot_paths:
                    fig = plt.figure(figsize=(11, 8.5))
                    img = plt.imread(plot_paths['synteny_dotplot'])
                    plt.imshow(img)
                    plt.axis('off')
                    plt.title('Synteny Dotplot Analysis', fontsize=16, fontweight='bold', pad=20)
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close()
                
                # Page 3: Quality comparison
                if 'quality_comparison' in plot_paths:
                    fig = plt.figure(figsize=(11, 8.5))
                    img = plt.imread(plot_paths['quality_comparison'])
                    plt.imshow(img)
                    plt.axis('off')
                    plt.title('Assembly Quality Comparison', fontsize=16, fontweight='bold', pad=20)
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close()
                
                # Add metadata page
                fig, ax = plt.subplots(figsize=(11, 8.5))
                ax.axis('off')
                
                summary_text = self._generate_summary_text(results)
                ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=12,
                       verticalalignment='top', fontfamily='monospace',
                       bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
                
                ax.set_title('Analysis Summary', fontsize=16, fontweight='bold')
                pdf.savefig(fig, bbox_inches='tight')
                plt.close()
            
            logger.info(f"  PDF report saved: {output_path}")
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Failed to create PDF report: {e}")
            return None
    
    def _generate_summary_text(self, results: Dict[str, Any]) -> str:
        """Generate summary text for PDF report."""
        lines = []
        lines.append("GENOME COMPARISON ANALYSIS SUMMARY")
        lines.append("=" * 50)
        lines.append("")
        
        # Basic statistics
        if 'ortholog_df' in results:
            lines.append(f"Ortholog pairs identified: {len(results['ortholog_df'])}")
            if len(results['ortholog_df']) > 0 and 'similarity' in results['ortholog_df'].columns:
                avg_sim = results['ortholog_df']['similarity'].mean()
                lines.append(f"Average similarity: {avg_sim:.3f}")
        
        if 'synteny_df' in results:
            lines.append(f"Synteny blocks found: {len(results['synteny_df'])}")
        
        if 'inversion_df' in results:
            lines.append(f"Inversions detected: {len(results['inversion_df'])}")
        
        if 'rearrangement_df' in results:
            lines.append(f"Rearrangements found: {len(results['rearrangement_df'])}")
        
        lines.append("")
        
        # Quality assessment
        if 'first_quality' in results and 'second_quality' in results:
            lines.append("Assembly Quality:")
            lines.append(f"  First genome: {results['first_quality']['quality_class']} "
                        f"(score: {results['first_quality']['quality_score']:.3f})")
            lines.append(f"  Second genome: {results['second_quality']['quality_class']} "
                        f"(score: {results['second_quality']['quality_score']:.3f})")
        
        return "\n".join(lines)