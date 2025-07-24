"""
SyRI Integration for Genome Inversion Analyzer
Creates SyRI-compatible outputs and SyRI-style visualizations
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
import logging

logger = logging.getLogger(__name__)


class SyRIIntegrator:
    """
    Integration with SyRI for synteny visualization and analysis
    Creates SyRI-compatible formats and SyRI-style plots
    """
    
    def __init__(self, output_dir: Union[str, Path]):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up SyRI-compatible color scheme
        self.syri_colors = {
            'syn': '#1f77b4',      # Blue for syntenic regions
            'inv': '#ff7f0e',      # Orange for inversions
            'trans': '#2ca02c',    # Green for translocations
            'dup': '#d62728',      # Red for duplications
            'del': '#9467bd',      # Purple for deletions
            'ins': '#8c564b',      # Brown for insertions
            'snp': '#e377c2',      # Pink for SNPs
            'indel': '#7f7f7f'     # Gray for indels
        }
    
    def create_syri_compatible_output(self, 
                                    ortholog_df: pd.DataFrame,
                                    synteny_df: pd.DataFrame,
                                    inversion_df: pd.DataFrame,
                                    species1: str,
                                    species2: str) -> Path:
        """
        Create SyRI-compatible output file from inversion analysis results
        
        Args:
            ortholog_df: Ortholog pairs DataFrame
            synteny_df: Synteny blocks DataFrame
            inversion_df: Inversion events DataFrame
            species1: First species name
            species2: Second species name
            
        Returns:
            Path to SyRI-compatible output file
        """
        logger.info(f"Creating SyRI-compatible output for {species1} vs {species2}")
        
        syri_records = []
        
        # Process syntenic regions
        for _, synteny in synteny_df.iterrows():
            # Get orthologs in this synteny block
            block_orthologs = ortholog_df[
                (ortholog_df['first_chr'] == synteny.get('first_chr', '')) &
                (ortholog_df['second_chr'] == synteny.get('second_chr', ''))
            ]
            
            if len(block_orthologs) > 0:
                # Define synteny block boundaries
                first_start = block_orthologs['first_start'].min()
                first_end = block_orthologs['first_end'].max()
                second_start = block_orthologs['second_start'].min()
                second_end = block_orthologs['second_end'].max()
                
                # Determine if block is inverted
                strand_consistency = block_orthologs['first_strand'] == block_orthologs['second_strand']
                is_inverted = strand_consistency.mean() < 0.5
                
                syri_type = 'INVTR' if is_inverted else 'SYN'
                
                syri_record = {
                    'chr1': synteny.get('first_chr', ''),
                    'start1': int(first_start),
                    'end1': int(first_end),
                    'chr2': synteny.get('second_chr', ''),
                    'start2': int(second_start),
                    'end2': int(second_end),
                    'type': syri_type,
                    'confidence': synteny.get('confidence', 1.0),
                    'size': first_end - first_start
                }
                syri_records.append(syri_record)
        
        # Process specific inversions
        for _, inversion in inversion_df.iterrows():
            syri_record = {
                'chr1': inversion.get('first_chr', ''),
                'start1': inversion.get('first_start', 0),
                'end1': inversion.get('first_end', 0),
                'chr2': inversion.get('second_chr', ''),
                'start2': inversion.get('second_start', 0),
                'end2': inversion.get('second_end', 0),
                'type': 'INV',
                'confidence': inversion.get('confidence', 1.0),
                'size': inversion.get('size_genes', 0)
            }
            syri_records.append(syri_record)
        
        # Create SyRI-compatible DataFrame
        syri_df = pd.DataFrame(syri_records)
        
        # Save in SyRI format
        output_file = self.output_dir / f'syri_output_{species1}_vs_{species2}.tsv'
        
        # SyRI format: chr1 start1 end1 chr2 start2 end2 type confidence
        if not syri_df.empty:
            syri_df[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'confidence']].to_csv(
                output_file, sep='\t', index=False, header=False
            )
        
        logger.info(f"Created SyRI-compatible output: {output_file}")
        return output_file
    
    def create_syri_style_plot(self,
                             ortholog_df: pd.DataFrame,
                             synteny_df: pd.DataFrame,
                             inversion_df: pd.DataFrame,
                             species1: str,
                             species2: str,
                             chromosome_lengths: Dict[str, int] = None) -> Path:
        """
        Create SyRI-style synteny plot
        
        Args:
            ortholog_df: Ortholog pairs DataFrame
            synteny_df: Synteny blocks DataFrame 
            inversion_df: Inversion events DataFrame
            species1: First species name
            species2: Second species name
            chromosome_lengths: Optional chromosome length information
            
        Returns:
            Path to generated plot
        """
        logger.info(f"Creating SyRI-style plot for {species1} vs {species2}")
        
        # Set up the plot
        fig, ax = plt.subplots(figsize=(16, 10))
        
        # Get unique chromosomes
        chr1_list = sorted(ortholog_df['first_chr'].unique())
        chr2_list = sorted(ortholog_df['second_chr'].unique())
        
        # Create chromosome position mappings
        chr1_positions = self._create_chromosome_positions(chr1_list, chromosome_lengths, 'first')
        chr2_positions = self._create_chromosome_positions(chr2_list, chromosome_lengths, 'second')
        
        # Plot chromosome backbones
        self._plot_chromosome_backbones(ax, chr1_positions, chr2_positions, species1, species2)
        
        # Plot synteny blocks
        self._plot_synteny_blocks(ax, synteny_df, ortholog_df, chr1_positions, chr2_positions)
        
        # Plot inversions
        self._plot_inversions(ax, inversion_df, chr1_positions, chr2_positions)
        
        # Plot individual orthologs as connecting lines
        self._plot_ortholog_connections(ax, ortholog_df, chr1_positions, chr2_positions)
        
        # Customize plot
        ax.set_xlim(-0.5, max(max(chr1_positions.values()), max(chr2_positions.values())) + 0.5)
        ax.set_ylim(-1, len(chr1_list) + len(chr2_list) + 1)
        
        ax.set_xlabel('Genomic Position (Mb)', fontsize=12)
        ax.set_ylabel('Chromosomes', fontsize=12)
        ax.set_title(f'Synteny Analysis: {species1} vs {species2}', fontsize=14, fontweight='bold')
        
        # Add legend
        self._add_syri_legend(ax)
        
        # Save plot
        plot_file = self.output_dir / f'syri_style_plot_{species1}_vs_{species2}.png'
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Created SyRI-style plot: {plot_file}")
        return plot_file
    
    def create_circular_synteny_plot(self,
                                   ortholog_df: pd.DataFrame,
                                   inversion_df: pd.DataFrame,
                                   species1: str,
                                   species2: str) -> Path:
        """
        Create circular synteny plot (Circos-style)
        
        Args:
            ortholog_df: Ortholog pairs DataFrame
            inversion_df: Inversion events DataFrame
            species1: First species name
            species2: Second species name
            
        Returns:
            Path to generated circular plot
        """
        logger.info(f"Creating circular synteny plot for {species1} vs {species2}")
        
        fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))
        
        # Get chromosome lists
        chr1_list = sorted(ortholog_df['first_chr'].unique())
        chr2_list = sorted(ortholog_df['second_chr'].unique())
        
        # Create angular positions for chromosomes
        n_chr1 = len(chr1_list)
        n_chr2 = len(chr2_list)
        
        # First species on top half, second species on bottom half
        chr1_angles = np.linspace(0, np.pi, n_chr1)
        chr2_angles = np.linspace(np.pi, 2*np.pi, n_chr2)
        
        # Plot chromosome arcs
        radius = 1.0
        for i, chr_name in enumerate(chr1_list):
            ax.plot([chr1_angles[i], chr1_angles[i]], [0.8, radius], 'k-', linewidth=3)
            ax.text(chr1_angles[i], radius + 0.1, chr_name, ha='center', va='center')
        
        for i, chr_name in enumerate(chr2_list):
            ax.plot([chr2_angles[i], chr2_angles[i]], [0.8, radius], 'k-', linewidth=3)
            ax.text(chr2_angles[i], radius + 0.1, chr_name, ha='center', va='center')
        
        # Plot synteny connections
        for _, ortholog in ortholog_df.iterrows():
            chr1_idx = chr1_list.index(ortholog['first_chr'])
            chr2_idx = chr2_list.index(ortholog['second_chr'])
            
            angle1 = chr1_angles[chr1_idx]
            angle2 = chr2_angles[chr2_idx]
            
            # Determine color based on strand consistency
            same_strand = ortholog['first_strand'] == ortholog['second_strand']
            color = self.syri_colors['syn'] if same_strand else self.syri_colors['inv']
            alpha = 0.1
            
            # Draw connection arc
            angles = np.linspace(angle1, angle2, 100)
            radii = 0.7 * np.sin(np.pi * (angles - angle1) / (angle2 - angle1))
            
            ax.plot(angles, radii, color=color, alpha=alpha, linewidth=0.5)
        
        # Highlight inversions
        for _, inversion in inversion_df.iterrows():
            try:
                chr1_idx = chr1_list.index(inversion['first_chr'])
                chr2_idx = chr2_list.index(inversion['second_chr'])
                
                angle1 = chr1_angles[chr1_idx]
                angle2 = chr2_angles[chr2_idx]
                
                # Draw thick inversion connection
                angles = np.linspace(angle1, angle2, 100)
                radii = 0.6 * np.sin(np.pi * (angles - angle1) / (angle2 - angle1))
                
                ax.plot(angles, radii, color=self.syri_colors['inv'], linewidth=2, alpha=0.8)
            except (ValueError, KeyError):
                continue
        
        # Customize circular plot
        ax.set_ylim(0, 1.2)
        ax.set_title(f'Circular Synteny: {species1} vs {species2}', fontsize=14, fontweight='bold', pad=20)
        ax.grid(False)
        ax.set_rticks([])
        ax.set_thetagrids([])
        
        # Save plot
        plot_file = self.output_dir / f'circular_synteny_{species1}_vs_{species2}.png'
        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Created circular synteny plot: {plot_file}")
        return plot_file
    
    def create_chromosome_comparison_plot(self,
                                        ortholog_df: pd.DataFrame,
                                        inversion_df: pd.DataFrame,
                                        species1: str,
                                        species2: str) -> Path:
        """
        Create detailed chromosome-by-chromosome comparison plot
        
        Args:
            ortholog_df: Ortholog pairs DataFrame
            inversion_df: Inversion events DataFrame
            species1: First species name
            species2: Second species name
            
        Returns:
            Path to generated plot
        """
        logger.info(f"Creating chromosome comparison plot for {species1} vs {species2}")
        
        # Get unique chromosome pairs
        chr_pairs = ortholog_df.groupby(['first_chr', 'second_chr']).size().reset_index(name='count')
        chr_pairs = chr_pairs.sort_values('count', ascending=False)
        
        n_pairs = min(len(chr_pairs), 9)  # Show top 9 pairs
        
        fig, axes = plt.subplots(3, 3, figsize=(15, 15))
        axes = axes.flatten()
        
        for i, (_, pair) in enumerate(chr_pairs.head(n_pairs).iterrows()):
            ax = axes[i]
            
            chr1 = pair['first_chr']
            chr2 = pair['second_chr']
            
            # Get orthologs for this chromosome pair
            pair_orthologs = ortholog_df[
                (ortholog_df['first_chr'] == chr1) & 
                (ortholog_df['second_chr'] == chr2)
            ]
            
            # Get inversions for this pair
            pair_inversions = inversion_df[
                (inversion_df['first_chr'] == chr1) & 
                (inversion_df['second_chr'] == chr2)
            ]
            
            # Create dot plot for this chromosome pair
            self._create_chromosome_dotplot(ax, pair_orthologs, pair_inversions, chr1, chr2)
            
            ax.set_title(f'{chr1} vs {chr2}\n({len(pair_orthologs)} orthologs)', fontsize=10)
        
        # Hide unused subplots
        for i in range(n_pairs, len(axes)):
            axes[i].set_visible(False)
        
        plt.suptitle(f'Chromosome Comparisons: {species1} vs {species2}', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        # Save plot
        plot_file = self.output_dir / f'chromosome_comparison_{species1}_vs_{species2}.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Created chromosome comparison plot: {plot_file}")
        return plot_file
    
    # Helper methods
    def _create_chromosome_positions(self, chr_list: List[str], 
                                   chromosome_lengths: Dict[str, int] = None,
                                   genome: str = 'first') -> Dict[str, float]:
        """Create cumulative positions for chromosomes"""
        positions = {}
        current_pos = 0
        
        for chr_name in chr_list:
            positions[chr_name] = current_pos
            
            # Use actual length if available, otherwise estimate
            if chromosome_lengths and chr_name in chromosome_lengths:
                chr_length = chromosome_lengths[chr_name] / 1_000_000  # Convert to Mb
            else:
                chr_length = 10  # Default 10 Mb spacing
            
            current_pos += chr_length + 2  # Add 2 Mb gap between chromosomes
        
        return positions
    
    def _plot_chromosome_backbones(self, ax, chr1_positions: Dict, chr2_positions: Dict,
                                 species1: str, species2: str):
        """Plot chromosome backbone representations"""
        y1_base = len(chr1_positions) + 2
        y2_base = 1
        
        # Plot species 1 chromosomes (top)
        for chr_name, x_pos in chr1_positions.items():
            ax.plot([x_pos, x_pos + 8], [y1_base, y1_base], 'k-', linewidth=4)
            ax.text(x_pos + 4, y1_base + 0.3, chr_name, ha='center', va='bottom', fontsize=8)
        
        # Plot species 2 chromosomes (bottom)
        for chr_name, x_pos in chr2_positions.items():
            ax.plot([x_pos, x_pos + 8], [y2_base, y2_base], 'k-', linewidth=4)
            ax.text(x_pos + 4, y2_base - 0.3, chr_name, ha='center', va='top', fontsize=8)
        
        # Add species labels
        ax.text(-2, y1_base, species1, ha='right', va='center', fontsize=12, fontweight='bold')
        ax.text(-2, y2_base, species2, ha='right', va='center', fontsize=12, fontweight='bold')
    
    def _plot_synteny_blocks(self, ax, synteny_df: pd.DataFrame, ortholog_df: pd.DataFrame,
                        chr1_positions: Dict, chr2_positions: Dict):
        """Plot synteny blocks as colored regions"""
        y1_base = len(chr1_positions) + 2
        y2_base = 1
        
        for _, synteny in synteny_df.iterrows():
            chr1 = synteny.get('first_chr', '')
            chr2 = synteny.get('second_chr', '')
            
            if chr1 in chr1_positions and chr2 in chr2_positions:
                # Get block boundaries from orthologs
                block_orthologs = ortholog_df[
                    (ortholog_df['first_chr'] == chr1) & 
                    (ortholog_df['second_chr'] == chr2)
                ]
                
                if len(block_orthologs) > 0:
                    x1_start = chr1_positions[chr1] + block_orthologs['first_start'].min() / 1_000_000
                    x1_end = chr1_positions[chr1] + block_orthologs['first_end'].max() / 1_000_000
                    x2_start = chr2_positions[chr2] + block_orthologs['second_start'].min() / 1_000_000
                    x2_end = chr2_positions[chr2] + block_orthologs['second_end'].max() / 1_000_000
                    
                    # Determine block color based on synteny type
                    synteny_type = synteny.get('synteny_type', 'colinear')
                    color = self.syri_colors['inv'] if synteny_type == 'inverted' else self.syri_colors['syn']
                    
                    # Draw synteny block regions
                    import matplotlib.patches as patches
                    rect1 = patches.Rectangle((x1_start, y1_base - 0.1), x1_end - x1_start, 0.2, 
                                            facecolor=color, alpha=0.6, edgecolor='none')
                    rect2 = patches.Rectangle((x2_start, y2_base - 0.1), x2_end - x2_start, 0.2,
                                            facecolor=color, alpha=0.6, edgecolor='none')
                    
                    ax.add_patch(rect1)
                    ax.add_patch(rect2)

    def integrate_with_syri(output_dir: Union[str, Path],
                            pairwise_results: Dict[str, Any],
                            species_pairs: List[Tuple[str, str]] = None) -> Dict[str, Any]:
        """
        Complete SyRI integration for multi-species analysis
        
        Args:
            output_dir: Output directory for SyRI files
            pairwise_results: Dictionary of pairwise analysis results
            species_pairs: Optional list of specific pairs to process
            
        Returns:
            Dictionary with SyRI integration results
        """
        logger.info("Starting complete SyRI integration")
        
        # Initialize SyRI integrator
        syri_dir = Path(output_dir) / 'syri_integration'
        syri_integrator = SyRIIntegrator(syri_dir)
        
        # Process all species pairs or specified ones
        if species_pairs is None:
            species_pairs = []
            for pair_key in pairwise_results.keys():
                if '_vs_' in pair_key and 'error' not in pairwise_results[pair_key]:
                    species1, species2 = pair_key.split('_vs_')
                    species_pairs.append((species1, species2))
        
        syri_results = {}
        
        for species1, species2 in species_pairs:
            try:
                logger.info(f"Processing SyRI integration: {species1} vs {species2}")
                
                pair_key = f"{species1}_vs_{species2}"
                pair_data = pairwise_results.get(pair_key, {})
                
                if 'error' in pair_data:
                    logger.warning(f"Skipping {pair_key} due to analysis error")
                    continue
                
                # Extract data
                ortholog_df = pair_data.get('ortholog_df', pd.DataFrame())
                synteny_df = pair_data.get('synteny_df', pd.DataFrame())
                inversion_df = pair_data.get('inversion_df', pd.DataFrame())
                
                # Create SyRI outputs
                syri_output = syri_integrator.create_syri_compatible_output(
                    ortholog_df, synteny_df, inversion_df, species1, species2
                )
                
                # Create visualizations
                syri_plot = syri_integrator.create_syri_style_plot(
                    ortholog_df, synteny_df, inversion_df, species1, species2
                )
                
                circular_plot = syri_integrator.create_circular_synteny_plot(
                    ortholog_df, inversion_df, species1, species2
                )
                
                chromosome_plot = syri_integrator.create_chromosome_comparison_plot(
                    ortholog_df, inversion_df, species1, species2
                )
                
                syri_results[pair_key] = {
                    'syri_output_file': str(syri_output),
                    'syri_style_plot': str(syri_plot),
                    'circular_plot': str(circular_plot),
                    'chromosome_plot': str(chromosome_plot),
                    'species_pair': (species1, species2)
                }
                
                logger.info(f"Completed SyRI integration: {pair_key}")
                
            except Exception as e:
                logger.error(f"SyRI integration failed for {species1} vs {species2}: {e}")
                syri_results[f"{species1}_vs_{species2}"] = {'error': str(e)}
        
        # Create summary report
        summary_report = syri_integrator.create_summary_report(species_pairs, pairwise_results)
        syri_results['summary_report'] = str(summary_report)
        
        logger.info(f"SyRI integration completed for {len(species_pairs)} species pairs")
        return syri_results

    def _plot_inversions(self, ax, inversion_df: pd.DataFrame, 
                        chr1_positions: Dict, chr2_positions: Dict):
        """Plot specific inversion events"""
        y1_base = len(chr1_positions) + 2
        y2_base = 1
        
        for _, inversion in inversion_df.iterrows():
            chr1 = inversion.get('first_chr', '')
            chr2 = inversion.get('second_chr', '')
            
            if chr1 in chr1_positions and chr2 in chr2_positions:
                x1 = chr1_positions[chr1] + inversion.get('first_start', 0) / 1_000_000
                x2 = chr2_positions[chr2] + inversion.get('second_start', 0) / 1_000_000
                
                # Draw inversion connection with distinctive style
                ax.plot([x1, x2], [y1_base, y2_base], 
                       color=self.syri_colors['inv'], linewidth=2, alpha=0.8, linestyle='--')
                
                # Add inversion markers
                ax.scatter([x1], [y1_base], color=self.syri_colors['inv'], s=50, marker='v', zorder=5)
                ax.scatter([x2], [y2_base], color=self.syri_colors['inv'], s=50, marker='^', zorder=5)
    
    def _plot_ortholog_connections(self, ax, ortholog_df: pd.DataFrame,
                                 chr1_positions: Dict, chr2_positions: Dict):
        """Plot individual ortholog connections"""
        from matplotlib.collections import LineCollection
        
        y1_base = len(chr1_positions) + 2
        y2_base = 1
        
        # Sample orthologs to avoid overcrowding
        if len(ortholog_df) > 1000:
            sample_orthologs = ortholog_df.sample(n=1000, random_state=42)
        else:
            sample_orthologs = ortholog_df
        
        lines = []
        colors = []
        
        for _, ortholog in sample_orthologs.iterrows():
            chr1 = ortholog['first_chr']
            chr2 = ortholog['second_chr']
            
            if chr1 in chr1_positions and chr2 in chr2_positions:
                x1 = chr1_positions[chr1] + ortholog['first_start'] / 1_000_000
                x2 = chr2_positions[chr2] + ortholog['second_start'] / 1_000_000
                
                lines.append([(x1, y1_base), (x2, y2_base)])
                
                # Color by strand consistency
                same_strand = ortholog['first_strand'] == ortholog['second_strand']
                color = self.syri_colors['syn'] if same_strand else self.syri_colors['inv']
                colors.append(color)
        
        # Plot all connections at once for efficiency
        lc = LineCollection(lines, colors=colors, alpha=0.3, linewidths=0.5)
        ax.add_collection(lc)
    
    def _add_syri_legend(self, ax):
        """Add SyRI-style legend to plot"""
        import matplotlib.pyplot as plt
        
        legend_elements = [
            plt.Line2D([0], [0], color=self.syri_colors['syn'], lw=2, label='Syntenic'),
            plt.Line2D([0], [0], color=self.syri_colors['inv'], lw=2, label='Inverted'),
            plt.Line2D([0], [0], color=self.syri_colors['inv'], lw=2, linestyle='--', label='Inversion Event')
        ]
        
        ax.legend(handles=legend_elements, loc='upper right', frameon=True, 
                 fancybox=True, shadow=True)
    
    def _create_chromosome_dotplot(self, ax, orthologs: pd.DataFrame, inversions: pd.DataFrame,
                                 chr1: str, chr2: str):
        """Create dot plot for chromosome pair comparison"""
        if len(orthologs) == 0:
            ax.text(0.5, 0.5, 'No orthologs', ha='center', va='center', transform=ax.transAxes)
            return
        
        # Plot orthologs as dots
        same_strand = orthologs['first_strand'] == orthologs['second_strand']
        
        # Syntenic orthologs
        syn_orthologs = orthologs[same_strand]
        if len(syn_orthologs) > 0:
            ax.scatter(syn_orthologs['first_start'], syn_orthologs['second_start'],
                      c=self.syri_colors['syn'], alpha=0.6, s=10, label='Syntenic')
        
        # Inverted orthologs
        inv_orthologs = orthologs[~same_strand]
        if len(inv_orthologs) > 0:
            ax.scatter(inv_orthologs['first_start'], inv_orthologs['second_start'],
                      c=self.syri_colors['inv'], alpha=0.6, s=10, label='Inverted')
        
        # Highlight specific inversions
        if len(inversions) > 0:
            ax.scatter(inversions['first_start'], inversions['second_start'],
                      c=self.syri_colors['inv'], s=50, marker='D', 
                      edgecolors='black', linewidth=0.5, label='Inversion Event')
        
        ax.set_xlabel(f'{chr1} position')
        ax.set_ylabel(f'{chr2} position')
        
        # Add diagonal line for perfect synteny
        if len(orthologs) > 0:
            min_pos = min(orthologs['first_start'].min(), orthologs['second_start'].min())
            max_pos = max(orthologs['first_end'].max(), orthologs['second_end'].max())
            ax.plot([min_pos, max_pos], [min_pos, max_pos], 'k--', alpha=0.3, linewidth=1)

        # Create summary report
        summary_report = syri_integrator.create_summary_report(species_pairs, pairwise_results)
        syri_results['summary_report'] = str(summary_report)

        logger.info(f"SyRI integration completed for {len(species_pairs)} species pairs")
        return syri_results