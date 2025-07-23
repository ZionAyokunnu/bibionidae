# =============================================================================
# Network Plots (network_plots.py)
# =============================================================================

"""
Network visualization system for complex rearrangement relationships.
Creates network plots showing chromosome relationships and rearrangement patterns.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Tuple
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
from matplotlib.patches import FancyBboxPatch, ConnectionPatch

from ..logger import get_logger

logger = get_logger()

class RearrangementNetworkPlotter:
    """
    Network visualization system for chromosome rearrangements and complex relationships.
    Creates network plots showing chromosome evolution and rearrangement patterns.
    """
    
    def __init__(self, config):
        """
        Initialize network plotter.
        
        Args:
            config: Configuration object with plotting parameters
        """
        self.config = config
        self.figure_size = config.get('network_plot_size', (12, 10))
        self.dpi = config.get('dpi', 300)
        self.node_size_scaling = config.get('node_size_scaling', 1000)
        self.edge_width_scaling = config.get('edge_width_scaling', 5)
        
        # Set style
        plt.style.use('default')
        plt.rcParams['figure.dpi'] = self.dpi
        plt.rcParams['font.size'] = config.get('font_size', 12)
        
        if not NETWORKX_AVAILABLE:
            logger.warning("NetworkX not available - some network plots will be simplified")
        
        logger.info("Network plotter initialized")
    
    def create_chromosome_relationship_network(self, ortholog_df: pd.DataFrame,
                                            rearrangement_df: pd.DataFrame = None,
                                            output_path: str = None) -> str:
        """
        Create network plot showing chromosome relationships.
        
        Args:
            ortholog_df: DataFrame with ortholog pairs
            rearrangement_df: Optional DataFrame with rearrangements
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        if len(ortholog_df) == 0:
            logger.warning("No ortholog data for network plot")
            return None
        
        logger.info("Creating chromosome relationship network...")
        
        if NETWORKX_AVAILABLE:
            return self._create_networkx_plot(ortholog_df, rearrangement_df, output_path)
        else:
            return self._create_simple_network_plot(ortholog_df, rearrangement_df, output_path)
    
    def _create_networkx_plot(self, ortholog_df: pd.DataFrame,
                             rearrangement_df: pd.DataFrame = None,
                             output_path: str = None) -> str:
        """Create network plot using NetworkX."""
        # Create network graph
        G = self._build_chromosome_network(ortholog_df, rearrangement_df)
        
        # Create figure
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        # Generate layout
        pos = self._generate_network_layout(G, ortholog_df)
        
        # Draw network
        self._draw_chromosome_network(ax, G, pos, ortholog_df, rearrangement_df)
        
        # Customize plot
        ax.set_title('Chromosome Relationship Network', fontsize=14, fontweight='bold', pad=20)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Add legend
        self._add_network_legend(ax, rearrangement_df)
        
        # Save plot
        if output_path is None:
            output_path = "chromosome_network.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Chromosome network saved: {output_path}")
        return output_path
    
    def _create_simple_network_plot(self, ortholog_df: pd.DataFrame,
                                   rearrangement_df: pd.DataFrame = None,
                                   output_path: str = None) -> str:
        """Create simplified network plot without NetworkX."""
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        # Get chromosome relationships
        chr_relationships = ortholog_df.groupby(['first_chr', 'second_chr']).agg({
            'busco_id': 'count',
            'similarity': 'mean',
            'confidence': 'mean'
        }).reset_index()
        
        # Get unique chromosomes
        first_chrs = sorted(ortholog_df['first_chr'].unique())
        second_chrs = sorted(ortholog_df['second_chr'].unique())
        
        # Position chromosomes
        y_positions_1 = np.linspace(0.8, 0.2, len(first_chrs))
        y_positions_2 = np.linspace(0.8, 0.2, len(second_chrs))
        
        # Draw first genome chromosomes
        for i, chr_name in enumerate(first_chrs):
            gene_count = len(ortholog_df[ortholog_df['first_chr'] == chr_name])
            size = max(0.02, gene_count / 100 * 0.1)
            
            circle = plt.Circle((0.2, y_positions_1[i]), size, color='#1f77b4', alpha=0.7)
            ax.add_patch(circle)
            ax.text(0.1, y_positions_1[i], chr_name[:10], ha='right', va='center', fontsize=10)
        
        # Draw second genome chromosomes
        for i, chr_name in enumerate(second_chrs):
            gene_count = len(ortholog_df[ortholog_df['second_chr'] == chr_name])
            size = max(0.02, gene_count / 100 * 0.1)
            
            circle = plt.Circle((0.8, y_positions_2[i]), size, color='#ff7f0e', alpha=0.7)
            ax.add_patch(circle)
            ax.text(0.9, y_positions_2[i], chr_name[:10], ha='left', va='center', fontsize=10)
        
        # Draw connections
        for _, row in chr_relationships.iterrows():
        # Draw connections
            for _, row in chr_relationships.iterrows():
                first_idx = first_chrs.index(row['first_chr'])
                second_idx = second_chrs.index(row['second_chr'])
                
                # Line width based on gene count
                width = max(0.5, row['busco_id'] / 20)
                alpha = min(1.0, row['confidence']) if 'confidence' in row else 0.6
                
                ax.plot([0.2, 0.8], [y_positions_1[first_idx], y_positions_2[second_idx]], 
                    'gray', linewidth=width, alpha=alpha)            

        # Customize plot
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title('Chromosome Relationship Network (Simplified)', fontsize=14, fontweight='bold')
        ax.axis('off')
        
        # Add genome labels
        ax.text(0.2, 0.95, 'Genome 1', ha='center', va='bottom', fontsize=14, fontweight='bold', color='#1f77b4')
        ax.text(0.8, 0.95, 'Genome 2', ha='center', va='bottom', fontsize=14, fontweight='bold', color='#ff7f0e')
        
        # Add simple legend
        ax.text(0.5, 0.05, 'Line thickness = number of ortholog genes\nLine transparency = confidence',
               ha='center', va='bottom', fontsize=10, 
               bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        # Save plot
        if output_path is None:
            output_path = "chromosome_network_simple.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Simple chromosome network saved: {output_path}")
        return output_path
    
    def create_rearrangement_flow_diagram(self, rearrangement_df: pd.DataFrame,
                                        output_path: str = None) -> str:
        """
        Create flow diagram showing major rearrangement events.
        
        Args:
            rearrangement_df: DataFrame with rearrangement data
            output_path: Optional output file path
            
        Returns:
            Path to saved plot file
        """
        if len(rearrangement_df) == 0:
            logger.warning("No rearrangement data for flow diagram")
            return None
        
        logger.info("Creating rearrangement flow diagram...")
        
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        # Separate rearrangement types
        splits = rearrangement_df[rearrangement_df['type'] == 'chromosome_split'] if 'type' in rearrangement_df.columns else pd.DataFrame()
        fusions = rearrangement_df[rearrangement_df['type'] == 'chromosome_fusion'] if 'type' in rearrangement_df.columns else pd.DataFrame()
        translocations = rearrangement_df[rearrangement_df['type'] == 'reciprocal_translocation'] if 'type' in rearrangement_df.columns else pd.DataFrame()
        
        # Draw flow elements
        total_events = len(splits) + len(fusions) + len(translocations)
        if total_events == 0:
            ax.text(0.5, 0.5, 'No rearrangement events to display', ha='center', va='center',
                   transform=ax.transAxes, fontsize=14)
            ax.set_title('Rearrangement Flow Diagram')
            ax.axis('off')
        else:
            y_positions = np.linspace(0.8, 0.2, total_events)
            current_y = 0
            
            # Draw splits
            for i, (_, split) in enumerate(splits.iterrows()):
                self._draw_split_event(ax, split, y_positions[current_y])
                current_y += 1
            
            # Draw fusions
            for i, (_, fusion) in enumerate(fusions.iterrows()):
                self._draw_fusion_event(ax, fusion, y_positions[current_y])
                current_y += 1
            
            # Draw translocations
            for i, (_, translocation) in enumerate(translocations.iterrows()):
                self._draw_translocation_event(ax, translocation, y_positions[current_y])
                current_y += 1
        
        # Customize plot
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title('Chromosome Rearrangement Flow Diagram', fontsize=14, fontweight='bold')
        ax.axis('off')
        
        # Add legend
        self._add_flow_legend(ax)
        
        # Save plot
        if output_path is None:
            output_path = "rearrangement_flow.png"
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        
        logger.info(f"  Rearrangement flow diagram saved: {output_path}")
        return output_path
    
    def _build_chromosome_network(self, ortholog_df: pd.DataFrame,
                                rearrangement_df: pd.DataFrame = None) -> 'nx.Graph':
        """Build network graph from chromosome relationships."""
        if not NETWORKX_AVAILABLE:
            return None
            
        G = nx.Graph()
        
        # Add nodes for all chromosomes
        first_chrs = set(ortholog_df['first_chr'].unique())
        second_chrs = set(ortholog_df['second_chr'].unique())
        
        for chr_name in first_chrs:
            G.add_node(f"G1_{chr_name}", genome=1, chromosome=chr_name, type='chromosome')
        
        for chr_name in second_chrs:
            G.add_node(f"G2_{chr_name}", genome=2, chromosome=chr_name, type='chromosome')
        
        # Add edges based on ortholog relationships
        chr_relationships = ortholog_df.groupby(['first_chr', 'second_chr']).agg({
            'busco_id': 'count',
            'similarity': 'mean',
            'confidence': 'mean'
        }).reset_index()
        
        for _, row in chr_relationships.iterrows():
            node1 = f"G1_{row['first_chr']}"
            node2 = f"G2_{row['second_chr']}"
            
            if node1 in G.nodes and node2 in G.nodes:
                G.add_edge(node1, node2,
                          weight=row['busco_id'],
                          similarity=row['similarity'],
                          confidence=row['confidence'],
                          type='synteny')
        
        # Add rearrangement edges if available
        if rearrangement_df is not None and len(rearrangement_df) > 0:
            for _, rearr in rearrangement_df.iterrows():
                self._add_rearrangement_edges(G, rearr)
        
        return G
    
    def _add_rearrangement_edges(self, G: 'nx.Graph', rearrangement: pd.Series):
        """Add rearrangement edges to network graph."""
        if not NETWORKX_AVAILABLE:
            return
            
        rearr_type = rearrangement.get('type', '')
        
        if rearr_type == 'chromosome_split':
            # Add split relationships
            source_chr = rearrangement.get('source_chr', '')
            target_chrs = rearrangement.get('target_chrs', [])
            
            if isinstance(target_chrs, str):
                try:
                    target_chrs = eval(target_chrs) if target_chrs.startswith('[') else [target_chrs]
                except:
                    target_chrs = [target_chrs]
            
            source_node = f"G1_{source_chr}"
            for target_chr in target_chrs:
                target_node = f"G2_{target_chr}"
                if source_node in G.nodes and target_node in G.nodes:
                    G.add_edge(source_node, target_node, type='split',
                              confidence=rearrangement.get('confidence', 0.5))
        
        elif rearr_type == 'chromosome_fusion':
            # Add fusion relationships
            target_chr = rearrangement.get('target_chr', '')
            source_chrs = rearrangement.get('source_chrs', [])
            
            if isinstance(source_chrs, str):
                try:
                    source_chrs = eval(source_chrs) if source_chrs.startswith('[') else [source_chrs]
                except:
                    source_chrs = [source_chrs]
            
            target_node = f"G2_{target_chr}"
            for source_chr in source_chrs:
                source_node = f"G1_{source_chr}"
                if source_node in G.nodes and target_node in G.nodes:
                    G.add_edge(source_node, target_node, type='fusion',
                              confidence=rearrangement.get('confidence', 0.5))
        
        elif rearr_type == 'reciprocal_translocation':
            # Add translocation relationships
            chr1 = rearrangement.get('chromosome_1', '')
            chr2 = rearrangement.get('chromosome_2', '')
            
            chr1_g1 = f"G1_{chr1}"
            chr1_g2 = f"G2_{chr1}"
            chr2_g1 = f"G1_{chr2}"
            chr2_g2 = f"G2_{chr2}"
            
            # Add cross-connections
            if all(node in G.nodes for node in [chr1_g1, chr2_g2, chr2_g1, chr1_g2]):
                G.add_edge(chr1_g1, chr2_g2, type='translocation',
                          confidence=rearrangement.get('confidence', 0.5))
                G.add_edge(chr2_g1, chr1_g2, type='translocation',
                          confidence=rearrangement.get('confidence', 0.5))
    
    def _generate_network_layout(self, G: 'nx.Graph', ortholog_df: pd.DataFrame) -> Dict:
        """Generate layout positions for network nodes."""
        if not NETWORKX_AVAILABLE or G is None:
            return {}
            
        # Separate genomes
        genome1_nodes = [n for n in G.nodes if G.nodes[n]['genome'] == 1]
        genome2_nodes = [n for n in G.nodes if G.nodes[n]['genome'] == 2]
        
        # Position genome 1 nodes on the left
        pos = {}
        
        # Sort chromosomes by size (number of genes)
        chr_sizes = ortholog_df.groupby('first_chr')['busco_id'].count().sort_values(ascending=False)
        sorted_g1_chrs = [f"G1_{chr_name}" for chr_name in chr_sizes.index if f"G1_{chr_name}" in genome1_nodes]
        
        chr_sizes = ortholog_df.groupby('second_chr')['busco_id'].count().sort_values(ascending=False)
        sorted_g2_chrs = [f"G2_{chr_name}" for chr_name in chr_sizes.index if f"G2_{chr_name}" in genome2_nodes]
        
        # Position nodes
        for i, node in enumerate(sorted_g1_chrs):
            pos[node] = (-1, len(sorted_g1_chrs)/2 - i)
        
        for i, node in enumerate(sorted_g2_chrs):
            pos[node] = (1, len(sorted_g2_chrs)/2 - i)
        
        return pos
    
    def _draw_chromosome_network(self, ax, G: 'nx.Graph', pos: Dict,
                               ortholog_df: pd.DataFrame, rearrangement_df: pd.DataFrame = None):
        """Draw the chromosome network with appropriate styling."""
        if not NETWORKX_AVAILABLE or G is None:
            return
            
        # Calculate node sizes based on gene count
        node_sizes = []
        node_colors = []
        
        for node in G.nodes:
            genome = G.nodes[node]['genome']
            chr_name = G.nodes[node]['chromosome']
            
            if genome == 1:
                gene_count = len(ortholog_df[ortholog_df['first_chr'] == chr_name])
                color = '#1f77b4'  # Blue for genome 1
            else:
                gene_count = len(ortholog_df[ortholog_df['second_chr'] == chr_name])
                color = '#ff7f0e'  # Orange for genome 2
            
            node_sizes.append(max(100, gene_count * self.node_size_scaling / 10))
            node_colors.append(color)
        
        # Draw nodes
        nx.draw_networkx_nodes(G, pos, ax=ax,
                              node_size=node_sizes,
                              node_color=node_colors,
                              alpha=0.7,
                              edgecolors='black',
                              linewidths=1)
        
        # Draw edges with different styles for different types
        synteny_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'synteny']
        rearr_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') != 'synteny']
        
        # Draw synteny edges
        if synteny_edges:
            edge_weights = [G[u][v]['weight'] for u, v in synteny_edges]
            max_weight = max(edge_weights) if edge_weights else 1
            edge_widths = [w / max_weight * self.edge_width_scaling for w in edge_weights]
            
            nx.draw_networkx_edges(G, pos, ax=ax,
                                  edgelist=synteny_edges,
                                  width=edge_widths,
                                  edge_color='gray',
                                  alpha=0.6)
        
        # Draw rearrangement edges
        if rearr_edges:
            nx.draw_networkx_edges(G, pos, ax=ax,
                                  edgelist=rearr_edges,
                                  width=2,
                                  edge_color='red',
                                  style='dashed',
                                  alpha=0.8)
        
        # Draw labels
        labels = {node: G.nodes[node]['chromosome'] for node in G.nodes}
        nx.draw_networkx_labels(G, pos, labels, ax=ax,
                               font_size=10, font_weight='bold')
        
        # Add genome labels
        if pos:
            max_y = max(pos.values(), key=lambda x: x[1])[1] if pos else 0.5
            ax.text(-1, max_y + 0.5, 'Genome 1',
                   ha='center', va='bottom', fontsize=14, fontweight='bold', color='#1f77b4')
            ax.text(1, max_y + 0.5, 'Genome 2',
                   ha='center', va='bottom', fontsize=14, fontweight='bold', color='#ff7f0e')
    
    def _add_network_legend(self, ax, rearrangement_df: pd.DataFrame = None):
        """Add legend to network plot."""
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4',
                      markersize=10, label='Genome 1 Chromosomes'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e',
                      markersize=10, label='Genome 2 Chromosomes'),
            plt.Line2D([0], [0], color='gray', linewidth=3, label='Synteny Relationships')
        ]
        
        if rearrangement_df is not None and len(rearrangement_df) > 0:
            legend_elements.append(
                plt.Line2D([0], [0], color='red', linewidth=2, linestyle='--',
                          label='Rearrangements')
            )
        
        ax.legend(handles=legend_elements, loc='upper right', frameon=True,
                 fancybox=True, shadow=True)
    
    def _draw_split_event(self, ax, split: pd.Series, y_pos: float):
        """Draw chromosome split event."""
        # Source chromosome
        source_chr = split.get('source_chr', 'Unknown')
        target_chrs = split.get('target_chrs', [])
        total_genes = split.get('total_genes', 0)
        
        # Handle target_chrs format
        if isinstance(target_chrs, str):
            try:
                target_chrs = eval(target_chrs) if target_chrs.startswith('[') else [target_chrs]
            except:
                target_chrs = [target_chrs]
        
        source_box = FancyBboxPatch((0.1, y_pos-0.02), 0.15, 0.04,
                                   boxstyle="round,pad=0.01", 
                                   facecolor='#1f77b4', alpha=0.7)
        ax.add_patch(source_box)
        ax.text(0.175, y_pos, str(source_chr)[:8], ha='center', va='center', 
               fontweight='bold', fontsize=8)
        
        # Arrow
        ax.annotate('', xy=(0.35, y_pos), xytext=(0.25, y_pos),
                   arrowprops=dict(arrowstyle='->', lw=2, color='black'))
        
        # Target chromosomes
        n_targets = len(target_chrs)
        for i, target_chr in enumerate(target_chrs):
            target_y = y_pos + (i - (n_targets-1)/2) * 0.03
            target_box = FancyBboxPatch((0.45, target_y-0.015), 0.12, 0.03,
                                       boxstyle="round,pad=0.01",
                                       facecolor='#ff7f0e', alpha=0.7)
            ax.add_patch(target_box)
            ax.text(0.51, target_y, str(target_chr)[:8], ha='center', va='center',
                   fontweight='bold', fontsize=7)
        
        # Label
        ax.text(0.7, y_pos, f'Split ({total_genes} genes)',
               ha='left', va='center', fontsize=9, style='italic')
    
    def _draw_fusion_event(self, ax, fusion: pd.Series, y_pos: float):
        """Draw chromosome fusion event."""
        # Get fusion data
        source_chrs = fusion.get('source_chrs', [])
        target_chr = fusion.get('target_chr', 'Unknown')
        total_genes = fusion.get('total_genes', 0)
        
        # Handle source_chrs format
        if isinstance(source_chrs, str):
            try:
                source_chrs = eval(source_chrs) if source_chrs.startswith('[') else [source_chrs]
            except:
                source_chrs = [source_chrs]
        
        # Source chromosomes
        n_sources = len(source_chrs)
        for i, source_chr in enumerate(source_chrs):
            source_y = y_pos + (i - (n_sources-1)/2) * 0.03
            source_box = FancyBboxPatch((0.1, source_y-0.015), 0.12, 0.03,
                                       boxstyle="round,pad=0.01",
                                       facecolor='#1f77b4', alpha=0.7)
            ax.add_patch(source_box)
            ax.text(0.16, source_y, str(source_chr)[:8], ha='center', va='center',
                   fontweight='bold', fontsize=7)
        
        # Arrow
        ax.annotate('', xy=(0.35, y_pos), xytext=(0.25, y_pos),
                   arrowprops=dict(arrowstyle='->', lw=2, color='black'))
        
        # Target chromosome
        target_box = FancyBboxPatch((0.45, y_pos-0.02), 0.15, 0.04,
                                   boxstyle="round,pad=0.01",
                                   facecolor='#ff7f0e', alpha=0.7)
        ax.add_patch(target_box)
        ax.text(0.525, y_pos, str(target_chr)[:8], ha='center', va='center',
               fontweight='bold', fontsize=8)
        
        # Label
        ax.text(0.7, y_pos, f'Fusion ({total_genes} genes)',
               ha='left', va='center', fontsize=9, style='italic')
    
    def _draw_translocation_event(self, ax, translocation: pd.Series, y_pos: float):
        """Draw reciprocal translocation event."""
        chr1 = translocation.get('chromosome_1', 'Chr1')
        chr2 = translocation.get('chromosome_2', 'Chr2')
        total_genes = translocation.get('total_genes', 0)
        
        # Chromosome 1
        chr1_box = FancyBboxPatch((0.1, y_pos+0.01), 0.12, 0.03,
                                 boxstyle="round,pad=0.01",
                                 facecolor='#2ca02c', alpha=0.7)
        ax.add_patch(chr1_box)
        ax.text(0.16, y_pos+0.025, str(chr1)[:8], 
               ha='center', va='center', fontweight='bold', fontsize=7)
        
        # Chromosome 2
        chr2_box = FancyBboxPatch((0.1, y_pos-0.025), 0.12, 0.03,
                                 boxstyle="round,pad=0.01",
                                 facecolor='#2ca02c', alpha=0.7)
        ax.add_patch(chr2_box)
        ax.text(0.16, y_pos-0.01, str(chr2)[:8], 
               ha='center', va='center', fontweight='bold', fontsize=7)
        
        # Reciprocal arrows
        ax.annotate('', xy=(0.35, y_pos+0.025), xytext=(0.25, y_pos+0.025),
                   arrowprops=dict(arrowstyle='<->', lw=2, color='red'))
        ax.annotate('', xy=(0.35, y_pos-0.01), xytext=(0.25, y_pos-0.01),
                   arrowprops=dict(arrowstyle='<->', lw=2, color='red'))
        
        # Crossed connection
        ax.plot([0.3, 0.3], [y_pos+0.025, y_pos-0.01], 'r--', lw=1, alpha=0.7)
        
        # Label
        ax.text(0.45, y_pos, f'Reciprocal Translocation\n({total_genes} genes)',
               ha='left', va='center', fontsize=9, style='italic')
    
    def _add_flow_legend(self, ax):
        """Add legend to flow diagram."""
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor='#1f77b4', alpha=0.7, label='Source Chromosomes'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#ff7f0e', alpha=0.7, label='Target Chromosomes'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#2ca02c', alpha=0.7, label='Translocation Chromosomes'),
            plt.Line2D([0], [0], color='black', linewidth=2, label='Split/Fusion'),
            plt.Line2D([0], [0], color='red', linewidth=2, linestyle='--', label='Translocation')
        ]
        
        ax.legend(handles=legend_elements, loc='upper right', frameon=True,
                 fancybox=True, shadow=True)