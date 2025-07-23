"""
Visualization module for the Genome Inversion Analyzer
Creates comprehensive plots for synteny, inversions, and quality assessment
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)


def create_enhanced_visualizations(output_dir, results_dict, config):
    """
    Create comprehensive visualizations including dotplots, synteny maps, and summary plots
    """    
    plots_dir = os.path.join(output_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Set style and parameters
    plt.style.use('default')
    plt.rcParams['figure.dpi'] = config.get('dpi', 300)
    plt.rcParams['font.size'] = config.get('font_size', 12)
    
    logger.info("  Creating enhanced visualizations...")
    
    # 1. BUSCO Synteny Dotplot
    if config.get('generate_dotplots', True) and 'ortholog_df' in results_dict:
        create_busco_synteny_dotplot(results_dict['ortholog_df'], plots_dir, config)
    
    # 2. Ortholog Quality Distribution
    if 'ortholog_df' in results_dict and not results_dict['ortholog_df'].empty:
        create_ortholog_quality_plots(results_dict['ortholog_df'], plots_dir, config)
    
    # 3. Synteny Block Visualization
    if 'synteny_df' in results_dict and not results_dict['synteny_df'].empty:
        create_synteny_block_plots(results_dict['synteny_df'], plots_dir, config)
    
    # 4. Inversion Landscape Plot
    if 'inversion_df' in results_dict and not results_dict['inversion_df'].empty:
        create_inversion_landscape_plot(results_dict['inversion_df'], results_dict['ortholog_df'], plots_dir, config)
    
    # 5. Chromosome Rearrangement Network
    if 'rearrangement_df' in results_dict and not results_dict['rearrangement_df'].empty:
        create_rearrangement_network_plot(results_dict['rearrangement_df'], plots_dir, config)
    
    # 6. Summary Dashboard
    create_analysis_summary_dashboard(results_dict, plots_dir, config)
    
    logger.info(f"  ✓ Enhanced visualizations created in: {plots_dir}")


def add_synteny_block_lines(df, ax, config):
    """Add lines connecting synteny blocks on dotplot"""
    # Group by chromosome pairs and draw connecting lines for syntenic regions
    for (first_chr, second_chr), group in df.groupby(['first_chr', 'second_chr']):
        if len(group) >= 3:  # Only draw lines for substantial blocks
            syntenic_group = group[group['first_strand'] == group['second_strand']]
            if len(syntenic_group) >= 2:
                # Sort by first genome order
                sorted_group = syntenic_group.sort_values('first_order')
                ax.plot(sorted_group['first_order'], sorted_group['second_order'], 
                       'gray', alpha=0.5, linewidth=0.5, zorder=0)


def create_busco_synteny_dotplot(ortholog_df, plots_dir, config):
    """Create BUSCO-based synteny dotplot showing gene order relationships"""
    if len(ortholog_df) == 0:
        return
    
    plt.figure(figsize=config.get('dotplot_size', (12, 10)))
    
    # Prepare data for plotting
    df = ortholog_df.copy()
    
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
    df['first_order'] = df['busco_id'].map(first_order)
    df['second_order'] = df['busco_id'].map(second_order)
    
    # Remove genes without order (shouldn't happen but safety check)
    df = df.dropna(subset=['first_order', 'second_order'])
    
    if len(df) == 0:
        return
    
    # Determine point colors based on strand orientation
    syntenic_mask = df['first_strand'] == df['second_strand']
    
    # Plot syntenic points
    syntenic_data = df[syntenic_mask]
    if len(syntenic_data) > 0:
        alpha_vals = syntenic_data['confidence'].values if config.get('confidence_alpha', True) and 'confidence' in df.columns else 0.7
        plt.scatter(syntenic_data['first_order'], syntenic_data['second_order'], 
                   c=config.get('synteny_color', '#1f77b4'), alpha=alpha_vals, s=20, label='Syntenic')
    
    # Plot inverted points
    inverted_data = df[~syntenic_mask]
    if len(inverted_data) > 0:
        alpha_vals = inverted_data['confidence'].values if config.get('confidence_alpha', True) and 'confidence' in df.columns else 0.7
        plt.scatter(inverted_data['first_order'], inverted_data['second_order'], 
                   c=config.get('inversion_color', '#d62728'), alpha=alpha_vals, s=20, label='Inverted')
    
    # Add synteny block lines if enabled
    if config.get('show_synteny_blocks', True):
        add_synteny_block_lines(df, plt.gca(), config)
    
    # Formatting
    plt.xlabel('Gene Order in First Genome')
    plt.ylabel('Gene Order in Second Genome')
    plt.title('BUSCO Synteny Dotplot\n(Blue=Syntenic, Red=Inverted)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'busco_synteny_dotplot.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info("    ✓ BUSCO synteny dotplot created")


def create_ortholog_quality_plots(ortholog_df, plots_dir, config):
    """Create plots showing ortholog quality metrics"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Ortholog Quality Assessment', fontsize=14, fontweight='bold')
    
    # 1. Similarity distribution
    if 'similarity' in ortholog_df.columns:
        axes[0, 0].hist(ortholog_df['similarity'], bins=30, alpha=0.7, edgecolor='black', color='skyblue')
        axes[0, 0].axvline(ortholog_df['similarity'].mean(), color='red', linestyle='--', label=f'Mean: {ortholog_df["similarity"].mean():.3f}')
        axes[0, 0].set_xlabel('Similarity Score')
        axes[0, 0].set_ylabel('Number of Orthologs')
        axes[0, 0].set_title('Similarity Score Distribution')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Confidence distribution
    if 'confidence' in ortholog_df.columns:
        axes[0, 1].hist(ortholog_df['confidence'], bins=30, alpha=0.7, edgecolor='black', color='lightgreen')
        axes[0, 1].axvline(ortholog_df['confidence'].mean(), color='red', linestyle='--', label=f'Mean: {ortholog_df["confidence"].mean():.3f}')
        axes[0, 1].set_xlabel('Confidence Score')
        axes[0, 1].set_ylabel('Number of Orthologs')
        axes[0, 1].set_title('Confidence Score Distribution')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Alignment method usage
    if 'alignment_method' in ortholog_df.columns:
        method_counts = ortholog_df['alignment_method'].value_counts()
        colors = plt.cm.Set3(range(len(method_counts)))
        wedges, texts, autotexts = axes[1, 0].pie(method_counts.values, labels=method_counts.index, 
                                                  autopct='%1.1f%%', colors=colors)
        axes[1, 0].set_title('Alignment Methods Used')
    
    # 4. Gene length vs similarity scatter
    if 'similarity' in ortholog_df.columns and 'first_length' in ortholog_df.columns:
        scatter = axes[1, 1].scatter(ortholog_df['first_length'], ortholog_df['similarity'], 
                                   alpha=0.6, s=10, c=ortholog_df.get('confidence', 'blue'))
        axes[1, 1].set_xlabel('Gene Length (bp)')
        axes[1, 1].set_ylabel('Similarity Score')
        axes[1, 1].set_title('Gene Length vs Similarity')
        axes[1, 1].grid(True, alpha=0.3)
        
        if 'confidence' in ortholog_df.columns:
            cbar = plt.colorbar(scatter, ax=axes[1, 1])
            cbar.set_label('Confidence')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'ortholog_quality_assessment.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info("    ✓ Ortholog quality plots created")


def create_synteny_block_plots(synteny_df, plots_dir, config):
    """Create synteny block analysis plots"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Synteny Block Analysis', fontsize=14, fontweight='bold')
    
    # 1. Block size distribution
    axes[0, 0].hist(synteny_df['block_size'], bins=20, alpha=0.7, edgecolor='black', color='orange')
    axes[0, 0].set_xlabel('Block Size (genes)')
    axes[0, 0].set_ylabel('Number of Blocks')
    axes[0, 0].set_title('Synteny Block Size Distribution')
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. Synteny types
    if 'synteny_type' in synteny_df.columns:
        type_counts = synteny_df['synteny_type'].value_counts()
        colors = plt.cm.Pastel1(range(len(type_counts)))
        axes[0, 1].pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', colors=colors)
        axes[0, 1].set_title('Synteny Types')
    
    # 3. Position correlation distribution
    if 'position_correlation' in synteny_df.columns:
        axes[1, 0].hist(synteny_df['position_correlation'], bins=20, alpha=0.7, edgecolor='black', color='lightcoral')
        axes[1, 0].set_xlabel('Position Correlation')
        axes[1, 0].set_ylabel('Number of Blocks')
        axes[1, 0].set_title('Position Correlation Distribution')
        axes[1, 0].grid(True, alpha=0.3)
    
    # 4. Strand consistency distribution
    if 'strand_consistency' in synteny_df.columns:
        axes[1, 1].hist(synteny_df['strand_consistency'], bins=20, alpha=0.7, edgecolor='black', color='mediumpurple')
        axes[1, 1].set_xlabel('Strand Consistency')
        axes[1, 1].set_ylabel('Number of Blocks')
        axes[1, 1].set_title('Strand Consistency Distribution')
        axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'synteny_block_analysis.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info("    ✓ Synteny block plots created")


def create_inversion_landscape_plot(inversion_df, ortholog_df, plots_dir, config):
    """Create inversion landscape showing inversions across chromosomes"""
    if len(inversion_df) == 0:
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    fig.suptitle('Inversion Landscape Analysis', fontsize=14, fontweight='bold')
    
    # 1. Inversion size distribution
    if 'size_genes' in inversion_df.columns:
        ax1.hist(inversion_df['size_genes'], bins=20, alpha=0.7, edgecolor='black', color='crimson')
        ax1.set_xlabel('Inversion Size (genes)')
        ax1.set_ylabel('Number of Inversions')
        ax1.set_title('Inversion Size Distribution')
        ax1.grid(True, alpha=0.3)
    
    # 2. Inversion types
    if 'inversion_type' in inversion_df.columns:
        type_counts = inversion_df['inversion_type'].value_counts()
        bars = ax2.bar(range(len(type_counts)), type_counts.values, color=plt.cm.Reds(0.7))
        ax2.set_xticks(range(len(type_counts)))
        ax2.set_xticklabels(type_counts.index, rotation=45, ha='right')
        ax2.set_ylabel('Number of Inversions')
        ax2.set_title('Inversion Types')
        ax2.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, type_counts.values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{value}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'inversion_landscape.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info("    ✓ Inversion landscape plot created")


def create_rearrangement_network_plot(rearrangement_df, plots_dir, config):
    """Create network plot showing chromosome rearrangements"""
    if len(rearrangement_df) == 0:
        return
    
    try:
        import networkx as nx
        
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        fig.suptitle('Chromosome Rearrangement Analysis', fontsize=14, fontweight='bold')
        
        # 1. Rearrangement types pie chart
        if 'type' in rearrangement_df.columns:
            type_counts = rearrangement_df['type'].value_counts()
            colors = plt.cm.Set2(range(len(type_counts)))
            axes[0].pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', colors=colors)
            axes[0].set_title('Rearrangement Types')
        
        # 2. Simple network representation
        G = nx.Graph()
        
        for _, rearr in rearrangement_df.iterrows():
            if rearr['type'] == 'chromosome_split':
                # Add edges from source to all targets
                source = rearr['first_chr']
                targets = eval(rearr['second_chrs']) if isinstance(rearr['second_chrs'], str) else rearr['second_chrs']
                for target in targets:
                    G.add_edge(source, target, weight=1, type='split')
            
            elif rearr['type'] == 'chromosome_fusion':
                # Add edges from all sources to target
                sources = eval(rearr['first_chrs']) if isinstance(rearr['first_chrs'], str) else rearr['first_chrs']
                target = rearr['second_chr']
                for source in sources:
                    G.add_edge(source, target, weight=1, type='fusion')
        
        if len(G.nodes()) > 0:
            pos = nx.spring_layout(G, k=1, iterations=50)
            
            # Draw nodes
            nx.draw_networkx_nodes(G, pos, ax=axes[1], node_color='lightblue', 
                                 node_size=500, alpha=0.8)
            
            # Draw edges with different colors for different types
            split_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'split']
            fusion_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('type') == 'fusion']
            
            if split_edges:
                nx.draw_networkx_edges(G, pos, edgelist=split_edges, ax=axes[1], 
                                     edge_color='red', style='dashed', alpha=0.6)
            if fusion_edges:
                nx.draw_networkx_edges(G, pos, edgelist=fusion_edges, ax=axes[1], 
                                     edge_color='blue', style='solid', alpha=0.6)
            
            # Draw labels
            nx.draw_networkx_labels(G, pos, ax=axes[1], font_size=8)
            
            axes[1].set_title('Chromosome Rearrangement Network\n(Red=Split, Blue=Fusion)')
            axes[1].axis('off')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'rearrangement_network.png'), dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info("    ✓ Rearrangement network plot created")
        
    except ImportError:
        logger.warning("    NetworkX not available, skipping rearrangement network plot")


def create_analysis_summary_dashboard(results_dict, plots_dir, config):
    """Create summary dashboard with key metrics"""
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
    fig.suptitle('Genome Comparison Analysis Dashboard', fontsize=16, fontweight='bold')
    
    # 1. Overall statistics (top row, spans 2 columns)
    ax1 = fig.add_subplot(gs[0, :2])
    create_statistics_summary_plot(results_dict, ax1)
    
    # 2. Quality scores (top row, spans 2 columns)
    ax2 = fig.add_subplot(gs[0, 2:])
    create_quality_summary_plot(results_dict, ax2)
    
    # 3. Chromosome mapping overview (middle row)
    ax3 = fig.add_subplot(gs[1, :])
    create_chromosome_mapping_overview(results_dict, ax3)
    
    # 4. Bottom row: Individual component summaries
    ax4 = fig.add_subplot(gs[2, 0])
    create_synteny_summary_plot(results_dict, ax4)
    
    ax5 = fig.add_subplot(gs[2, 1])
    create_inversion_summary_plot(results_dict, ax5)
    
    ax6 = fig.add_subplot(gs[2, 2])
    create_rearrangement_summary_plot(results_dict, ax6)
    
    ax7 = fig.add_subplot(gs[2, 3])
    create_method_summary_plot(results_dict, ax7)
    
    plt.savefig(os.path.join(plots_dir, 'analysis_dashboard.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info("    ✓ Analysis summary dashboard created")


def create_statistics_summary_plot(results_dict, ax):
    """Create overall statistics summary"""
    stats = []
    labels = []
    
    if 'ortholog_df' in results_dict and not results_dict['ortholog_df'].empty:
        stats.extend([
            len(results_dict['ortholog_df']),
            results_dict['ortholog_df']['similarity'].mean() if 'similarity' in results_dict['ortholog_df'].columns else 0,
            results_dict['ortholog_df']['confidence'].mean() if 'confidence' in results_dict['ortholog_df'].columns else 0
        ])
        labels.extend(['Orthologs', 'Avg Similarity', 'Avg Confidence'])
    
    if 'synteny_df' in results_dict:
        stats.append(len(results_dict['synteny_df']))
        labels.append('Synteny Blocks')
    
    if 'inversion_df' in results_dict:
        stats.append(len(results_dict['inversion_df']))
        labels.append('Inversions')
    
    if 'rearrangement_df' in results_dict:
        stats.append(len(results_dict['rearrangement_df']))
        labels.append('Rearrangements')
    
    if stats:
        # Normalize stats for display (except counts)
        display_stats = []
        for i, (stat, label) in enumerate(zip(stats, labels)):
            if 'Avg' in label:
                display_stats.append(stat * 100)  # Convert to percentage
            else:
                display_stats.append(stat)
        
        colors = plt.cm.Set3(range(len(stats)))
        bars = ax.bar(labels, display_stats, color=colors, alpha=0.8)
        
        # Add value labels on bars
        for bar, stat, label in zip(bars, stats, labels):
            height = bar.get_height()
            if 'Avg' in label:
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                       f'{stat:.3f}', ha='center', va='bottom', fontsize=10)
            else:
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                       f'{int(stat)}', ha='center', va='bottom', fontsize=10)
        
        ax.set_title('Analysis Summary Statistics')
        ax.set_ylabel('Count / Score')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')


def create_quality_summary_plot(results_dict, ax):
    """Create quality summary plot"""
    if 'first_quality' in results_dict and 'second_quality' in results_dict:
        genomes = ['First Genome', 'Second Genome']
        scores = [results_dict['first_quality']['quality_score'], 
                 results_dict['second_quality']['quality_score']]
        classes = [results_dict['first_quality']['quality_class'], 
                  results_dict['second_quality']['quality_class']]
        
        # Color by quality class
        color_map = {'high': 'green', 'medium': 'orange', 'low': 'red', 'fragmented': 'darkred'}
        colors = [color_map.get(cls, 'gray') for cls in classes]
        
        bars = ax.bar(genomes, scores, color=colors, alpha=0.7)
        
        # Add score labels and quality class
        for bar, score, cls in zip(bars, scores, classes):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.02,
                   f'{score:.3f}\n({cls})', ha='center', va='bottom', fontsize=10)
        
        ax.set_title('Assembly Quality Scores')
        ax.set_ylabel('Quality Score')
        ax.set_ylim(0, 1.1)
        ax.grid(True, alpha=0.3)


def create_chromosome_mapping_overview(results_dict, ax):
    """Create chromosome mapping overview"""
    if 'ortholog_df' not in results_dict or results_dict['ortholog_df'].empty:
        ax.text(0.5, 0.5, 'No ortholog data available', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Chromosome Mapping Overview')
        return
    
    df = results_dict['ortholog_df']
    
    # Create chromosome mapping matrix
    first_chrs = sorted(df['first_chr'].unique())
    second_chrs = sorted(df['second_chr'].unique())
    
    # Create mapping matrix
    mapping_matrix = np.zeros((len(first_chrs), len(second_chrs)))
    
    for i, first_chr in enumerate(first_chrs):
        for j, second_chr in enumerate(second_chrs):
            count = len(df[(df['first_chr'] == first_chr) & (df['second_chr'] == second_chr)])
            mapping_matrix[i, j] = count
    
    # Create heatmap
    im = ax.imshow(mapping_matrix, cmap='YlOrRd', aspect='auto')
    
    # Set ticks and labels
    ax.set_xticks(range(len(second_chrs)))
    ax.set_yticks(range(len(first_chrs)))
    ax.set_xticklabels([chr[:10] for chr in second_chrs], rotation=45, ha='right')
    ax.set_yticklabels([chr[:10] for chr in first_chrs])
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.02)
    cbar.set_label('Number of Orthologs')
    
    ax.set_title('Chromosome Mapping Matrix')
    ax.set_xlabel('Second Genome Chromosomes')
    ax.set_ylabel('First Genome Chromosomes')


def create_synteny_summary_plot(results_dict, ax):
    """Create synteny summary plot"""
    if 'synteny_df' not in results_dict or results_dict['synteny_df'].empty:
        ax.text(0.5, 0.5, 'No synteny\ndata', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Synteny Summary')
        return
    
    synteny_df = results_dict['synteny_df']
    if 'synteny_type' in synteny_df.columns:
        type_counts = synteny_df['synteny_type'].value_counts()
        colors = plt.cm.Pastel1(range(len(type_counts)))
        ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.0f%%', 
               colors=colors, textprops={'fontsize': 8})
    ax.set_title('Synteny Types', fontsize=10)


def create_inversion_summary_plot(results_dict, ax):
    """Create inversion summary plot"""
    if 'inversion_df' not in results_dict or results_dict['inversion_df'].empty:
        ax.text(0.5, 0.5, 'No inversion\ndata', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Inversion Summary')
        return
    
    inversion_df = results_dict['inversion_df']
    if 'size_genes' in inversion_df.columns:
        ax.hist(inversion_df['size_genes'], bins=10, alpha=0.7, color='crimson', edgecolor='black')
        ax.set_xlabel('Size (genes)', fontsize=8)
        ax.set_ylabel('Count', fontsize=8)
    ax.set_title('Inversion Sizes', fontsize=10)


def create_rearrangement_summary_plot(results_dict, ax):
    """Create rearrangement summary plot"""
    if 'rearrangement_df' not in results_dict or results_dict['rearrangement_df'].empty:
        ax.text(0.5, 0.5, 'No rearrangement\ndata', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Rearrangement Summary')
        return
    
    rearr_df = results_dict['rearrangement_df']
    if 'type' in rearr_df.columns:
        type_counts = rearr_df['type'].value_counts()
        colors = plt.cm.Set2(range(len(type_counts)))
        ax.pie(type_counts.values, labels=type_counts.index, autopct='%1.0f%%', colors=colors, textprops={'fontsize': 8})
    ax.set_title('Rearrangement Types', fontsize=10)


def create_method_summary_plot(results_dict, ax):
    """Create alignment method summary plot"""
    if 'ortholog_df' not in results_dict or results_dict['ortholog_df'].empty:
        ax.text(0.5, 0.5, 'No method\ndata', ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Method Summary')
        return
    
    ortholog_df = results_dict['ortholog_df']
    if 'alignment_method' in ortholog_df.columns:
        method_counts = ortholog_df['alignment_method'].value_counts()
        colors = plt.cm.Set3(range(len(method_counts)))
        ax.pie(method_counts.values, labels=method_counts.index,
               autopct='%1.0f%%', colors=colors, textprops={'fontsize': 8})
    ax.set_title('Alignment Methods', fontsize=10)