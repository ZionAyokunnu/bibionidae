# =============================================================================
# Enhanced Visualization Integration (main_script_addition.py)
# =============================================================================

"""
Enhanced visualization functions to integrate with the main analysis script.
These functions replace the placeholder visualization code in the original script.
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

def create_enhanced_visualizations(output_dir, results_dict, config):
    """
    Create comprehensive visualizations using the visualization module.
    
    Args:
        output_dir: Output directory path
        results_dict: Dictionary containing analysis results
        config: Configuration dictionary
    """
    from visualization import AnalysisDashboard, SyntenyDotplot, SummaryPlotter, QualityPlotter, RearrangementNetworkPlotter
    from logger import get_logger
    
    logger = get_logger()
    
    plots_dir = Path(output_dir) / 'plots'
    plots_dir.mkdir(parents=True, exist_ok=True)
    
    # Set style and parameters
    plt.style.use('default')
    plt.rcParams['figure.dpi'] = config.get('dpi', 300)
    plt.rcParams['font.size'] = config.get('font_size', 12)
    
    logger.info("  Creating enhanced visualizations...")
    
    # Initialize visualization components
    dashboard = AnalysisDashboard(config)
    dotplot_maker = SyntenyDotplot(config)
    summary_plotter = SummaryPlotter(config)
    quality_plotter = QualityPlotter(config)
    network_plotter = RearrangementNetworkPlotter(config)
    
    plot_paths = {}
    
    try:
        # 1. Main Dashboard
        logger.info("    Creating main analysis dashboard...")
        dashboard_path = dashboard.create_comprehensive_dashboard(
            results_dict, plots_dir
        )
        if dashboard_path and 'main_dashboard' in dashboard_path:
            plot_paths.update(dashboard_path)
        
        # 2. BUSCO Synteny Dotplot
        if config.get('generate_dotplots', True) and 'ortholog_df' in results_dict:
            logger.info("    Creating BUSCO synteny dotplot...")
            dotplot_path = dotplot_maker.create_busco_synteny_dotplot(
                results_dict['ortholog_df'], 
                str(plots_dir / 'busco_synteny_dotplot.png')
            )
            if dotplot_path:
                plot_paths['synteny_dotplot'] = dotplot_path
        
        # 3. Ortholog Quality Distribution
        if 'ortholog_df' in results_dict and not results_dict['ortholog_df'].empty:
            logger.info("    Creating ortholog quality plots...")
            quality_path = summary_plotter.create_ortholog_quality_plots(
                results_dict['ortholog_df'], 
                str(plots_dir / 'ortholog_quality_assessment.png')
            )
            if quality_path:
                plot_paths['ortholog_quality'] = quality_path
        
        # 4. Synteny Block Visualization
        if 'synteny_df' in results_dict and not results_dict['synteny_df'].empty:
            logger.info("    Creating synteny block plots...")
            synteny_path = summary_plotter.create_synteny_analysis_plots(
                results_dict['synteny_df'], 
                str(plots_dir / 'synteny_block_analysis.png')
            )
            if synteny_path:
                plot_paths['synteny_analysis'] = synteny_path
        
        # 5. Assembly Quality Comparison
        if 'first_quality' in results_dict and 'second_quality' in results_dict:
            logger.info("    Creating assembly quality comparison...")
            quality_comp_path = quality_plotter.create_assembly_quality_comparison(
                results_dict['first_quality'], 
                results_dict['second_quality'],
                str(plots_dir / 'assembly_quality_comparison.png')
            )
            if quality_comp_path:
                plot_paths['quality_comparison'] = quality_comp_path
        
        # 6. Inversion Analysis
        if 'inversion_df' in results_dict and not results_dict['inversion_df'].empty:
            logger.info("    Creating inversion analysis plots...")
            inversion_path = create_inversion_landscape_plot(
                results_dict['inversion_df'], 
                results_dict.get('ortholog_df', pd.DataFrame()), 
                plots_dir, config
            )
            if inversion_path:
                plot_paths['inversion_landscape'] = inversion_path
        
        # 7. Chromosome Rearrangement Network
        if 'rearrangement_df' in results_dict and not results_dict['rearrangement_df'].empty:
            logger.info("    Creating rearrangement network plot...")
            network_path = network_plotter.create_chromosome_relationship_network(
                results_dict.get('ortholog_df', pd.DataFrame()),
                results_dict['rearrangement_df'],
                str(plots_dir / 'chromosome_network.png')
            )
            if network_path:
                plot_paths['chromosome_network'] = network_path
        
        # 8. Multi-chromosome dotplot
        if config.get('generate_multi_chromosome_plots', False) and 'ortholog_df' in results_dict:
            logger.info("    Creating multi-chromosome dotplot...")
            multi_plot_path = dotplot_maker.create_multi_chromosome_dotplot(
                results_dict['ortholog_df'],
                str(plots_dir / 'multi_chromosome_dotplot.png')
            )
            if multi_plot_path:
                plot_paths['multi_chromosome'] = multi_plot_path
        
        # 9. BUSCO Assessment (if BUSCO dataframes available)
        if 'first_busco_df' in results_dict and 'second_busco_df' in results_dict:
            logger.info("    Creating BUSCO assessment plots...")
            busco_path = quality_plotter.create_busco_assessment_plot(
                results_dict['first_busco_df'],
                results_dict['second_busco_df'],
                str(plots_dir / 'busco_assessment.png')
            )
            if busco_path:
                plot_paths['busco_assessment'] = busco_path
        
        logger.info(f"  ✓ Enhanced visualizations created: {len(plot_paths)} plots generated")
        logger.info(f"    Output directory: {plots_dir}")
        
        # Log individual plot paths
        for plot_type, path in plot_paths.items():
            logger.info(f"      {plot_type}: {Path(path).name}")
    
    except Exception as e:
        logger.error(f"Error creating visualizations: {e}")
        logger.error("Falling back to basic visualization functions...")
        
        # Fallback to basic visualizations
        create_basic_visualizations_fallback(plots_dir, results_dict, config)
    
    return plot_paths

def create_basic_visualizations_fallback(plots_dir, results_dict, config):
    """
    Fallback visualization functions using basic matplotlib.
    Used when the main visualization module fails.
    """
    from logger import get_logger
    logger = get_logger()
    
    logger.info("    Creating basic fallback visualizations...")
    
    # 1. Basic synteny dotplot
    if 'ortholog_df' in results_dict and not results_dict['ortholog_df'].empty:
        create_basic_synteny_dotplot(results_dict['ortholog_df'], plots_dir, config)
    
    # 2. Basic quality comparison
    if 'first_quality' in results_dict and 'second_quality' in results_dict:
        create_basic_quality_comparison(results_dict, plots_dir, config)
    
    # 3. Basic summary statistics
    create_basic_summary_plots(results_dict, plots_dir, config)
    
    logger.info("    ✓ Basic fallback visualizations completed")

def create_basic_synteny_dotplot(ortholog_df, plots_dir, config):
    """Create basic synteny dotplot using matplotlib."""
    try:
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
        
        # Remove genes without order
        df = df.dropna(subset=['first_order', 'second_order'])
        
        if len(df) == 0:
            return
        
        # Determine point colors based on strand orientation
        syntenic_mask = df['first_strand'] == df['second_strand']
        
        # Plot syntenic points
        syntenic_data = df[syntenic_mask]
        if len(syntenic_data) > 0:
            alpha_vals = syntenic_data.get('confidence', 0.7)
            if hasattr(alpha_vals, 'values'):
                alpha_vals = alpha_vals.values
            plt.scatter(syntenic_data['first_order'], syntenic_data['second_order'], 
                       c=config.get('synteny_color', '#1f77b4'), alpha=alpha_vals, s=20, label='Syntenic')
        
        # Plot inverted points
        inverted_data = df[~syntenic_mask]
        if len(inverted_data) > 0:
            alpha_vals = inverted_data.get('confidence', 0.7)
            if hasattr(alpha_vals, 'values'):
                alpha_vals = alpha_vals.values
            plt.scatter(inverted_data['first_order'], inverted_data['second_order'], 
                       c=config.get('inversion_color', '#d62728'), alpha=alpha_vals, s=20, label='Inverted')
        
        # Formatting
        plt.xlabel('Gene Order in First Genome')
        plt.ylabel('Gene Order in Second Genome')
        plt.title('BUSCO Synteny Dotplot\n(Blue=Syntenic, Red=Inverted)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(plots_dir / 'basic_synteny_dotplot.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error creating basic synteny dotplot: {e}")
        plt.close()

def create_basic_quality_comparison(results_dict, plots_dir, config):
    """Create basic quality comparison plot."""
    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Quality scores
        first_quality = results_dict['first_quality']
        second_quality = results_dict['second_quality']
        
        genomes = ['First Genome', 'Second Genome']
        scores = [first_quality['quality_score'], second_quality['quality_score']]
        classes = [first_quality['quality_class'], second_quality['quality_class']]
        
        # Color by quality class
        color_map = {'high': 'green', 'medium': 'orange', 'low': 'red', 'fragmented': 'darkred'}
        colors = [color_map.get(cls.lower(), 'gray') for cls in classes]
        
        bars = axes[0].bar(genomes, scores, color=colors, alpha=0.7)
        
        # Add score labels and quality class
        for bar, score, cls in zip(bars, scores, classes):
            height = bar.get_height()
            axes[0].text(bar.get_x() + bar.get_width()/2., height + 0.02,
                        f'{score:.3f}\n({cls})', ha='center', va='bottom', fontsize=10)
        
        axes[0].set_title('Assembly Quality Scores')
        axes[0].set_ylabel('Quality Score')
        axes[0].set_ylim(0, 1.1)
        axes[0].grid(True, alpha=0.3)
        
        # Basic metrics comparison
        metrics_to_compare = ['total_length', 'n_contigs', 'n50']
        metric_labels = ['Total Length (Mbp)', 'Contigs', 'N50 (Mbp)']
        
        first_metrics = first_quality.get('metrics', {})
        second_metrics = second_quality.get('metrics', {})
        
        x = np.arange(len(metrics_to_compare))
        width = 0.35
        
        first_values = []
        second_values = []
        
        for metric in metrics_to_compare:
            val1 = first_metrics.get(metric, 0)
            val2 = second_metrics.get(metric, 0)
            
            if metric in ['total_length', 'n50']:
                val1 = val1 / 1e6  # Convert to Mbp
                val2 = val2 / 1e6
            
            first_values.append(val1)
            second_values.append(val2)
        
        bars1 = axes[1].bar(x - width/2, first_values, width, label='First Genome', alpha=0.7)
        bars2 = axes[1].bar(x + width/2, second_values, width, label='Second Genome', alpha=0.7)
        
        axes[1].set_xlabel('Metrics')
        axes[1].set_ylabel('Values')
        axes[1].set_title('Assembly Statistics Comparison')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(metric_labels, rotation=45, ha='right')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'basic_quality_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error creating basic quality comparison: {e}")
        plt.close()

def create_basic_summary_plots(results_dict, plots_dir, config):
    """Create basic summary statistics plots."""
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Analysis Summary', fontsize=14, fontweight='bold')
        
        # 1. Overall statistics
        categories = []
        values = []
        colors = ['#1f77b4', '#ff7f0e', '#d62728', '#2ca02c']
        
        if 'ortholog_df' in results_dict:
            categories.append('Orthologs')
            values.append(len(results_dict['ortholog_df']))
        
        if 'synteny_df' in results_dict:
            categories.append('Synteny Blocks')
            values.append(len(results_dict['synteny_df']))
        
        if 'inversion_df' in results_dict:
            categories.append('Inversions')
            values.append(len(results_dict['inversion_df']))
        
        if 'rearrangement_df' in results_dict:
            categories.append('Rearrangements')
            values.append(len(results_dict['rearrangement_df']))
        
        if categories:
            bars = axes[0, 0].bar(categories, values, color=colors[:len(categories)], alpha=0.7)
            axes[0, 0].set_ylabel('Count')
            axes[0, 0].set_title('Overall Results')
            axes[0, 0].grid(True, alpha=0.3)
            
            # Add value labels
            for bar, value in zip(bars, values):
                height = bar.get_height()
                axes[0, 0].text(bar.get_x() + bar.get_width()/2., height + 0.01*max(values),
                               f'{value}', ha='center', va='bottom', fontweight='bold')
        else:
            axes[0, 0].text(0.5, 0.5, 'No data available', ha='center', va='center',
                           transform=axes[0, 0].transAxes, fontsize=12)
        
        # 2. Similarity distribution (if available)
        if 'ortholog_df' in results_dict and 'similarity' in results_dict['ortholog_df'].columns:
            similarity_data = results_dict['ortholog_df']['similarity']
            axes[0, 1].hist(similarity_data, bins=20, alpha=0.7, edgecolor='black', color='skyblue')
            axes[0, 1].axvline(similarity_data.mean(), color='red', linestyle='--', 
                              label=f'Mean: {similarity_data.mean():.3f}')
            axes[0, 1].set_xlabel('Similarity Score')
            axes[0, 1].set_ylabel('Count')
            axes[0, 1].set_title('Similarity Distribution')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        else:
            axes[0, 1].text(0.5, 0.5, 'No similarity data', ha='center', va='center',
                           transform=axes[0, 1].transAxes, fontsize=12)
            axes[0, 1].set_title('Similarity Distribution')
        
        # 3. Block size distribution (if available)
        if 'synteny_df' in results_dict and 'block_size' in results_dict['synteny_df'].columns:
            block_sizes = results_dict['synteny_df']['block_size']
            axes[1, 0].hist(block_sizes, bins=15, alpha=0.7, edgecolor='black', color='orange')
            axes[1, 0].axvline(block_sizes.mean(), color='red', linestyle='--',
                              label=f'Mean: {block_sizes.mean():.1f}')
            axes[1, 0].set_xlabel('Block Size (genes)')
            axes[1, 0].set_ylabel('Count')
            axes[1, 0].set_title('Synteny Block Sizes')
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        else:
            axes[1, 0].text(0.5, 0.5, 'No synteny data', ha='center', va='center',
                           transform=axes[1, 0].transAxes, fontsize=12)
            axes[1, 0].set_title('Synteny Block Sizes')
        
        # 4. Method usage (if available)
        if 'ortholog_df' in results_dict and 'alignment_method' in results_dict['ortholog_df'].columns:
            method_counts = results_dict['ortholog_df']['alignment_method'].value_counts()
            colors_pie = plt.cm.Set3(range(len(method_counts)))
            axes[1, 1].pie(method_counts.values, labels=method_counts.index, 
                          autopct='%1.1f%%', colors=colors_pie)
            axes[1, 1].set_title('Alignment Methods')
        else:
            axes[1, 1].text(0.5, 0.5, 'No method data', ha='center', va='center',
                           transform=axes[1, 1].transAxes, fontsize=12)
            axes[1, 1].set_title('Alignment Methods')
        
        plt.tight_layout()
        plt.savefig(plots_dir / 'basic_summary_plots.png', dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error creating basic summary plots: {e}")
        plt.close()

def create_inversion_landscape_plot(inversion_df, ortholog_df, plots_dir, config):
    """Create inversion landscape showing inversions across chromosomes."""
    if len(inversion_df) == 0:
        return None
    
    try:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
        fig.suptitle('Inversion Landscape Analysis', fontsize=14, fontweight='bold')
        
        # 1. Inversion size distribution
        if 'size_genes' in inversion_df.columns:
            ax1.hist(inversion_df['size_genes'], bins=20, alpha=0.7, 
                    edgecolor='black', color='crimson')
            ax1.set_xlabel('Inversion Size (genes)')
            ax1.set_ylabel('Number of Inversions')
            ax1.set_title('Inversion Size Distribution')
            ax1.grid(True, alpha=0.3)
            
            # Add mean line
            mean_size = inversion_df['size_genes'].mean()
            ax1.axvline(mean_size, color='blue', linestyle='--',
                       label=f'Mean: {mean_size:.1f}')
            ax1.legend()
        
        # 2. Inversion types or chromosome distribution
        if 'inversion_type' in inversion_df.columns:
            type_counts = inversion_df['inversion_type'].value_counts()
            bars = ax2.bar(range(len(type_counts)), type_counts.values, 
                          color=plt.cm.Reds(np.linspace(0.4, 0.8, len(type_counts))))
            ax2.set_xticks(range(len(type_counts)))
            ax2.set_xticklabels([label.replace('_', '\n') for label in type_counts.index], 
                              rotation=45, ha='right')
            ax2.set_ylabel('Number of Inversions')
            ax2.set_title('Inversion Types')
            ax2.grid(True, alpha=0.3)
            
            # Add value labels on bars
            for bar, value in zip(bars, type_counts.values):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{value}', ha='center', va='bottom')
        elif 'first_chr' in inversion_df.columns:
            # Show chromosome distribution instead
            chr_counts = inversion_df['first_chr'].value_counts()
            colors = plt.cm.Set3(range(len(chr_counts)))
            ax2.pie(chr_counts.values, labels=chr_counts.index,
                   autopct='%1.1f%%', colors=colors)
            ax2.set_title('Inversions by Chromosome')
        
        plt.tight_layout()
        output_path = plots_dir / 'inversion_landscape.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return str(output_path)
        
    except Exception as e:
        print(f"Error creating inversion landscape plot: {e}")
        plt.close()
        return None