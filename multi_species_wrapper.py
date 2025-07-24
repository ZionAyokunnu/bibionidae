#!/usr/bin/env python3
"""
Enhanced Multi-Species Wrapper with Complete Visualization Suite
Boolean switches for all plot types
"""

import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns

# Import your working modules
sys.path.insert(0, str(Path(__file__).parent))

from genome_inversion_analyser.main import run_complete_enhanced_analysis_with_registry
from genome_inversion_analyser.config import ENHANCED_HYBRID_CONFIG
from genome_inversion_analyser.registry import FileRegistry
from genome_inversion_analyser.phylogenetic import PhylogeneticIntegrator

# VISUALIZATION CONTROL SETTINGS
VISUALIZATION_CONFIG = {
    'create_traditional_dendrograms': True,      # Classic phylogenetic trees
    'create_syri_style_plots': True,            # SyRI-style synteny plots
    'create_multi_species_matrix': True,        # Comparison matrix heatmaps
    'create_circular_phylogeny': True,          # Circular phylogenetic plots
    'create_inversion_landscape': True,         # Multi-species inversion overview
    'create_summary_dashboard': True,           # Complete analysis dashboard
    'plot_format': 'png',                       # 'png', 'pdf', 'svg'
    'plot_dpi': 300,                           # Resolution
    'figure_size_large': (16, 12),            # Large plots
    'figure_size_medium': (12, 8),            # Medium plots
    'figure_size_small': (8, 6),              # Small plots
}

def load_species_config():
    """Load our species configuration"""
    with open('species_config.json', 'r') as f:
        config = json.load(f)
    return config

def create_pairwise_config(species1, species2):
    """Create config for one pairwise comparison"""
    pairwise_config = ENHANCED_HYBRID_CONFIG.copy()
    pair_name = f"{species1['name']}_vs_{species2['name']}"
    
    pairwise_config.update({
        'first_fasta_path': species1['fasta'],
        'second_fasta_path': species2['fasta'], 
        'first_busco_path': species1['busco'],
        'second_busco_path': species2['busco'],
        'base_output_dir': f'pairwise_results/{pair_name}',
        'synteny_analysis_csv': f'pairwise_results/{pair_name}/synteny_analysis.csv',
        'inversion_summary_csv': f'pairwise_results/{pair_name}/inversion_summary.csv',
        'chromosome_rearrangements_csv': f'pairwise_results/{pair_name}/chromosome_rearrangements.csv',
        'paralog_analysis_csv': f'pairwise_results/{pair_name}/paralog_analysis.csv',
        'quality_report_csv': f'pairwise_results/{pair_name}/quality_report.csv',
    })
    
    return pairwise_config

def create_traditional_dendrograms(distance_matrices, trees, species_names, output_dir):
    """Create traditional phylogenetic tree dendrograms"""
    if not VISUALIZATION_CONFIG['create_traditional_dendrograms']:
        return
    
    print("  🌳 Creating traditional dendrograms...")
    
    try:
        from scipy.cluster.hierarchy import dendrogram
        
        dendro_dir = output_dir / 'phylogenetic_trees'
        dendro_dir.mkdir(exist_ok=True)
        
        for method_name, distance_matrix in distance_matrices.items():
            for linkage_method in ['ward', 'average', 'complete']:
                tree_key = f"{method_name}_{linkage_method}"
                
                if tree_key in trees:
                    fig, ax = plt.subplots(figsize=VISUALIZATION_CONFIG['figure_size_medium'])
                    
                    # Create dendrogram
                    tree_data = trees[tree_key]
                    linkage_matrix = np.array(tree_data['linkage_matrix'])
                    
                    dendrogram(linkage_matrix, 
                             labels=species_names,
                             orientation='top',
                             distance_sort='descending',
                             ax=ax)
                    
                    ax.set_title(f'Bibionidae Phylogeny\n({method_name} distance, {linkage_method} linkage)', 
                               fontsize=14, fontweight='bold')
                    ax.set_xlabel('Species', fontsize=12)
                    ax.set_ylabel('Distance', fontsize=12)
                    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
                    
                    # Add grid
                    ax.grid(True, alpha=0.3)
                    
                    plt.tight_layout()
                    
                    tree_file = dendro_dir / f'dendrogram_{tree_key}.{VISUALIZATION_CONFIG["plot_format"]}'
                    plt.savefig(tree_file, dpi=VISUALIZATION_CONFIG['plot_dpi'], bbox_inches='tight')
                    plt.close()
                    
                    print(f"    • {tree_file.name}")
        
        print(f"  ✅ Dendrograms saved to {dendro_dir}")
        
    except Exception as e:
        print(f"  ❌ Dendrogram creation failed: {e}")

def create_syri_style_plots(all_results, output_dir):
    """Create SyRI-style synteny plots for each species pair"""
    if not VISUALIZATION_CONFIG['create_syri_style_plots']:
        return
    
    print("  📊 Creating SyRI-style synteny plots...")
    
    try:
        from genome_inversion_analyser.visualization.syri_integration import SyRIIntegrator
        
        syri_dir = output_dir / 'syri_plots'
        syri_integrator = SyRIIntegrator(syri_dir)
        
        plots_created = 0
        for pair_name, pair_results in all_results.items():
            if 'error' not in pair_results:
                try:
                    species1, species2 = pair_results['species_pair']
                    
                    # Create SyRI-compatible output
                    syri_output = syri_integrator.create_syri_compatible_output(
                        pair_results['ortholog_df'],
                        pair_results['synteny_df'],
                        pair_results['inversion_df'],
                        species1, species2
                    )
                    
                    # Create SyRI-style plot
                    syri_plot = syri_integrator.create_syri_style_plot(
                        pair_results['ortholog_df'],
                        pair_results['synteny_df'], 
                        pair_results['inversion_df'],
                        species1, species2
                    )
                    
                    # Create circular plot
                    circular_plot = syri_integrator.create_circular_synteny_plot(
                        pair_results['ortholog_df'],
                        pair_results['inversion_df'],
                        species1, species2
                    )
                    
                    plots_created += 1
                    print(f"    • {pair_name}: synteny + circular plots")
                    
                except Exception as e:
                    print(f"    ❌ {pair_name}: {e}")
        
        print(f"  ✅ SyRI plots created for {plots_created} species pairs")
        
    except Exception as e:
        print(f"  ❌ SyRI plot creation failed: {e}")

def create_multi_species_matrix(distance_matrices, species_names, all_results, output_dir):
    """Create multi-species comparison matrices"""
    if not VISUALIZATION_CONFIG['create_multi_species_matrix']:
        return
    
    print("  🔢 Creating multi-species comparison matrices...")
    
    try:
        matrix_dir = output_dir / 'comparison_matrices'
        matrix_dir.mkdir(exist_ok=True)
        
        # 1. Distance matrices heatmaps
        for method_name, distance_matrix in distance_matrices.items():
            fig, ax = plt.subplots(figsize=VISUALIZATION_CONFIG['figure_size_medium'])
            
            # Create distance matrix heatmap
            sns.heatmap(distance_matrix, 
                       xticklabels=species_names,
                       yticklabels=species_names,
                       annot=True,
                       fmt='.3f',
                       cmap='viridis',
                       ax=ax)
            
            ax.set_title(f'Species Distance Matrix ({method_name})', fontsize=14, fontweight='bold')
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
            plt.setp(ax.get_yticklabels(), rotation=0)
            
            plt.tight_layout()
            
            matrix_file = matrix_dir / f'distance_matrix_{method_name}.{VISUALIZATION_CONFIG["plot_format"]}'
            plt.savefig(matrix_file, dpi=VISUALIZATION_CONFIG['plot_dpi'], bbox_inches='tight')
            plt.close()
            
            print(f"    • {matrix_file.name}")
        
        # 2. Inversion count matrix
        n_species = len(species_names)
        inversion_matrix = np.zeros((n_species, n_species))
        
        for pair_name, pair_results in all_results.items():
            if 'error' not in pair_results:
                species1, species2 = pair_results['species_pair']
                i = species_names.index(species1)
                j = species_names.index(species2)
                inv_count = pair_results['inversion_count']
                inversion_matrix[i, j] = inv_count
                inversion_matrix[j, i] = inv_count
        
        fig, ax = plt.subplots(figsize=VISUALIZATION_CONFIG['figure_size_medium'])
        sns.heatmap(inversion_matrix,
                   xticklabels=species_names,
                   yticklabels=species_names,
                   annot=True,
                   fmt='d',
                   cmap='Reds',
                   ax=ax)
        
        ax.set_title('Inversion Count Matrix', fontsize=14, fontweight='bold')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        plt.setp(ax.get_yticklabels(), rotation=0)
        
        plt.tight_layout()
        
        inv_matrix_file = matrix_dir / f'inversion_matrix.{VISUALIZATION_CONFIG["plot_format"]}'
        plt.savefig(inv_matrix_file, dpi=VISUALIZATION_CONFIG['plot_dpi'], bbox_inches='tight')
        plt.close()
        
        print(f"    • {inv_matrix_file.name}")
        print(f"  ✅ Comparison matrices saved to {matrix_dir}")
        
    except Exception as e:
        print(f"  ❌ Matrix creation failed: {e}")

def create_circular_phylogeny(distance_matrices, species_names, output_dir):
    """Create circular phylogenetic plots"""
    if not VISUALIZATION_CONFIG['create_circular_phylogeny']:
        return
    
    print("  ⭕ Creating circular phylogenetic plots...")
    
    try:
        circular_dir = output_dir / 'circular_plots'
        circular_dir.mkdir(exist_ok=True)
        
        for method_name, distance_matrix in distance_matrices.items():
            fig, ax = plt.subplots(figsize=VISUALIZATION_CONFIG['figure_size_medium'], 
                                 subplot_kw=dict(projection='polar'))
            
            n_species = len(species_names)
            angles = np.linspace(0, 2*np.pi, n_species, endpoint=False)
            
            # Plot species as points on circle
            radius = 1.0
            for i, (angle, species) in enumerate(zip(angles, species_names)):
                ax.plot([angle], [radius], 'o', markersize=15, color=f'C{i}')
                ax.text(angle, radius + 0.1, species, ha='center', va='center',
                       fontsize=10, fontweight='bold')
            
            # Draw connections based on distances
            for i in range(n_species):
                for j in range(i+1, n_species):
                    distance = distance_matrix[i, j]
                    # Normalize distance to line width (invert so closer = thicker)
                    max_dist = np.max(distance_matrix)
                    line_width = 5 * (1 - distance / max_dist) if max_dist > 0 else 1
                    
                    # Draw connection
                    ax.plot([angles[i], angles[j]], [radius, radius], 
                           color='gray', alpha=0.6, linewidth=line_width)
            
            ax.set_ylim(0, 1.3)
            ax.set_title(f'Circular Phylogeny ({method_name})', fontsize=14, fontweight='bold', pad=20)
            ax.grid(False)
            ax.set_rticks([])
            ax.set_thetagrids([])
            
            circular_file = circular_dir / f'circular_phylogeny_{method_name}.{VISUALIZATION_CONFIG["plot_format"]}'
            plt.savefig(circular_file, dpi=VISUALIZATION_CONFIG['plot_dpi'], bbox_inches='tight')
            plt.close()
            
            print(f"    • {circular_file.name}")
        
        print(f"  ✅ Circular plots saved to {circular_dir}")
        
    except Exception as e:
        print(f"  ❌ Circular plot creation failed: {e}")

def create_inversion_landscape(all_results, species_stats, output_dir):
    """Create multi-species inversion landscape overview"""
    if not VISUALIZATION_CONFIG['create_inversion_landscape']:
        return
    
    print("  🏔️ Creating inversion landscape overview...")
    
    try:
        landscape_dir = output_dir / 'inversion_landscape'
        landscape_dir.mkdir(exist_ok=True)
        
        # Collect inversion data across all species
        species_inv_data = {}
        for species_name in species_stats.keys():
            species_inv_data[species_name] = {
                'total_inversions': 0,
                'avg_inversion_size': 0,
                'genome_size': species_stats[species_name]['quality']['metrics']['total_length'],
                'quality_class': species_stats[species_name]['quality']['quality_class']
            }
        
        # Count inversions per species
        for pair_results in all_results.values():
            if 'error' not in pair_results:
                species1, species2 = pair_results['species_pair']
                inv_count = pair_results['inversion_count']
                
                species_inv_data[species1]['total_inversions'] += inv_count
                species_inv_data[species2]['total_inversions'] += inv_count
                
                if not pair_results['inversion_df'].empty:
                    avg_size = pair_results['inversion_df']['size_genes'].mean()
                    species_inv_data[species1]['avg_inversion_size'] = avg_size
                    species_inv_data[species2]['avg_inversion_size'] = avg_size
        
        # Create landscape plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=VISUALIZATION_CONFIG['figure_size_large'])
        
        species_names = list(species_inv_data.keys())
        
        # Plot 1: Total inversions
        inv_counts = [species_inv_data[sp]['total_inversions'] for sp in species_names]
        ax1.bar(species_names, inv_counts, color=['C0', 'C1', 'C2', 'C3'][:len(species_names)])
        ax1.set_title('Total Inversions per Species', fontweight='bold')
        ax1.set_ylabel('Inversion Count')
        plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')
        
        # Plot 2: Inversion rate per Mb
        rates = [inv_counts[i] / (species_inv_data[sp]['genome_size'] / 1_000_000) 
                for i, sp in enumerate(species_names)]
        ax2.bar(species_names, rates, color=['C0', 'C1', 'C2', 'C3'][:len(species_names)])
        ax2.set_title('Inversion Rate per Mb', fontweight='bold')
        ax2.set_ylabel('Inversions per Mb')
        plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')
        
        # Plot 3: Genome size
        genome_sizes = [species_inv_data[sp]['genome_size'] / 1_000_000 for sp in species_names]
        ax3.bar(species_names, genome_sizes, color=['C0', 'C1', 'C2', 'C3'][:len(species_names)])
        ax3.set_title('Genome Size (Mb)', fontweight='bold')
        ax3.set_ylabel('Size (Mb)')
        plt.setp(ax3.get_xticklabels(), rotation=45, ha='right')
        
        # Plot 4: Quality comparison
        quality_map = {'high': 3, 'medium': 2, 'low': 1, 'fragmented': 0}
        quality_scores = [quality_map.get(species_inv_data[sp]['quality_class'], 1) for sp in species_names]
        ax4.bar(species_names, quality_scores, color=['C0', 'C1', 'C2', 'C3'][:len(species_names)])
        ax4.set_title('Assembly Quality', fontweight='bold')
        ax4.set_ylabel('Quality Score')
        ax4.set_yticks([0, 1, 2, 3])
        ax4.set_yticklabels(['Fragmented', 'Low', 'Medium', 'High'])
        plt.setp(ax4.get_xticklabels(), rotation=45, ha='right')
        
        plt.suptitle('Bibionidae Inversion Landscape Overview', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        landscape_file = landscape_dir / f'inversion_landscape.{VISUALIZATION_CONFIG["plot_format"]}'
        plt.savefig(landscape_file, dpi=VISUALIZATION_CONFIG['plot_dpi'], bbox_inches='tight')
        plt.close()
        
        print(f"    • {landscape_file.name}")
        print(f"  ✅ Inversion landscape saved to {landscape_dir}")
        
    except Exception as e:
        print(f"  ❌ Inversion landscape creation failed: {e}")

def create_summary_dashboard(all_results, species_stats, distance_matrices, trees, output_dir):
    """Create comprehensive analysis dashboard"""
    if not VISUALIZATION_CONFIG['create_summary_dashboard']:
        return
    
    print("  📊 Creating summary dashboard...")
    
    try:
        dashboard_file = output_dir / f'analysis_dashboard.{VISUALIZATION_CONFIG["plot_format"]}'
        
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
        
        # Summary statistics
        ax1 = fig.add_subplot(gs[0, :2])
        successful_pairs = len([r for r in all_results.values() if 'error' not in r])
        total_pairs = len(all_results)
        total_inversions = sum(r.get('inversion_count', 0) for r in all_results.values() if 'error' not in r)
        
        summary_text = f"""
        BIBIONIDAE MULTI-SPECIES ANALYSIS SUMMARY
        
        Species Analyzed: {len(species_stats)}
        Successful Pairwise Comparisons: {successful_pairs}/{total_pairs}
        Total Inversions Detected: {total_inversions:,}
        Distance Methods: {len(distance_matrices)}
        Phylogenetic Trees: {len(trees)}
        """
        
        ax1.text(0.1, 0.5, summary_text, transform=ax1.transAxes, fontsize=12,
                verticalalignment='center', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1)
        ax1.axis('off')
        ax1.set_title('Analysis Summary', fontsize=14, fontweight='bold')
        
        # Species quality overview
        ax2 = fig.add_subplot(gs[0, 2:])
        species_names = list(species_stats.keys())
        quality_classes = [species_stats[sp]['quality']['quality_class'] for sp in species_names]
        quality_counts = pd.Series(quality_classes).value_counts()
        
        ax2.pie(quality_counts.values, labels=quality_counts.index, autopct='%1.0f%%')
        ax2.set_title('Assembly Quality Distribution', fontsize=14, fontweight='bold')
        
        # Add more dashboard elements...
        # (Distance matrix, tree preview, etc.)
        
        plt.suptitle('Bibionidae Multi-Species Analysis Dashboard', fontsize=18, fontweight='bold')
        plt.savefig(dashboard_file, dpi=VISUALIZATION_CONFIG['plot_dpi'], bbox_inches='tight')
        plt.close()
        
        print(f"    • {dashboard_file.name}")
        
    except Exception as e:
        print(f"  ❌ Dashboard creation failed: {e}")

def run_multi_species_analysis():
    """Run multi-species analysis using your existing pairwise code"""
    
    print("🧬 MULTI-SPECIES PHYLOGENETIC ANALYSIS")
    print("Using your proven pairwise analysis code")
    print("=" * 60)
    
    # Load configuration
    config = load_species_config()
    species_data = config['species_data']
    
    print(f"📊 Species to analyze: {len(species_data)}")
    for i, species in enumerate(species_data, 1):
        print(f"  {i}. {species['name']}")
    
    # [File verification code - same as before]
    print(f"\n🔍 Verifying files exist...")
    all_files_exist = True
    for species in species_data:
        fasta_path = Path(species['fasta'])
        busco_path = Path(species['busco'])
        
        if not fasta_path.exists():
            print(f"❌ FASTA not found: {fasta_path}")
            all_files_exist = False
        if not busco_path.exists():
            print(f"❌ BUSCO not found: {busco_path}")
            all_files_exist = False
    
    if not all_files_exist:
        print("❌ Some files missing. Please check paths in species_config.json")
        return False
    
    print("✅ All files found!")
    
    # [Pairwise analysis code - same as before]
    results_dir = Path("pairwise_results")
    results_dir.mkdir(exist_ok=True)
    
    species_pairs = list(combinations(species_data, 2))
    total_comparisons = len(species_pairs)
    
    print(f"\n🔄 Running {total_comparisons} pairwise comparisons:")
    
    all_results = {}
    species_stats = {}
    
    for i, (species1, species2) in enumerate(species_pairs, 1):
        pair_name = f"{species1['name']}_vs_{species2['name']}"
        print(f"\n📊 [{i}/{total_comparisons}] {pair_name}")
        
        try:
            pairwise_config = create_pairwise_config(species1, species2)
            print(f"  • Running your pairwise analysis...")
            
            # DIRECT call to your main function
            results = run_complete_enhanced_analysis_with_registry(pairwise_config)
            
            if results is None:
                print(f"  ❌ Analysis returned None")
                continue
            
            # [Results processing - same as before]
            inversion_count = len(results['inversion_df'])
            synteny_blocks = len(results['synteny_df'])
            ortholog_count = len(results['ortholog_df'])
            avg_similarity = results['ortholog_df']['similarity'].mean() if len(results['ortholog_df']) > 0 else 0
            
            print(f"  • Results: {ortholog_count} orthologs, {inversion_count} inversions, {synteny_blocks} synteny blocks")
            
            all_results[pair_name] = {
                'species_pair': (species1['name'], species2['name']),
                'inversion_df': results['inversion_df'],
                'ortholog_df': results['ortholog_df'],
                'synteny_df': results['synteny_df'],
                'inversion_count': inversion_count,
                'avg_similarity': avg_similarity,
                'full_results': results
            }
            
            # Store species stats
            for species_info, quality_key in [(species1, 'first_quality'), (species2, 'second_quality')]:
                if species_info['name'] not in species_stats:
                    species_stats[species_info['name']] = {
                        'quality': results[quality_key],
                        'metadata': species_info.get('metadata', {}),
                        'fasta': species_info['fasta'],
                        'busco': species_info['busco']
                    }
            
            print(f"  ✅ {pair_name} completed successfully")
            
        except Exception as e:
            print(f"  ❌ Failed: {e}")
            all_results[pair_name] = {'error': str(e)}
            continue
    
    successful_comparisons = len([r for r in all_results.values() if 'error' not in r])
    print(f"\n🎉 Completed {successful_comparisons}/{total_comparisons} pairwise comparisons!")
    
    if successful_comparisons < 1:
        print("❌ No successful comparisons. Cannot build phylogenetic tree.")
        return False
    
    # [Phylogenetic analysis - same as before]
    print(f"\n🌳 Building phylogenetic trees...")
    
    try:
        phylo_output_dir = Path("bibionidae_phylogenetic_analysis")
        registry = FileRegistry(phylo_output_dir, project_name="bibionidae_phylogeny")
        phylo_analyzer = PhylogeneticIntegrator(registry)
        
        # [Add species data - same as before]
        for species_name, stats in species_stats.items():
            species_inversions = []
            species_orthologs = []
            
            for pair_name, pair_results in all_results.items():
                if 'error' not in pair_results and species_name in pair_results['species_pair']:
                    inv_df = pair_results['inversion_df'].copy()
                    if not inv_df.empty:
                        inv_df['comparison'] = pair_name
                        species_inversions.extend(inv_df.to_dict('records'))
            
            species_inversion_df = pd.DataFrame(species_inversions) if species_inversions else pd.DataFrame()
            
            genome_size = stats['quality']['metrics'].get('total_length', 1000000)
            contextual_metrics = {
                'rate_metrics': {
                    'inversions_per_mb': len(species_inversions) / (genome_size / 1_000_000)
                }
            }
            
            phylo_analyzer.add_species_data(
                species_name,
                species_inversion_df,
                stats['quality']['metrics'],
                contextual_metrics
            )
            
            print(f"  • Added {species_name}: {len(species_inversions)} inversions")
        
        # Generate phylogenetic analysis
        print(f"  • Computing distance matrices...")
        distance_matrices = phylo_analyzer.compute_distance_matrices()
        
        print(f"  • Inferring phylogenetic trees...")
        trees = phylo_analyzer.infer_phylogenetic_trees()
        
        print(f"  • Analyzing phylogenetic signal...")
        signal_results = phylo_analyzer.analyze_phylogenetic_signal()
        
        print(f"  • Generating comprehensive summary...")
        phylo_summary = phylo_analyzer.generate_phylogenetic_summary()
        
        # NEW: Create all visualizations with boolean controls
        print(f"\n🎨 Creating comprehensive visualizations...")
        species_names = list(species_stats.keys())
        
        create_traditional_dendrograms(distance_matrices, trees, species_names, phylo_output_dir)
        create_syri_style_plots(all_results, phylo_output_dir)
        create_multi_species_matrix(distance_matrices, species_names, all_results, phylo_output_dir)
        create_circular_phylogeny(distance_matrices, species_names, phylo_output_dir)
        create_inversion_landscape(all_results, species_stats, phylo_output_dir)
        create_summary_dashboard(all_results, species_stats, distance_matrices, trees, phylo_output_dir)
        
        print(f"\n🎉 COMPLETE ANALYSIS FINISHED!")
        print("=" * 60)
        
        # [Results summary - same as before]
        print(f"\n📊 COMPREHENSIVE RESULTS:")
        print(f"  • Species analyzed: {len(species_stats)}")
        print(f"  • Successful pairwise comparisons: {successful_comparisons}")
        print(f"  • Distance matrices computed: {len(distance_matrices)}")
        print(f"  • Phylogenetic trees inferred: {len(trees)}")
        
        print(f"\n🎨 VISUALIZATIONS CREATED:")
        print(f"  • Traditional dendrograms: {VISUALIZATION_CONFIG['create_traditional_dendrograms']}")
        print(f"  • SyRI-style plots: {VISUALIZATION_CONFIG['create_syri_style_plots']}")
        print(f"  • Multi-species matrices: {VISUALIZATION_CONFIG['create_multi_species_matrix']}")
        print(f"  • Circular phylogeny: {VISUALIZATION_CONFIG['create_circular_phylogeny']}")
        print(f"  • Inversion landscape: {VISUALIZATION_CONFIG['create_inversion_landscape']}")
        print(f"  • Summary dashboard: {VISUALIZATION_CONFIG['create_summary_dashboard']}")
        
        print(f"\n📁 OUTPUT LOCATIONS:")
        print(f"  • Pairwise results: pairwise_results/")
        print(f"  • Phylogenetic analysis: {phylo_output_dir}/")
        print(f"  • All visualizations: {phylo_output_dir}/*/")
        
        print(f"\n🎯 SUCCESS! Complete Bibionidae analysis with all visualizations!")
        
        return True
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False