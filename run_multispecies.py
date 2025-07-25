#!/usr/bin/env python3
"""
Fixed Multi-Species Analysis - Based on Your Working Version
Just adds the missing SyRI integration calls
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
        'base_output_dir': f'pairwise_results/{pair_name}',  # Keep simple path that works
        'synteny_analysis_csv': f'pairwise_results/{pair_name}/synteny_analysis.csv',
        'inversion_summary_csv': f'pairwise_results/{pair_name}/inversion_summary.csv',
        'chromosome_rearrangements_csv': f'pairwise_results/{pair_name}/chromosome_rearrangements.csv',
        'paralog_analysis_csv': f'pairwise_results/{pair_name}/paralog_analysis.csv',
        'quality_report_csv': f'pairwise_results/{pair_name}/quality_report.csv',
    })
    
    return pairwise_config

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
    
    # File verification
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
    
    # Clean old results
    results_dir = Path("pairwise_results")
    if results_dir.exists():
        import shutil
        shutil.rmtree(results_dir)
        print("🗑️ Cleaned old pairwise results")
    
    phylo_dir = Path("bibionidae_phylogenetic_analysis")
    if phylo_dir.exists():
        import shutil
        shutil.rmtree(phylo_dir)
        print("🗑️ Cleaned old phylogenetic results")
    
    # Pairwise analysis - YOUR WORKING CODE
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
            
            # Call your main function - THIS WORKS
            results = run_complete_enhanced_analysis_with_registry(pairwise_config)
            
            if results is None:
                print(f"  ❌ Analysis returned None")
                continue
            
            # Process results
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
            import traceback
            traceback.print_exc()
            continue
    
    successful_comparisons = len([r for r in all_results.values() if 'full_results' in r])
    print(f"\n🎉 Completed {successful_comparisons}/{total_comparisons} pairwise comparisons!")
    
    if successful_comparisons < 1:
        print("❌ No successful comparisons. Cannot build phylogenetic tree.")
        return False
    
    # Phylogenetic analysis - WORKING CODE
    print(f"\n🌳 Building phylogenetic trees...")
    
    try:
        phylo_output_dir = Path("bibionidae_phylogenetic_analysis")
        registry = FileRegistry(phylo_output_dir, project_name="bibionidae_phylogeny")
        phylo_analyzer = PhylogeneticIntegrator(registry, {'min_species': 2})  # Lower threshold
        
        # Add species data
        for species_name, stats in species_stats.items():
            species_inversions = []
            
            for pair_name, pair_results in all_results.items():
                if 'full_results' in pair_results and species_name in pair_results['species_pair']:
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
        distance_matrices = phylo_analyzer.compute_distance_matrices(['jaccard', 'euclidean'])
        
        print(f"  • Inferring phylogenetic trees...")
        trees = phylo_analyzer.infer_phylogenetic_trees(['ward'])
        
        print(f"  ✅ Created {len(distance_matrices)} distance matrices, {len(trees)} trees")
        
        # CREATE SYRI PLOTS - THE MISSING PIECE
        print(f"\n📊 Creating SyRI visualizations...")
        
        try:
            from genome_inversion_analyser.visualization.syri_integration import SyRIIntegrator
            
            syri_dir = phylo_output_dir / 'syri_plots'
            syri_integrator = SyRIIntegrator(syri_dir)
            
            syri_plots_created = 0
            for pair_name, pair_results in all_results.items():
                if 'full_results' in pair_results:
                    try:
                        species1, species2 = pair_results['species_pair']
                        
                        print(f"  • Creating SyRI plot for {pair_name}...")
                        
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
                        
                        syri_plots_created += 1
                        print(f"    ✅ {pair_name}: SyRI plots created")
                        
                    except Exception as e:
                        print(f"    ❌ {pair_name}: {e}")
            
            print(f"  ✅ Created SyRI plots for {syri_plots_created} species pairs")
            
        except Exception as e:
            print(f"  ❌ SyRI integration failed: {e}")
            import traceback
            traceback.print_exc()

        # CREATE PUBLICATION PLOTS
        print(f"\n🎨 Creating publication-quality plots...")

        try:
            from genome_inversion_analyser.config import PUBLICATION_CONFIG
            from genome_inversion_analyser.visualization.publication_plots import create_publication_plots
            
            # Update config with publication settings
            pub_config = {
                'publication_config': PUBLICATION_CONFIG
            }
            
            # Create publication plots
            pub_results = create_publication_plots(
                all_results,  # Pass all pairwise results
                None,  # Not a single ortholog_df
                None,  # Not a single inversion_df
                registry,
                pub_config,
                phylo_output_dir,
                species_stats
            )
            
            if pub_results:
                print(f"  ✅ Publication plots created:")
                for plot_type, plots in pub_results.items():
                    if isinstance(plots, dict):
                        print(f"    • {plot_type}: {len(plots)} files")
                        for name, path in plots.items():
                            print(f"      ✓ {name}: {path}")
                    else:
                        print(f"    • {plot_type}: {plots}")
            else:
                print(f"  ⚠️ Publication plots returned empty")
                
        except Exception as e:
            print(f"  ❌ Publication plots failed: {e}")
            import traceback
            traceback.print_exc()
        
        # CREATE PHYLOGENETIC TREE PLOTS
        print(f"\n🌳 Creating phylogenetic tree plots...")
        
        if distance_matrices and trees:
            from scipy.cluster.hierarchy import dendrogram
            
            tree_dir = phylo_output_dir / 'tree_plots'
            tree_dir.mkdir(exist_ok=True)
            
            # Create tree for each distance method
            species_names = list(species_stats.keys())
            
            for method_name, dist_matrix in distance_matrices.items():
                tree_key = f"{method_name}_ward"
                if tree_key in trees:
                    tree_data = trees[tree_key]
                    
                    plt.figure(figsize=(12, 8))
                    dendrogram(tree_data['linkage_matrix'], 
                              labels=species_names,
                              orientation='top')
                    
                    plt.title(f'Bibionidae Phylogeny ({method_name})', fontsize=14, fontweight='bold')
                    plt.xlabel('Species')
                    plt.ylabel('Distance')
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    
                    tree_file = tree_dir / f'tree_{method_name}.png'
                    plt.savefig(tree_file, dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    print(f"  • Tree saved: {tree_file}")
        
        print(f"\n🎉 COMPLETE ANALYSIS FINISHED!")
        print("=" * 60)
        
        print(f"\n📊 FINAL RESULTS:")
        print(f"  • Species analyzed: {len(species_stats)}")
        print(f"  • Successful pairwise comparisons: {successful_comparisons}")
        print(f"  • Distance matrices computed: {len(distance_matrices)}")
        print(f"  • Phylogenetic trees inferred: {len(trees)}")
        print(f"  • SyRI plots created: {syri_plots_created}")
        
        print(f"\n📁 OUTPUT LOCATIONS:")
        print(f"  • Pairwise results: pairwise_results/")
        print(f"  • Phylogenetic analysis: {phylo_output_dir}/")
        print(f"  • SyRI plots: {phylo_output_dir}/syri_plots/")
        print(f"  • Tree plots: {phylo_output_dir}/tree_plots/")
        
        print(f"\n🎯 SUCCESS! Complete Bibionidae analysis with SyRI plots and phylogenetic trees!")
        
        return True
        
    except Exception as e:
        print(f"❌ Phylogenetic analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🧬 FIXED Multi-Species Analysis")
    print("Based on your working code + SyRI integration")
    print("=" * 70)
    
    success = run_multi_species_analysis()
    
    if success:
        print(f"\n🚀 COMPLETE SUCCESS!")
        print(f"Your Bibionidae analysis with SyRI plots is ready!")
    else:
        print(f"\n❌ Analysis incomplete. Check error messages above.")