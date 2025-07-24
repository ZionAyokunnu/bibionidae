#!/usr/bin/env python3
"""
CLEAN Multi-Species Analysis - No Redundancies
Combines all working components into one streamlined script
"""

import sys
import json
import pandas as pd
from pathlib import Path
from itertools import combinations

sys.path.insert(0, str(Path(__file__).parent))

print("🧬 CLEAN MULTI-SPECIES ANALYSIS")
print("=" * 50)

# Import your proven working components
from genome_inversion_analyser.main import run_complete_enhanced_analysis_with_registry
from genome_inversion_analyser.config import ENHANCED_HYBRID_CONFIG
from genome_inversion_analyser.registry import FileRegistry
from genome_inversion_analyser.phylogenetic import PhylogeneticIntegrator
from genome_inversion_analyser.visualization.syri_integration import SyRIIntegrator

def main():
    """Single clean function - no redundancy"""
    
    # Load species config
    with open('species_config.json', 'r') as f:
        config = json.load(f)
    
    species_data = config['species_data']
    print(f"Species: {[s['name'] for s in species_data]}")
    
    # Clean output directories
    output_base = Path("clean_bibionidae_analysis")
    if output_base.exists():
        import shutil
        shutil.rmtree(output_base)
    
    # 1. RUN PAIRWISE ANALYSES
    print(f"\n🔄 PHASE 1: Pairwise Analyses")
    
    species_pairs = list(combinations(species_data, 2))
    all_results = {}
    species_stats = {}
    
    for i, (sp1, sp2) in enumerate(species_pairs, 1):
        pair_name = f"{sp1['name']}_vs_{sp2['name']}"
        print(f"  [{i}/{len(species_pairs)}] {pair_name}")
        
        try:
            # Create config
            pair_config = ENHANCED_HYBRID_CONFIG.copy()
            pair_config.update({
                'first_fasta_path': sp1['fasta'],
                'second_fasta_path': sp2['fasta'],
                'first_busco_path': sp1['busco'],
                'second_busco_path': sp2['busco'],
                'base_output_dir': str(output_base / 'pairwise' / pair_name),
                'enable_debug_output': False
            })
            
            # Run analysis
            results = run_complete_enhanced_analysis_with_registry(pair_config)
            
            if results:
                print(f"    ✅ {len(results['ortholog_df'])} orthologs, {len(results['inversion_df'])} inversions")
                
                all_results[pair_name] = {
                    'species_pair': (sp1['name'], sp2['name']),
                    'ortholog_df': results['ortholog_df'],
                    'inversion_df': results['inversion_df'],
                    'synteny_df': results['synteny_df'],
                    'results': results
                }
                
                # Store species stats
                for sp_info, quality_key in [(sp1, 'first_quality'), (sp2, 'second_quality')]:
                    if sp_info['name'] not in species_stats:
                        species_stats[sp_info['name']] = {
                            'quality': results[quality_key],
                            'metadata': sp_info
                        }
            else:
                print(f"    ❌ Failed")
                
        except Exception as e:
            print(f"    ❌ Error: {e}")
    
    print(f"\n✅ Completed {len([r for r in all_results.values() if 'results' in r])} pairwise analyses")
    
    # 2. PHYLOGENETIC ANALYSIS
    print(f"\n🌳 PHASE 2: Phylogenetic Analysis")
    
    phylo_dir = output_base / 'phylogenetic'
    registry = FileRegistry(phylo_dir, project_name="clean_phylogeny")
    phylo_analyzer = PhylogeneticIntegrator(registry, {'min_species': 2})
    
    # Add species data
    for species_name, stats in species_stats.items():
        all_inversions = []
        for pair_name, pair_data in all_results.items():
            if 'results' in pair_data and species_name in pair_data['species_pair']:
                inv_df = pair_data['inversion_df']
                if not inv_df.empty:
                    all_inversions.extend(inv_df.to_dict('records'))
        
        species_inv_df = pd.DataFrame(all_inversions) if all_inversions else pd.DataFrame()
        
        phylo_analyzer.add_species_data(
            species_name,
            species_inv_df,
            stats['quality']['metrics'],
            {'inversion_rate': len(all_inversions)}
        )
        
        print(f"  • {species_name}: {len(all_inversions)} inversions")
    
    # Compute phylogenetic analysis
    distance_matrices = phylo_analyzer.compute_distance_matrices(['jaccard', 'euclidean'])
    trees = phylo_analyzer.infer_phylogenetic_trees(['ward'])
    
    print(f"  ✅ Created {len(distance_matrices)} distance matrices, {len(trees)} trees")
    
    # 3. SYRI INTEGRATION
    print(f"\n📊 PHASE 3: SyRI Integration")
    
    syri_dir = output_base / 'syri_plots'
    syri_integrator = SyRIIntegrator(syri_dir)
    
    syri_plots_created = 0
    for pair_name, pair_data in all_results.items():
        if 'results' in pair_data:
            try:
                sp1, sp2 = pair_data['species_pair']
                
                # Create SyRI output
                syri_output = syri_integrator.create_syri_compatible_output(
                    pair_data['ortholog_df'],
                    pair_data['synteny_df'],
                    pair_data['inversion_df'],
                    sp1, sp2
                )
                
                # Create SyRI plot
                syri_plot = syri_integrator.create_syri_style_plot(
                    pair_data['ortholog_df'],
                    pair_data['synteny_df'],
                    pair_data['inversion_df'],
                    sp1, sp2
                )
                
                syri_plots_created += 1
                print(f"  • {pair_name}: SyRI plot created")
                
            except Exception as e:
                print(f"  ❌ {pair_name}: {e}")
    
    print(f"  ✅ Created {syri_plots_created} SyRI plots")
    
    # 4. TREE VISUALIZATION
    print(f"\n🌳 PHASE 4: Tree Visualization")
    
    if distance_matrices and trees:
        import matplotlib.pyplot as plt
        from scipy.cluster.hierarchy import dendrogram
        
        tree_dir = output_base / 'trees'
        tree_dir.mkdir(exist_ok=True)
        
        # Create tree for each distance method
        for method_name, dist_matrix in distance_matrices.items():
            tree_key = f"{method_name}_ward"
            if tree_key in trees:
                tree_data = trees[tree_key]
                
                plt.figure(figsize=(12, 8))
                dendrogram(tree_data['linkage_matrix'], 
                          labels=tree_data['species_names'],
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
    
    # 5. SUMMARY
    print(f"\n🎉 ANALYSIS COMPLETE!")
    print(f"📁 ALL OUTPUTS IN: {output_base}/")
    print(f"  • Pairwise results: {output_base}/pairwise/")
    print(f"  • Phylogenetic analysis: {output_base}/phylogenetic/")
    print(f"  • SyRI plots: {output_base}/syri_plots/")
    print(f"  • Tree plots: {output_base}/trees/")

if __name__ == "__main__":
    main()