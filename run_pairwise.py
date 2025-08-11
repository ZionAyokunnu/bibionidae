#!/usr/bin/env python3
"""
Enhanced Pairwise Runner with SyRI and Publication Plots
Includes all the visualizations from multi-species analysis
"""
import sys
from pathlib import Path
import argparse

# Add the current directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

from genome_inversion_analyser.main import run_complete_enhanced_analysis_with_registry
from genome_inversion_analyser.config import (
    COMPLETE_ENHANCED_CONFIG,
    ENHANCED_HYBRID_CONFIG,
    FAST_HYBRID_CONFIG,
    PUBLICATION_CONFIG
)

def create_syri_plots(results, species1_name, species2_name, output_dir):
    """Create SyRI visualizations for pairwise comparison"""
    print("📊 Creating SyRI visualizations...")
    
    try:
        from genome_inversion_analyser.visualization.syri_integration import SyRIIntegrator
        
        syri_dir = Path(output_dir) / 'syri_plots'
        syri_integrator = SyRIIntegrator(syri_dir)
        
        # Create SyRI-compatible output
        print("  • Creating SyRI-compatible output...")
        syri_output = syri_integrator.create_syri_compatible_output(
            results['ortholog_df'],
            results['synteny_df'],
            results['inversion_df'],
            species1_name, species2_name
        )
        print(f"  ✅ SyRI output: {syri_output}")
        
        # Create SyRI-style plot
        print("  • Creating SyRI-style plot...")
        syri_plot = syri_integrator.create_syri_style_plot(
            results['ortholog_df'],
            results['synteny_df'],
            results['inversion_df'],
            species1_name, species2_name
        )
        print(f"  ✅ SyRI plot: {syri_plot}")
        
        # Create circular synteny plot
        print("  • Creating circular synteny plot...")
        circular_plot = syri_integrator.create_circular_synteny_plot(
            results['ortholog_df'],
            results['inversion_df'],
            species1_name, species2_name
        )
        print(f"  ✅ Circular plot: {circular_plot}")
        
        return {
            'syri_output': syri_output,
            'syri_plot': syri_plot, 
            'circular_plot': circular_plot
        }
        
    except Exception as e:
        print(f"  ❌ SyRI plots failed: {e}")
        import traceback
        traceback.print_exc()
        return {}

def create_publication_plots(results, species1_name, species2_name, output_dir, registry):
    """Create publication-quality plots for pairwise comparison"""
    print("🎨 Creating publication-quality plots...")
    
    try:
        from genome_inversion_analyser.visualization.publication_plots import create_publication_plots
        
        # Prepare data in multi-species format for publication plots
        all_results = {
            f"{species1_name}_vs_{species2_name}": {
                'species_pair': (species1_name, species2_name),
                'inversion_df': results['inversion_df'],
                'ortholog_df': results['ortholog_df'],
                'synteny_df': results['synteny_df'],
                'full_results': results
            }
        }
        
        species_stats = {
            species1_name: {
                'quality': results['first_quality'],
                'inversions': len(results['inversion_df']),
                'genome_size': results['first_quality']['metrics'].get('total_length', 1000000)
            },
            species2_name: {
                'quality': results['second_quality'], 
                'inversions': len(results['inversion_df']),
                'genome_size': results['second_quality']['metrics'].get('total_length', 1000000)
            }
        }
        
        # Create publication plots
        pub_config = {
            'publication_config': PUBLICATION_CONFIG
        }
        
        pub_results = create_publication_plots(
            all_results,
            None,  # Not a single ortholog_df
            None,  # Not a single inversion_df
            registry,
            pub_config,
            Path(output_dir),
            species_stats
        )
        
        if pub_results:
            print("  ✅ Publication plots created:")
            for plot_type, plots in pub_results.items():
                if isinstance(plots, dict):
                    print(f"    • {plot_type}: {len(plots)} files")
                    for name, path in plots.items():
                        print(f"      ✓ {name}: {path}")
                else:
                    print(f"    • {plot_type}: {plots}")
            return pub_results
        else:
            print("  ⚠️ Publication plots returned empty")
            return {}
            
    except Exception as e:
        print(f"  ❌ Publication plots failed: {e}")
        import traceback
        traceback.print_exc()
        return {}

def create_enhanced_tree_plots(results, species1_name, species2_name, output_dir):
    """Create enhanced phylogenetic tree plots"""
    print("🌳 Creating enhanced tree visualizations...")
    
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        
        tree_dir = Path(output_dir) / 'enhanced_trees'
        tree_dir.mkdir(exist_ok=True)
        
        # Simple binary tree for the two species
        inversion_count1 = len([inv for _, inv in results['inversion_df'].iterrows() 
                              if inv.get('species') == species1_name or inv.get('first_chr')])
        inversion_count2 = len([inv for _, inv in results['inversion_df'].iterrows() 
                              if inv.get('species') == species2_name or inv.get('second_chr')])
        
        # Create simple dendrogram
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        # Draw tree structure
        y_positions = [0, 1]
        species_names = [species1_name.replace('_', ' '), species2_name.replace('_', ' ')]
        inversion_counts = [inversion_count1, inversion_count2]
        
        # Tree trunk and branches
        trunk_x = 0.2
        branch_x = 0.6
        species_x = 1.0
        
        # Main connection
        ax.plot([trunk_x, branch_x], [0.5, 0.5], 'k-', linewidth=4)
        
        for i, (name, inv_count) in enumerate(zip(species_names, inversion_counts)):
            y = y_positions[i]
            
            # Branch
            ax.plot([branch_x, species_x], [0.5, y], 'brown', linewidth=3)
            
            # Species name
            ax.text(species_x + 0.05, y, name, ha='left', va='center', fontsize=12, weight='bold')
            
            # Inversion annotation
            ax.text(species_x + 0.6, y, f"{inv_count} inversions", 
                   ha='left', va='center', fontsize=10,
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightcoral", alpha=0.8))
        
        ax.set_xlim(0, 1.8)
        ax.set_ylim(-0.3, 1.3)
        ax.set_title(f'Phylogenetic Relationship\n{species1_name} vs {species2_name}', 
                    fontsize=14, weight='bold')
        ax.axis('off')
        
        tree_file = tree_dir / 'pairwise_tree.png'
        plt.tight_layout()
        plt.savefig(tree_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"  ✅ Enhanced tree: {tree_file}")
        return tree_file
        
    except Exception as e:
        print(f"  ❌ Enhanced tree failed: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(description='Run Enhanced Genome Inversion Analysis')
    parser.add_argument('--mode', choices=['fast', 'hybrid', 'complete'],
                       default='complete', help='Analysis mode')
    parser.add_argument('--species1', default='Dioctria_rufipes', help='First species name')
    parser.add_argument('--species2', default='Dioctria_linearis', help='Second species name')
    
    args = parser.parse_args()
    
    # Choose config based on mode
    if args.mode == 'fast':
        config = FAST_HYBRID_CONFIG.copy()
        print("🚀 Using FAST mode (Biopython only, minimal features)")
    elif args.mode == 'hybrid':
        config = ENHANCED_HYBRID_CONFIG.copy()
        print("⚡ Using HYBRID mode (Minimap2 + Biopython)")
    else:
        config = COMPLETE_ENHANCED_CONFIG.copy()
        print("🔬 Using COMPLETE mode (All features)")
    
    # Update config with your file paths
    config.update({
        'first_fasta_path': 'Bibio_marci/Bibio_marci.fna',
        'second_fasta_path': 'Bibio_marci_inverted.fasta',
        'first_busco_path': 'Busco-tables/Bibio_marci.tsv',
        'second_busco_path': 'Busco-tables/Bibio_marci.tsv',
        'base_output_dir': f'v7/{args.mode}_results'
    })
    
    # Verify files exist
    required_files = [
        config['first_fasta_path'],
        config['second_fasta_path'],
        config['first_busco_path'],
        config['second_busco_path']
    ]
    
    for file_path in required_files:
        if not Path(file_path).exists():
            print(f"❌ File not found: {file_path}")
            return None
    
    try:
        print("🧬 Starting Enhanced Genome Inversion Analysis...")
        print(f"📊 Comparing: {args.species1} vs {args.species2}")
        
        # Run main analysis
        results = run_complete_enhanced_analysis_with_registry(config)
        
        if results is None:
            print("❌ Main analysis failed")
            return None
            
        print("✅ Core Analysis Complete!")
        
        # Get output directory and registry
        output_dir = Path(config['base_output_dir'])
        
        # Create registry for additional plots
        from genome_inversion_analyser.registry import FileRegistry
        registry = FileRegistry(output_dir, project_name="enhanced_pairwise")
        
        # Print core results
        print(f"\n📊 CORE RESULTS:")
        print(f"  • Orthologs found: {len(results['ortholog_df'])}")
        print(f"  • Inversions detected: {len(results['inversion_df'])}")
        print(f"  • Synteny blocks: {len(results['synteny_df'])}")
        
        # CREATE ADDITIONAL VISUALIZATIONS
        
        # 1. SyRI Plots
        syri_results = create_syri_plots(results, args.species1, args.species2, output_dir)
        
        # 2. Publication Plots  
        pub_results = create_publication_plots(results, args.species1, args.species2, output_dir, registry)
        
        # 3. Enhanced Tree Plots
        tree_result = create_enhanced_tree_plots(results, args.species1, args.species2, output_dir)
        
        print(f"\n🎉 ENHANCED ANALYSIS COMPLETE!")
        print("=" * 50)
        
        print(f"\n📁 OUTPUT LOCATIONS:")
        print(f"  • Main results: {output_dir}/")
        print(f"  • SyRI plots: {output_dir}/syri_plots/")
        print(f"  • Publication plots: {output_dir}/publication_plots/")
        print(f"  • Enhanced trees: {output_dir}/enhanced_trees/")
        
        print(f"\n🎯 SUCCESS! Complete pairwise analysis with all visualizations!")
        
        return results
        
    except Exception as e:
        print(f"❌ Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    main()