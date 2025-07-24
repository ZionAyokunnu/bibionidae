#!/usr/bin/env python3
"""
Multi-Species Wrapper for Existing Pairwise Analysis
Uses your proven main.py logic in loops for phylogenetic analysis
"""

import sys
import json
import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations
import shutil

# Import your working modules exactly as they are
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
    """Create config for one pairwise comparison using your exact config structure"""
    
    # Start with your exact ENHANCED_HYBRID_CONFIG
    pairwise_config = ENHANCED_HYBRID_CONFIG.copy()
    
    # Update paths for this specific pair
    pair_name = f"{species1['name']}_vs_{species2['name']}"
    
    pairwise_config.update({
        # Input files for this pair
        'first_fasta_path': species1['fasta'],
        'second_fasta_path': species2['fasta'], 
        'first_busco_path': species1['busco'],
        'second_busco_path': species2['busco'],
        
        # Create unique output paths for this pair
        'base_output_dir': f'pairwise_results/{pair_name}',
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
        print(f"     FASTA: {species['fasta']}")
        print(f"     BUSCO: {species['busco']}")
    
    # Verify all files exist first
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
    
    # Create output directory
    results_dir = Path("pairwise_results")
    results_dir.mkdir(exist_ok=True)
    
    # Generate all pairwise combinations
    species_pairs = list(combinations(species_data, 2))
    total_comparisons = len(species_pairs)
    
    print(f"\n🔄 Running {total_comparisons} pairwise comparisons:")
    
    # Store results for phylogenetic analysis
    all_results = {}
    species_stats = {}
    
    for i, (species1, species2) in enumerate(species_pairs, 1):
        pair_name = f"{species1['name']}_vs_{species2['name']}"
        print(f"\n📊 [{i}/{total_comparisons}] {pair_name}")
        
        try:
            # Create pairwise config using your exact structure
            pairwise_config = create_pairwise_config(species1, species2)
            
            print(f"  • Running your pairwise analysis...")
            
            # Call your exact main function
            results = run_complete_enhanced_analysis_with_registry(pairwise_config)
            
            if results is None:
                print(f"  ❌ Analysis returned None")
                continue
            
            # Extract key information for phylogenetic analysis
            inversion_count = len(results['inversion_df'])
            synteny_blocks = len(results['synteny_df'])
            ortholog_count = len(results['ortholog_df'])
            avg_similarity = results['ortholog_df']['similarity'].mean() if len(results['ortholog_df']) > 0 else 0
            
            print(f"  • Results: {ortholog_count} orthologs, {inversion_count} inversions, {synteny_blocks} synteny blocks")
            print(f"  • Average similarity: {avg_similarity:.3f}")
            
            # Store results
            all_results[pair_name] = {
                'species_pair': (species1['name'], species2['name']),
                'inversion_df': results['inversion_df'],
                'ortholog_df': results['ortholog_df'],
                'synteny_df': results['synteny_df'],
                'inversion_count': inversion_count,
                'avg_similarity': avg_similarity,
                'full_results': results
            }
            
            # Store individual species stats (for first occurrence)
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
    
    # Now build phylogenetic analysis
    print(f"\n🌳 Building phylogenetic trees...")
    
    try:
        # Set up phylogenetic analysis
        phylo_output_dir = Path("bibionidae_phylogenetic_analysis")
        registry = FileRegistry(phylo_output_dir, project_name="bibionidae_phylogeny")
        phylo_analyzer = PhylogeneticIntegrator(registry)
        
        # Add each species to phylogenetic analyzer
        for species_name, stats in species_stats.items():
            
            # Collect all inversions involving this species
            species_inversions = []
            species_orthologs = []
            
            for pair_name, pair_results in all_results.items():
                if 'error' not in pair_results and species_name in pair_results['species_pair']:
                    # Add inversions
                    inv_df = pair_results['inversion_df'].copy()
                    if not inv_df.empty:
                        inv_df['comparison'] = pair_name
                        species_inversions.extend(inv_df.to_dict('records'))
                    
                    # Add orthologs for better analysis
                    orth_df = pair_results['ortholog_df'].copy()
                    if not orth_df.empty:
                        orth_df['comparison'] = pair_name
                        species_orthologs.extend(orth_df.to_dict('records'))
            
            # Convert to DataFrame
            species_inversion_df = pd.DataFrame(species_inversions) if species_inversions else pd.DataFrame()
            species_ortholog_df = pd.DataFrame(species_orthologs) if species_orthologs else pd.DataFrame()
            
            # Calculate contextual metrics
            genome_size = stats['quality']['metrics'].get('total_length', 1000000)
            contextual_metrics = {
                'rate_metrics': {
                    'inversions_per_mb': len(species_inversions) / (genome_size / 1_000_000)
                },
                'species_stats': {
                    'genome_size': genome_size,
                    'quality_class': stats['quality']['quality_class'],
                    'n_orthologs': len(species_orthologs)
                }
            }
            
            # Add to phylogenetic analysis
            phylo_analyzer.add_species_data(
                species_name,
                species_inversion_df,
                stats['quality']['metrics'],
                contextual_metrics
            )
            
            print(f"  • Added {species_name}: {len(species_inversions)} inversions, {contextual_metrics['rate_metrics']['inversions_per_mb']:.3f} inv/Mb")
        
        # Generate phylogenetic analysis
        print(f"  • Computing distance matrices...")
        distance_matrices = phylo_analyzer.compute_distance_matrices()
        
        print(f"  • Inferring phylogenetic trees...")
        trees = phylo_analyzer.infer_phylogenetic_trees()
        
        print(f"  • Analyzing phylogenetic signal...")
        signal_results = phylo_analyzer.analyze_phylogenetic_signal()
        
        print(f"  • Generating comprehensive summary...")
        phylo_summary = phylo_analyzer.generate_phylogenetic_summary()
        
        print(f"\n🎉 PHYLOGENETIC ANALYSIS COMPLETE!")
        print("=" * 60)
        
        # Print detailed results summary
        print(f"\n📊 COMPREHENSIVE RESULTS:")
        print(f"  • Species analyzed: {len(species_stats)}")
        print(f"  • Successful pairwise comparisons: {successful_comparisons}")
        print(f"  • Distance matrices computed: {len(distance_matrices)}")
        print(f"  • Phylogenetic trees inferred: {len(trees)}")
        
        print(f"\n🧬 SPECIES SUMMARY:")
        for species_name, stats in species_stats.items():
            quality = stats['quality']
            print(f"  • {species_name}:")
            print(f"    - Quality: {quality['quality_class']} (score: {quality['quality_score']:.3f})")
            print(f"    - Genome size: {quality['metrics']['total_length']:,} bp")
        
        print(f"\n🌳 PHYLOGENETIC TREES:")
        for tree_name in trees.keys():
            print(f"  • {tree_name}.json")
        
        print(f"\n📁 OUTPUT LOCATIONS:")
        print(f"  • Pairwise results: pairwise_results/")
        print(f"  • Phylogenetic analysis: {phylo_output_dir}/")
        print(f"  • Distance matrices: {phylo_output_dir}/exports/phylogenetic/")
        print(f"  • Registry and manifest: {registry.registry_file}")
        
        print(f"\n🎯 SUCCESS! Your Bibionidae phylogenetic analysis is complete!")
        print(f"📈 You now have evolutionary trees showing relationships between your 4 species!")
        
        return True
        
    except Exception as e:
        print(f"❌ Phylogenetic analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    print("🧬 Bibionidae Multi-Species Phylogenetic Analysis")
    print("Using your proven pairwise analysis pipeline")
    print("=" * 70)
    
    success = run_multi_species_analysis()
    
    if success:
        print(f"\n🚀 COMPLETE SUCCESS!")
        print(f"Your Bibionidae phylogenetic trees are ready for analysis!")
    else:
        print(f"\n❌ Analysis incomplete. Check error messages above.")