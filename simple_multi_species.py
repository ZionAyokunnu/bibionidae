#!/usr/bin/env python3
"""
SIMPLE Multi-Species Analysis - Step by Step Debug
"""

import sys
import json
import pandas as pd
from pathlib import Path
from itertools import combinations

# Add path and imports
sys.path.insert(0, str(Path(__file__).parent))

print("🧬 SIMPLE MULTI-SPECIES ANALYSIS")
print("=" * 50)

try:
    from genome_inversion_analyser.main import run_complete_enhanced_analysis_with_registry
    from genome_inversion_analyser.config import ENHANCED_HYBRID_CONFIG
    from genome_inversion_analyser.registry import FileRegistry
    from genome_inversion_analyser.phylogenetic import PhylogeneticIntegrator
    print("✅ All imports successful")
except Exception as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

# Load config
print("\n📋 Loading species configuration...")
with open('species_config.json', 'r') as f:
    config = json.load(f)

species_data = config['species_data']
print(f"Species found: {[s['name'] for s in species_data]}")

# Test ONE pairwise analysis first
print(f"\n🔬 Testing ONE pairwise analysis...")

# Pick first two species
species1, species2 = species_data[0], species_data[1]
print(f"Testing: {species1['name']} vs {species2['name']}")

# Create config for this pair
pair_config = ENHANCED_HYBRID_CONFIG.copy()
pair_config.update({
    'first_fasta_path': species1['fasta'],
    'second_fasta_path': species2['fasta'],
    'first_busco_path': species1['busco'],
    'second_busco_path': species2['busco'],
    'base_output_dir': 'test_single_pair'
})

print("Config created, running analysis...")

try:
    # Run ONE pairwise analysis
    results = run_complete_enhanced_analysis_with_registry(pair_config)
    
    if results is None:
        print("❌ Analysis returned None")
    else:
        print(f"✅ Analysis successful!")
        print(f"  • Orthologs: {len(results['ortholog_df'])}")
        print(f"  • Inversions: {len(results['inversion_df'])}")
        print(f"  • Synteny blocks: {len(results['synteny_df'])}")
        
        # Now test phylogenetic integration
        print(f"\n🌳 Testing phylogenetic integration...")
        
        # Create phylogenetic analyzer
        phylo_dir = Path("test_phylo")
        registry = FileRegistry(phylo_dir, project_name="test_phylo")
        phylo_analyzer = PhylogeneticIntegrator(registry)
        
        # Add this pair as fake "species" for testing
        phylo_analyzer.add_species_data(
            species1['name'],
            results['inversion_df'],
            results['first_quality']['metrics'],
            {'rate_metrics': {'inversions_per_mb': len(results['inversion_df'])}}
        )
        
        phylo_analyzer.add_species_data(
            species2['name'], 
            results['inversion_df'],
            results['second_quality']['metrics'],
            {'rate_metrics': {'inversions_per_mb': len(results['inversion_df'])}}
        )
        
        print("✅ Species added to phylogenetic analyzer")
        
        # Test distance matrices
        print("Computing distance matrices...")
        distance_matrices = phylo_analyzer.compute_distance_matrices(['jaccard'])
        print(f"✅ Distance matrices: {list(distance_matrices.keys())}")
        
        # Test SyRI integration
        print(f"\n📊 Testing SyRI integration...")
        from genome_inversion_analyser.visualization.syri_integration import SyRIIntegrator
        
        syri_dir = Path("test_syri")
        syri_integrator = SyRIIntegrator(syri_dir)
        
        # Create SyRI output
        syri_output = syri_integrator.create_syri_compatible_output(
            results['ortholog_df'],
            results['synteny_df'],
            results['inversion_df'],
            species1['name'],
            species2['name']
        )
        print(f"✅ SyRI output created: {syri_output}")
        
        # Create SyRI plot
        syri_plot = syri_integrator.create_syri_style_plot(
            results['ortholog_df'],
            results['synteny_df'],
            results['inversion_df'],
            species1['name'],
            species2['name']
        )
        print(f"✅ SyRI plot created: {syri_plot}")
        
        print(f"\n🎉 ALL COMPONENTS WORKING!")
        print(f"Ready to scale to all species pairs")


        
except Exception as e:
    print(f"❌ Analysis failed: {e}")
    import traceback
    traceback.print_exc()

print(f"\n📁 Check these directories for outputs:")
print(f"  • test_single_pair/")
print(f"  • test_phylo/")
print(f"  • test_syri/")