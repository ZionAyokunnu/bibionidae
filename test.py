#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from genome_inversion_analyser.visualization.publication_plots import PublicationPlotGenerator
from genome_inversion_analyser.registry import FileRegistry

print("📊 Testing synteny plot fix...")

# Load real data
ortholog_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/ortholog_pairs.csv')
inversion_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/inversion_events.csv')

print(f"Loaded: {len(ortholog_df)} orthologs, {len(inversion_df)} inversions")

# Test config (disable external tool)
test_config = {
    'publication_config': {
        'synteny_visualization': {'enabled': True},
        'external_tools': {'synteny_plotter': None}
    }
}

# Test
output_dir = Path("test_synteny_fix")
registry = FileRegistry(output_dir, project_name="synteny_test")
plot_generator = PublicationPlotGenerator(registry, test_config)

try:
    plot_file = plot_generator._create_single_fallback_plot(
        ortholog_df, inversion_df, 'Dioctria_linearis', 'Dioctria_rufipes', output_dir
    )
    
    if plot_file and plot_file.exists():
        print(f"✅ Synteny plot created: {plot_file}")
    else:
        print("❌ Synteny plot creation failed")
        
except Exception as e:
    print(f"❌ Synteny test failed: {e}")
    import traceback
    traceback.print_exc()