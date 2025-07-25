#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from genome_inversion_analyser.visualization.publication_plots import PublicationPlotGenerator
from genome_inversion_analyser.registry import FileRegistry

print("📊 Testing fixed publication synteny plot...")

# Load real data
ortholog_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/ortholog_pairs.csv')
inversion_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/inversion_events.csv')

print(f"Loaded: {len(ortholog_df)} orthologs, {len(inversion_df)} inversions")

# Test config
test_config = {
    'publication_config': {
        'synteny_visualization': {'enabled': True},
        'external_tools': {'synteny_plotter': None}
    }
}

# Test
output_dir = Path("test_publication_synteny")
output_dir.mkdir(exist_ok=True)

registry = FileRegistry(output_dir, project_name="synteny_test")
plot_generator = PublicationPlotGenerator(registry, test_config)

try:
    print("Calling _create_single_fallback_plot...")
    plot_file = plot_generator._create_single_fallback_plot(
        ortholog_df, inversion_df, 'Dioctria_linearis', 'Dioctria_rufipes', output_dir
    )
    
    if plot_file and plot_file.exists():
        print(f"✅ Publication synteny plot created: {plot_file}")
        print(f"File size: {plot_file.stat().st_size} bytes")
    else:
        print("❌ Publication synteny plot creation failed - no file returned")
        
except Exception as e:
    print(f"❌ Publication synteny test failed: {e}")
    import traceback
    traceback.print_exc()