# #!/usr/bin/env python3

# import sys
# import pandas as pd
# from pathlib import Path

# print("🔍 Debugging column issue...")

# # Load the data
# ortholog_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/ortholog_pairs.csv')

# print("Original columns:", list(ortholog_df.columns))

# # Test renaming
# ortholog_df_mapped = ortholog_df.rename(columns={
#     'first_chr': 'first_chromosome',
#     'second_chr': 'second_chromosome'
# })

# print("After renaming:", list(ortholog_df_mapped.columns))

# # Check if renaming worked
# if 'first_chromosome' in ortholog_df_mapped.columns:
#     print("✅ Renaming worked!")
#     print("Sample data:")
#     print(ortholog_df_mapped[['first_chromosome', 'second_chromosome', 'first_start', 'second_start']].head(2))
# else:
#     print("❌ Renaming failed!")
    
# # Test accessing the columns
# try:
#     chr1_list = sorted(ortholog_df_mapped['first_chromosome'].unique())
#     chr2_list = sorted(ortholog_df_mapped['second_chromosome'].unique())
#     print(f"✅ Chromosomes found: {len(chr1_list)} and {len(chr2_list)}")
#     print(f"Chr1 sample: {chr1_list[:3]}")
#     print(f"Chr2 sample: {chr2_list[:3]}")
# except Exception as e:
#     print(f"❌ Column access failed: {e}")




    #!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Simple direct test
print("📊 Direct synteny plot test...")

# Load data
ortholog_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/ortholog_pairs.csv')
print(f"Loaded {len(ortholog_df)} orthologs")

# Test direct plotting
output_dir = Path("test_direct_synteny")
output_dir.mkdir(exist_ok=True)

try:
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get chromosomes
    chr1_list = sorted(ortholog_df['first_chr'].unique())
    chr2_list = sorted(ortholog_df['second_chr'].unique())
    
    print(f"Chromosomes: {len(chr1_list)} vs {len(chr2_list)}")
    
    # Simple plot
    for i, (_, row) in enumerate(ortholog_df.head(100).iterrows()):  # Just first 100
        y1 = 0.7
        y2 = 0.3
        x1 = 0.2 + 0.6 * (i / 100)
        x2 = 0.2 + 0.6 * (i / 100)
        
        color = plt.cm.viridis(row['similarity'])
        ax.plot([x1, x2], [y1, y2], color=color, alpha=0.6, linewidth=1)
    
    ax.set_title('Simple Synteny Test')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    plot_file = output_dir / 'simple_synteny_test.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Simple test successful: {plot_file}")
    
except Exception as e:
    print(f"❌ Simple test failed: {e}")
    import traceback
    traceback.print_exc()