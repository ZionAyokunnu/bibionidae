#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path

print("🔍 Debugging column issue...")

# Load the data
ortholog_df = pd.read_csv('pairwise_results/Dioctria_linearis_vs_Dioctria_rufipes/data/ortholog_pairs.csv')

print("Original columns:", list(ortholog_df.columns))

# Test renaming
ortholog_df_mapped = ortholog_df.rename(columns={
    'first_chr': 'first_chromosome',
    'second_chr': 'second_chromosome'
})

print("After renaming:", list(ortholog_df_mapped.columns))

# Check if renaming worked
if 'first_chromosome' in ortholog_df_mapped.columns:
    print("✅ Renaming worked!")
    print("Sample data:")
    print(ortholog_df_mapped[['first_chromosome', 'second_chromosome', 'first_start', 'second_start']].head(2))
else:
    print("❌ Renaming failed!")
    
# Test accessing the columns
try:
    chr1_list = sorted(ortholog_df_mapped['first_chromosome'].unique())
    chr2_list = sorted(ortholog_df_mapped['second_chromosome'].unique())
    print(f"✅ Chromosomes found: {len(chr1_list)} and {len(chr2_list)}")
    print(f"Chr1 sample: {chr1_list[:3]}")
    print(f"Chr2 sample: {chr2_list[:3]}")
except Exception as e:
    print(f"❌ Column access failed: {e}")