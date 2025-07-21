# # Debug version - Enhanced Multi-Species Synteny Analysis
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.patches import Rectangle
# import os
# import glob

# def create_debug_synteny_analysis():
#     """Debug version to see what's happening with connections"""
#     print("=== DEBUGGING SYNTENY CONNECTIONS ===\n")
    
#     # Load data (simplified)
#     all_species_data = {}
#     busco_files = glob.glob('busco-data/*.tsv')
    
#     for file_path in sorted(busco_files)[:8]:  # First 8 species
#         species_name = os.path.basename(file_path).replace('.tsv', '')
        
#         with open(file_path, 'r') as f:
#             lines = f.readlines()
        
#         header_idx = None
#         for i, line in enumerate(lines):
#             if line.startswith('# Busco id'):
#                 header_idx = i
#                 break
        
#         if header_idx is None:
#             continue
            
#         data_lines = []
#         for line in lines[header_idx + 1:]:
#             if line.startswith('#') or line.strip() == '':
#                 continue
#             parts = line.strip().split('\t')
#             if len(parts) >= 8:
#                 data_lines.append(parts)
        
#         df = pd.DataFrame(data_lines, columns=[
#             'busco_id', 'status', 'sequence', 'start', 'end', 
#             'strand', 'score', 'length'
#         ])
        
#         df = df[df['status'] == 'Complete'].copy()
#         df['start'] = pd.to_numeric(df['start'], errors='coerce')
#         df['end'] = pd.to_numeric(df['end'], errors='coerce')
#         df = df.dropna(subset=['start', 'end'])
        
#         all_species_data[species_name] = df
#         print(f"Loaded {species_name}: {len(df)} genes")
    
#     # Find common genes
#     species_names = list(all_species_data.keys())
#     common_busco_ids = set(all_species_data[species_names[0]]['busco_id'])
    
#     for species_name in species_names[1:]:
#         common_busco_ids = common_busco_ids.intersection(
#             set(all_species_data[species_name]['busco_id'])
#         )
    
#     print(f"\nCommon genes across all species: {len(common_busco_ids)}")
    
#     # Create common datasets
#     common_datasets = {}
#     for species_name in species_names:
#         common_df = all_species_data[species_name][
#             all_species_data[species_name]['busco_id'].isin(common_busco_ids)
#         ].copy()
#         common_datasets[species_name] = common_df
#         print(f"{species_name}: {len(common_df)} common genes")
    
#     # DEBUG: Check pairwise synteny between first two species
#     species1, species2 = species_names[0], species_names[1]
#     print(f"\n=== DEBUGGING {species1} vs {species2} ===")
    
#     # Find common genes between these two
#     common_pair = common_datasets[species1].merge(
#         common_datasets[species2], on='busco_id', suffixes=('_1', '_2')
#     )
#     print(f"Common genes between pair: {len(common_pair)}")
    
#     # Group by chromosome pairs
#     chr_pairs = common_pair.groupby(['sequence_1', 'sequence_2']).size().reset_index()
#     chr_pairs.columns = ['chr1', 'chr2', 'gene_count']
#     chr_pairs = chr_pairs.sort_values('gene_count', ascending=False)
    
#     print(f"\nTop chromosome pairs:")
#     print(chr_pairs.head(10))
    
#     # Find largest synteny blocks
#     largest_pair = chr_pairs.iloc[0]
#     print(f"\nLargest synteny block: {largest_pair['chr1']} -> {largest_pair['chr2']} ({largest_pair['gene_count']} genes)")
    
#     # Create visualization with GUARANTEED connections
#     fig, ax = plt.subplots(figsize=(16, 10))
    
#     # Draw chromosomes for first 6 species
#     key_species = species_names[:6]
    
#     for i, species in enumerate(key_species):
#         y_pos = len(key_species) - i - 1
        
#         # Get top 4 chromosomes
#         chroms = common_datasets[species].groupby('sequence').agg({
#             'start': 'min', 'end': 'max', 'busco_id': 'count'
#         }).reset_index()
#         chroms.columns = ['sequence', 'chr_start', 'chr_end', 'gene_count']
#         chroms = chroms.sort_values('gene_count', ascending=False).head(4)
        
#         # Layout
#         total_genes = chroms['gene_count'].sum()
#         chroms['prop_length'] = chroms['gene_count'] / total_genes * 300
#         chroms['cum_start'] = chroms['prop_length'].cumsum() - chroms['prop_length']
#         chroms['cum_end'] = chroms['cum_start'] + chroms['prop_length']
        
#         # Draw chromosomes
#         for _, chr_data in chroms.iterrows():
#             rect = Rectangle((chr_data['cum_start'], y_pos - 0.15), 
#                            chr_data['prop_length'], 0.3,
#                            facecolor='lightblue' if i % 2 == 0 else 'lightcoral', 
#                            edgecolor='black', linewidth=1.5, alpha=0.7)
#             ax.add_patch(rect)
            
#             # Labels
#             ax.text(chr_data['cum_start'] + chr_data['prop_length']/2, y_pos + 0.35,
#                    f'Chr{chroms.index[chroms.index == chr_data.name][0] + 1}', 
#                    ha='center', va='center', fontsize=9, fontweight='bold')
        
#         # Species label
#         ax.text(-30, y_pos, species.replace('_', ' '), ha='right', va='center',
#                fontsize=11, fontstyle='italic', fontweight='bold')
    
#     # FORCE draw connections between top chromosome pairs
#     print(f"\nForcing connections for top {min(5, len(chr_pairs))} chromosome pairs...")
    
#     reference_species = key_species[0]
#     ref_layout = common_datasets[reference_species].groupby('sequence').agg({
#         'start': 'min', 'end': 'max', 'busco_id': 'count'
#     }).reset_index()
#     ref_layout.columns = ['sequence', 'chr_start', 'chr_end', 'gene_count']
#     ref_layout = ref_layout.sort_values('gene_count', ascending=False).head(4)
#     total_genes = ref_layout['gene_count'].sum()
#     ref_layout['prop_length'] = ref_layout['gene_count'] / total_genes * 300
#     ref_layout['cum_start'] = ref_layout['prop_length'].cumsum() - ref_layout['prop_length']
    
#     # Draw connections for each species vs reference
#     colors = ['red', 'blue', 'green', 'purple', 'orange']
    
#     for target_idx, target_species in enumerate(key_species[1:], 1):
#         target_layout = common_datasets[target_species].groupby('sequence').agg({
#             'start': 'min', 'end': 'max', 'busco_id': 'count'
#         }).reset_index()
#         target_layout.columns = ['sequence', 'chr_start', 'chr_end', 'gene_count']
#         target_layout = target_layout.sort_values('gene_count', ascending=False).head(4)
#         total_genes = target_layout['gene_count'].sum()
#         target_layout['prop_length'] = target_layout['gene_count'] / total_genes * 300
#         target_layout['cum_start'] = target_layout['prop_length'].cumsum() - target_layout['prop_length']
        
#         # Find pairwise common genes
#         pair_common = common_datasets[reference_species].merge(
#             common_datasets[target_species], on='busco_id', suffixes=('_ref', '_target')
#         )
        
#         # Get top chromosome pair
#         top_pairs = pair_common.groupby(['sequence_ref', 'sequence_target']).size().reset_index()
#         top_pairs.columns = ['chr_ref', 'chr_target', 'count']
#         top_pairs = top_pairs.sort_values('count', ascending=False).head(3)
        
#         print(f"  {target_species}: Drawing {len(top_pairs)} connections")
        
#         ref_y = len(key_species) - 1
#         target_y = len(key_species) - target_idx - 1
        
#         for _, pair in top_pairs.iterrows():
#             # Find positions
#             ref_chr_data = ref_layout[ref_layout['sequence'] == pair['chr_ref']]
#             target_chr_data = target_layout[target_layout['sequence'] == pair['chr_target']]
            
#             if len(ref_chr_data) > 0 and len(target_chr_data) > 0:
#                 ref_pos = ref_chr_data['cum_start'].iloc[0] + ref_chr_data['prop_length'].iloc[0] / 2
#                 target_pos = target_chr_data['cum_start'].iloc[0] + target_chr_data['prop_length'].iloc[0] / 2
                
#                 # Draw thick connection
#                 ax.plot([ref_pos, target_pos], [ref_y, target_y], 
#                        color=colors[target_idx % len(colors)], 
#                        linewidth=max(2, pair['count'] / 50), 
#                        alpha=0.7)
                
#                 print(f"    Drew connection: {pair['chr_ref']} -> {pair['chr_target']} ({pair['count']} genes)")
    
#     # Formatting
#     ax.set_xlim(-80, 350)
#     ax.set_ylim(-0.5, len(key_species) - 0.5)
#     ax.axis('off')
#     ax.set_title(f'DEBUG: Forced Synteny Connections\n'
#                 f'(Line thickness ∝ gene count)',
#                 fontsize=16, fontweight='bold', pad=20)
    
#     plt.tight_layout()
#     plt.savefig('debug_synteny_connections.png', dpi=300, bbox_inches='tight')
#     plt.show()
    
#     print(f"\n✓ Debug plot saved: debug_synteny_connections.png")
#     print(f"✓ This should show connections between chromosomes")

# # Run debug version
# create_debug_synteny_analysis()


# WORKING Enhanced Multi-Species Synteny Analysis
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import os
import glob

def create_working_enhanced_synteny():
    """Working enhanced version based on debug insights"""
    print("=== WORKING ENHANCED SYNTENY ANALYSIS ===\n")
    
    # Load data (same as debug)
    all_species_data = {}
    busco_files = glob.glob('busco-data/*.tsv')
    
    for file_path in sorted(busco_files)[:8]:  # First 8 species
        species_name = os.path.basename(file_path).replace('.tsv', '')
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        header_idx = None
        for i, line in enumerate(lines):
            if line.startswith('# Busco id'):
                header_idx = i
                break
        
        if header_idx is None:
            continue
            
        data_lines = []
        for line in lines[header_idx + 1:]:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                data_lines.append(parts)
        
        df = pd.DataFrame(data_lines, columns=[
            'busco_id', 'status', 'sequence', 'start', 'end', 
            'strand', 'score', 'length'
        ])
        
        df = df[df['status'] == 'Complete'].copy()
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        df = df.dropna(subset=['start', 'end'])
        
        all_species_data[species_name] = df
        print(f"Loaded {species_name}: {len(df)} genes")
    
    # Find common genes
    species_names = list(all_species_data.keys())
    common_busco_ids = set(all_species_data[species_names[0]]['busco_id'])
    
    for species_name in species_names[1:]:
        common_busco_ids = common_busco_ids.intersection(
            set(all_species_data[species_name]['busco_id'])
        )
    
    print(f"\nCommon genes across all species: {len(common_busco_ids)}")
    
    # Create common datasets
    common_datasets = {}
    for species_name in species_names:
        common_df = all_species_data[species_name][
            all_species_data[species_name]['busco_id'].isin(common_busco_ids)
        ].copy()
        common_datasets[species_name] = common_df
    
    # Create enhanced visualization
    fig, ax = plt.subplots(figsize=(20, len(species_names) * 1.8))
    
    # Family colors
    family_colors = {
        'Aedes_aegypti': '#FF6B6B', 'Anopheles_gambia': '#FF6B6B', 'Culex_pipiens': '#FF6B6B',  # Culicidae
        'Bibio_marci': '#4ECDC4', 'Dilophus_febrilis': '#4ECDC4',  # Bibionidae
        'Drosophila_melanogaster': '#45B7D1', 'Drosophila_simulans': '#45B7D1',  # Drosophilidae
        'Bactrocera_dorsalis': '#FFA07A',  # Tephritidae
        'Calliphora_vicina': '#98D8C8',  # Calliphoridae  
        'Chironomus_riparius': '#F7DC6F'  # Chironomidae
    }
    
    # Draw chromosomes and create layouts
    layouts = {}
    for i, species in enumerate(species_names):
        y_pos = len(species_names) - i - 1
        
        # Get top chromosomes
        chroms = common_datasets[species].groupby('sequence').agg({
            'start': 'min', 'end': 'max', 'busco_id': 'count'
        }).reset_index()
        chroms.columns = ['sequence', 'chr_start', 'chr_end', 'gene_count']
        chroms = chroms.sort_values('gene_count', ascending=False).head(6)
        
        # Calculate proportional layout
        total_genes = chroms['gene_count'].sum()
        chroms['prop_length'] = chroms['gene_count'] / total_genes * 400
        chroms['cum_start'] = chroms['prop_length'].cumsum() - chroms['prop_length']
        chroms['cum_end'] = chroms['cum_start'] + chroms['prop_length']
        chroms['chr_name'] = [f'Chr{j+1}' for j in range(len(chroms))]
        
        layouts[species] = chroms
        
        # Draw chromosomes
        color = family_colors.get(species, '#CCCCCC')
        for _, chr_data in chroms.iterrows():
            rect = Rectangle((chr_data['cum_start'], y_pos - 0.12), 
                           chr_data['prop_length'], 0.24,
                           facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.8)
            ax.add_patch(rect)
            
            # Chromosome labels
            ax.text(chr_data['cum_start'] + chr_data['prop_length']/2, y_pos + 0.3,
                   chr_data['chr_name'], ha='center', va='center', 
                   fontsize=9, fontweight='bold')
        
        # Species labels
        ax.text(-40, y_pos, species.replace('_', ' '), ha='right', va='center',
               fontsize=11, fontstyle='italic', fontweight='bold')
    
    # Create ENHANCED SYNTENY BLOCKS with confidence scoring
    print("\nCreating enhanced synteny blocks...")
    
    reference_species = species_names[0]  # Aedes as reference
    ref_y = len(species_names) - 1
    
    # For each target species, create synteny blocks
    all_connections = []
    
    for target_idx, target_species in enumerate(species_names[1:], 1):
        target_y = len(species_names) - target_idx - 1
        
        # Find pairwise synteny
        pair_common = common_datasets[reference_species].merge(
            common_datasets[target_species], on='busco_id', suffixes=('_ref', '_target')
        )
        
        # Create synteny blocks (chromosome pairs with >50 genes)
        synteny_blocks = pair_common.groupby(['sequence_ref', 'sequence_target']).agg({
            'busco_id': 'count',
            'start_ref': ['min', 'max'],
            'start_target': ['min', 'max']
        }).reset_index()
        
        # Flatten column names
        synteny_blocks.columns = ['chr_ref', 'chr_target', 'gene_count', 
                                'ref_start', 'ref_end', 'target_start', 'target_end']
        
        # Filter for significant blocks (>50 genes)
        synteny_blocks = synteny_blocks[synteny_blocks['gene_count'] >= 50]
        
        # Calculate confidence (gene_count / max_possible)
        max_genes = synteny_blocks['gene_count'].max()
        synteny_blocks['confidence'] = synteny_blocks['gene_count'] / max_genes
        
        print(f"  {target_species}: {len(synteny_blocks)} synteny blocks")
        
        # Draw enhanced ribbons
        for _, block in synteny_blocks.iterrows():
            # Find chromosome positions
            ref_chr = layouts[reference_species][
                layouts[reference_species]['sequence'] == block['chr_ref']
            ]
            target_chr = layouts[target_species][
                layouts[target_species]['sequence'] == block['chr_target']
            ]
            
            if len(ref_chr) > 0 and len(target_chr) > 0:
                # Calculate positions within chromosomes
                ref_rel_start = (block['ref_start'] - ref_chr['chr_start'].iloc[0]) / (ref_chr['chr_end'].iloc[0] - ref_chr['chr_start'].iloc[0])
                ref_rel_end = (block['ref_end'] - ref_chr['chr_start'].iloc[0]) / (ref_chr['chr_end'].iloc[0] - ref_chr['chr_start'].iloc[0])
                
                target_rel_start = (block['target_start'] - target_chr['chr_start'].iloc[0]) / (target_chr['chr_end'].iloc[0] - target_chr['chr_start'].iloc[0])
                target_rel_end = (block['target_end'] - target_chr['chr_start'].iloc[0]) / (target_chr['chr_end'].iloc[0] - target_chr['chr_start'].iloc[0])
                
                # Calculate plot positions
                ref_start_pos = ref_chr['cum_start'].iloc[0] + ref_rel_start * ref_chr['prop_length'].iloc[0]
                ref_end_pos = ref_chr['cum_start'].iloc[0] + ref_rel_end * ref_chr['prop_length'].iloc[0]
                target_start_pos = target_chr['cum_start'].iloc[0] + target_rel_start * target_chr['prop_length'].iloc[0]
                target_end_pos = target_chr['cum_start'].iloc[0] + target_rel_end * target_chr['prop_length'].iloc[0]
                
                # Create thick curved ribbon
                n_points = 30
                t = np.linspace(0, 1, n_points)
                
                # Top edge (ref_start to target_start)
                top_x = ref_start_pos * (1-t) + target_start_pos * t
                top_y = ref_y * (1-t)**2 + (ref_y + target_y)/2 * 2*(1-t)*t + target_y * t**2
                
                # Bottom edge (ref_end to target_end)
                bottom_x = ref_end_pos * (1-t) + target_end_pos * t
                bottom_y = ref_y * (1-t)**2 + (ref_y + target_y)/2 * 2*(1-t)*t + target_y * t**2
                
                # Create ribbon polygon
                ribbon_x = np.concatenate([top_x, bottom_x[::-1]])
                ribbon_y = np.concatenate([top_y, bottom_y[::-1]])
                
                # Color and transparency based on confidence
                color = plt.cm.viridis(block['confidence'])
                alpha = 0.4 + (block['confidence'] * 0.5)
                
                # Draw ribbon
                ax.fill(ribbon_x, ribbon_y, color=color, alpha=alpha, 
                       edgecolor='black', linewidth=0.5)
                
                all_connections.append({
                    'species': target_species,
                    'genes': block['gene_count'],
                    'confidence': block['confidence']
                })
    
    # Add colorbar for confidence
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Synteny Block Confidence', fontsize=12)
    
    # Formatting
    ax.set_xlim(-100, 450)
    ax.set_ylim(-0.5, len(species_names) - 0.5)
    ax.axis('off')
    ax.set_title(f'Enhanced Multi-Species Synteny Analysis\n'
                f'Synteny Blocks with Confidence Scoring ({len(common_busco_ids)} common genes)', 
                fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig('working_enhanced_synteny.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Summary statistics
    print(f"\n=== ENHANCED SYNTENY SUMMARY ===")
    print(f"Total connections drawn: {len(all_connections)}")
    print(f"High confidence blocks (>0.8): {sum(1 for c in all_connections if c['confidence'] > 0.8)}")
    print(f"Medium confidence blocks (0.5-0.8): {sum(1 for c in all_connections if 0.5 <= c['confidence'] <= 0.8)}")
    print(f"Average genes per block: {np.mean([c['genes'] for c in all_connections]):.1f}")
    
    print(f"\n✓ Enhanced synteny plot saved: working_enhanced_synteny.png")

# Run the working enhanced version
create_working_enhanced_synteny()