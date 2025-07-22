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


# # WORKING Enhanced Multi-Species Synteny Analysis
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.patches import Rectangle
# import os
# import glob

# def create_working_enhanced_synteny():
#     """Working enhanced version based on debug insights"""
#     print("=== WORKING ENHANCED SYNTENY ANALYSIS ===\n")
    
#     # Load data (same as debug)
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
    
#     # Create enhanced visualization
#     fig, ax = plt.subplots(figsize=(20, len(species_names) * 1.8))
    
#     # Family colors
#     family_colors = {
#         'Aedes_aegypti': '#FF6B6B', 'Anopheles_gambia': '#FF6B6B', 'Culex_pipiens': '#FF6B6B',  # Culicidae
#         'Bibio_marci': '#4ECDC4', 'Dilophus_febrilis': '#4ECDC4',  # Bibionidae
#         'Drosophila_melanogaster': '#45B7D1', 'Drosophila_simulans': '#45B7D1',  # Drosophilidae
#         'Bactrocera_dorsalis': '#FFA07A',  # Tephritidae
#         'Calliphora_vicina': '#98D8C8',  # Calliphoridae  
#         'Chironomus_riparius': '#F7DC6F'  # Chironomidae
#     }
    
#     # Draw chromosomes and create layouts
#     layouts = {}
#     for i, species in enumerate(species_names):
#         y_pos = len(species_names) - i - 1
        
#         # Get top chromosomes
#         chroms = common_datasets[species].groupby('sequence').agg({
#             'start': 'min', 'end': 'max', 'busco_id': 'count'
#         }).reset_index()
#         chroms.columns = ['sequence', 'chr_start', 'chr_end', 'gene_count']
#         chroms = chroms.sort_values('gene_count', ascending=False).head(6)
        
#         # Calculate proportional layout
#         total_genes = chroms['gene_count'].sum()
#         chroms['prop_length'] = chroms['gene_count'] / total_genes * 400
#         chroms['cum_start'] = chroms['prop_length'].cumsum() - chroms['prop_length']
#         chroms['cum_end'] = chroms['cum_start'] + chroms['prop_length']
#         chroms['chr_name'] = [f'Chr{j+1}' for j in range(len(chroms))]
        
#         layouts[species] = chroms
        
#         # Draw chromosomes
#         color = family_colors.get(species, '#CCCCCC')
#         for _, chr_data in chroms.iterrows():
#             rect = Rectangle((chr_data['cum_start'], y_pos - 0.12), 
#                            chr_data['prop_length'], 0.24,
#                            facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.8)
#             ax.add_patch(rect)
            
#             # Chromosome labels
#             ax.text(chr_data['cum_start'] + chr_data['prop_length']/2, y_pos + 0.3,
#                    chr_data['chr_name'], ha='center', va='center', 
#                    fontsize=9, fontweight='bold')
        
#         # Species labels
#         ax.text(-40, y_pos, species.replace('_', ' '), ha='right', va='center',
#                fontsize=11, fontstyle='italic', fontweight='bold')
    
#     # Create ENHANCED SYNTENY BLOCKS with confidence scoring
#     print("\nCreating enhanced synteny blocks...")
    
#     reference_species = species_names[0]  # Aedes as reference
#     ref_y = len(species_names) - 1
    
#     # For each target species, create synteny blocks
#     all_connections = []
    
#     for target_idx, target_species in enumerate(species_names[1:], 1):
#         target_y = len(species_names) - target_idx - 1
        
#         # Find pairwise synteny
#         pair_common = common_datasets[reference_species].merge(
#             common_datasets[target_species], on='busco_id', suffixes=('_ref', '_target')
#         )
        
#         # Create synteny blocks (chromosome pairs with >50 genes)
#         synteny_blocks = pair_common.groupby(['sequence_ref', 'sequence_target']).agg({
#             'busco_id': 'count',
#             'start_ref': ['min', 'max'],
#             'start_target': ['min', 'max']
#         }).reset_index()
        
#         # Flatten column names
#         synteny_blocks.columns = ['chr_ref', 'chr_target', 'gene_count', 
#                                 'ref_start', 'ref_end', 'target_start', 'target_end']
        
#         # Filter for significant blocks (>50 genes)
#         synteny_blocks = synteny_blocks[synteny_blocks['gene_count'] >= 50]
        
#         # Calculate confidence (gene_count / max_possible)
#         max_genes = synteny_blocks['gene_count'].max()
#         synteny_blocks['confidence'] = synteny_blocks['gene_count'] / max_genes
        
#         print(f"  {target_species}: {len(synteny_blocks)} synteny blocks")
        
#         # Draw enhanced ribbons
#         for _, block in synteny_blocks.iterrows():
#             # Find chromosome positions
#             ref_chr = layouts[reference_species][
#                 layouts[reference_species]['sequence'] == block['chr_ref']
#             ]
#             target_chr = layouts[target_species][
#                 layouts[target_species]['sequence'] == block['chr_target']
#             ]
            
#             if len(ref_chr) > 0 and len(target_chr) > 0:
#                 # Calculate positions within chromosomes
#                 ref_rel_start = (block['ref_start'] - ref_chr['chr_start'].iloc[0]) / (ref_chr['chr_end'].iloc[0] - ref_chr['chr_start'].iloc[0])
#                 ref_rel_end = (block['ref_end'] - ref_chr['chr_start'].iloc[0]) / (ref_chr['chr_end'].iloc[0] - ref_chr['chr_start'].iloc[0])
                
#                 target_rel_start = (block['target_start'] - target_chr['chr_start'].iloc[0]) / (target_chr['chr_end'].iloc[0] - target_chr['chr_start'].iloc[0])
#                 target_rel_end = (block['target_end'] - target_chr['chr_start'].iloc[0]) / (target_chr['chr_end'].iloc[0] - target_chr['chr_start'].iloc[0])
                
#                 # Calculate plot positions
#                 ref_start_pos = ref_chr['cum_start'].iloc[0] + ref_rel_start * ref_chr['prop_length'].iloc[0]
#                 ref_end_pos = ref_chr['cum_start'].iloc[0] + ref_rel_end * ref_chr['prop_length'].iloc[0]
#                 target_start_pos = target_chr['cum_start'].iloc[0] + target_rel_start * target_chr['prop_length'].iloc[0]
#                 target_end_pos = target_chr['cum_start'].iloc[0] + target_rel_end * target_chr['prop_length'].iloc[0]
                
#                 # Create thick curved ribbon
#                 n_points = 30
#                 t = np.linspace(0, 1, n_points)
                
#                 # Top edge (ref_start to target_start)
#                 top_x = ref_start_pos * (1-t) + target_start_pos * t
#                 top_y = ref_y * (1-t)**2 + (ref_y + target_y)/2 * 2*(1-t)*t + target_y * t**2
                
#                 # Bottom edge (ref_end to target_end)
#                 bottom_x = ref_end_pos * (1-t) + target_end_pos * t
#                 bottom_y = ref_y * (1-t)**2 + (ref_y + target_y)/2 * 2*(1-t)*t + target_y * t**2
                
#                 # Create ribbon polygon
#                 ribbon_x = np.concatenate([top_x, bottom_x[::-1]])
#                 ribbon_y = np.concatenate([top_y, bottom_y[::-1]])
                
#                 # Color and transparency based on confidence
#                 color = plt.cm.viridis(block['confidence'])
#                 alpha = 0.4 + (block['confidence'] * 0.5)
                
#                 # Draw ribbon
#                 ax.fill(ribbon_x, ribbon_y, color=color, alpha=alpha, 
#                        edgecolor='black', linewidth=0.5)
                
#                 all_connections.append({
#                     'species': target_species,
#                     'genes': block['gene_count'],
#                     'confidence': block['confidence']
#                 })
    
#     # Add colorbar for confidence
#     sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
#     sm.set_array([])
#     cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
#     cbar.set_label('Synteny Block Confidence', fontsize=12)
    
#     # Formatting
#     ax.set_xlim(-100, 450)
#     ax.set_ylim(-0.5, len(species_names) - 0.5)
#     ax.axis('off')
#     ax.set_title(f'Enhanced Multi-Species Synteny Analysis\n'
#                 f'Synteny Blocks with Confidence Scoring ({len(common_busco_ids)} common genes)', 
#                 fontsize=16, fontweight='bold', pad=20)
    
#     plt.tight_layout()
#     plt.savefig('working_enhanced_synteny.png', dpi=300, bbox_inches='tight')
#     plt.show()
    
#     # Summary statistics
#     print(f"\n=== ENHANCED SYNTENY SUMMARY ===")
#     print(f"Total connections drawn: {len(all_connections)}")
#     print(f"High confidence blocks (>0.8): {sum(1 for c in all_connections if c['confidence'] > 0.8)}")
#     print(f"Medium confidence blocks (0.5-0.8): {sum(1 for c in all_connections if 0.5 <= c['confidence'] <= 0.8)}")
#     print(f"Average genes per block: {np.mean([c['genes'] for c in all_connections]):.1f}")
    
#     print(f"\n✓ Enhanced synteny plot saved: working_enhanced_synteny.png")

# # Run the working enhanced version
# create_working_enhanced_synteny()



# #!/usr/bin/env python3
# """
# BUSCO File Diagnostics and Quick Fix
# Analyze and fix BUSCO parsing issues
# """

# import pandas as pd
# import sys

# def diagnose_busco_file(busco_path):
#     """Diagnose BUSCO file structure and content"""
#     print(f"\n=== DIAGNOSING BUSCO FILE: {busco_path} ===")
    
#     try:
#         # Read first few lines to understand structure
#         with open(busco_path, 'r') as f:
#             lines = f.readlines()
        
#         # Find non-comment lines
#         non_comment_lines = [line for line in lines if not line.startswith('#')]
        
#         print(f"Total lines: {len(lines)}")
#         print(f"Comment lines: {len(lines) - len(non_comment_lines)}")
#         print(f"Data lines: {len(non_comment_lines)}")
        
#         if len(non_comment_lines) > 0:
#             # Show first few data lines
#             print(f"\nFirst 3 data lines:")
#             for i, line in enumerate(non_comment_lines[:3]):
#                 parts = line.strip().split('\t')
#                 print(f"Line {i+1}: {len(parts)} columns")
#                 print(f"  Raw: {repr(line.strip())}")
#                 print(f"  Parts: {parts}")
        
#         # Analyze column structure
#         valid_entries = []
#         invalid_entries = []
        
#         for line_num, line in enumerate(non_comment_lines, 1):
#             if line.strip():
#                 parts = line.strip().split('\t')
#                 if len(parts) >= 6:
#                     # Check if coordinates are valid
#                     try:
#                         start = parts[3]
#                         end = parts[4]
                        
#                         if start != 'N/A' and end != 'N/A':
#                             start_int = int(start)
#                             end_int = int(end)
#                             if start_int < end_int:
#                                 valid_entries.append(line_num)
#                             else:
#                                 invalid_entries.append((line_num, "start >= end"))
#                         else:
#                             invalid_entries.append((line_num, "N/A coordinates"))
#                     except ValueError as e:
#                         invalid_entries.append((line_num, f"Invalid number: {e}"))
#                 else:
#                     invalid_entries.append((line_num, f"Too few columns: {len(parts)}"))
        
#         print(f"\nCOORDINATE ANALYSIS:")
#         print(f"Valid coordinate entries: {len(valid_entries)}")
#         print(f"Invalid coordinate entries: {len(invalid_entries)}")
        
#         if len(invalid_entries) > 0:
#             print(f"\nFirst 5 invalid entries:")
#             for line_num, reason in invalid_entries[:5]:
#                 print(f"  Line {line_num}: {reason}")
        
#         # Check status distribution
#         status_counts = {}
#         for line in non_comment_lines:
#             if line.strip():
#                 parts = line.strip().split('\t')
#                 if len(parts) >= 2:
#                     status = parts[1]
#                     status_counts[status] = status_counts.get(status, 0) + 1
        
#         print(f"\nSTATUS DISTRIBUTION:")
#         for status, count in status_counts.items():
#             print(f"  {status}: {count}")
        
#         return len(valid_entries), len(invalid_entries)
        
#     except Exception as e:
#         print(f"ERROR reading file: {e}")
#         return 0, 0

# def create_fixed_busco_parser():
#     """Create a fixed BUSCO parser that handles real BUSCO output"""
    
#     code = '''
# def fixed_parse_busco_table(busco_path, config):
#     """Fixed BUSCO table parsing that handles real BUSCO output properly"""
#     logger.info(f"Parsing BUSCO table: {busco_path}")
    
#     with open(busco_path, 'r') as f:
#         lines = [line for line in f if not line.startswith('#')]
    
#     busco_data = []
#     parsing_stats = {
#         'total_lines': len(lines),
#         'valid_entries': 0,
#         'invalid_coordinates': 0,
#         'missing_columns': 0,
#         'na_coordinates': 0,
#         'other_errors': 0
#     }
    
#     for line_num, line in enumerate(lines, 1):
#         if line.strip():
#             parts = line.strip().split('\\t')
            
#             # Check minimum columns
#             if len(parts) < 6:
#                 parsing_stats['missing_columns'] += 1
#                 continue
            
#             try:
#                 # Extract basic info
#                 busco_id = parts[0]
#                 status = parts[1]
#                 sequence = parts[2]
#                 start_str = parts[3]
#                 end_str = parts[4]
#                 strand = parts[5] if len(parts) > 5 else '+'
                
#                 # Only process Complete BUSCOs with valid coordinates
#                 if status == 'Complete' and start_str != 'N/A' and end_str != 'N/A':
#                     try:
#                         gene_start = int(start_str)
#                         gene_end = int(end_str)
                        
#                         if gene_start >= gene_end:
#                             parsing_stats['invalid_coordinates'] += 1
#                             continue
                        
#                         # Extract additional info if available
#                         score = None
#                         length = None
#                         if len(parts) > 6 and parts[6] != 'N/A':
#                             try:
#                                 score = float(parts[6])
#                             except:
#                                 pass
#                         if len(parts) > 7 and parts[7] != 'N/A':
#                             try:
#                                 length = int(parts[7])
#                             except:
#                                 length = gene_end - gene_start
#                         else:
#                             length = gene_end - gene_start
                        
#                         entry = {
#                             'busco_id': busco_id,
#                             'status': status,
#                             'sequence': sequence,
#                             'gene_start': gene_start,
#                             'gene_end': gene_end,
#                             'strand': strand,
#                             'score': score,
#                             'length': length,
#                             'line_number': line_num
#                         }
                        
#                         busco_data.append(entry)
#                         parsing_stats['valid_entries'] += 1
                        
#                     except ValueError:
#                         parsing_stats['invalid_coordinates'] += 1
                        
#                 elif start_str == 'N/A' or end_str == 'N/A':
#                     parsing_stats['na_coordinates'] += 1
#                 else:
#                     parsing_stats['other_errors'] += 1
                    
#             except Exception as e:
#                 parsing_stats['other_errors'] += 1
    
#     busco_df = pd.DataFrame(busco_data)
    
#     # Add paralog detection if needed
#     if config.get('enable_paralog_detection', False) and len(busco_df) > 0:
#         busco_df = detect_and_annotate_paralogs(busco_df, config)
    
#     # Report parsing statistics
#     logger.info(f"  BUSCO Parsing Results:")
#     logger.info(f"    Total lines processed: {parsing_stats['total_lines']}")
#     logger.info(f"    Valid Complete BUSCOs: {parsing_stats['valid_entries']}")
#     logger.info(f"    N/A coordinates: {parsing_stats['na_coordinates']}")
#     logger.info(f"    Invalid coordinates: {parsing_stats['invalid_coordinates']}")
#     logger.info(f"    Missing columns: {parsing_stats['missing_columns']}")
#     logger.info(f"    Other errors: {parsing_stats['other_errors']}")
    
#     if parsing_stats['valid_entries'] == 0:
#         logger.error(f"  NO VALID BUSCO ENTRIES FOUND!")
#         logger.error(f"  This will cause downstream failures.")
    
#     logger.info(f"  Found {len(busco_df)} valid BUSCO entries")
#     return busco_df
#     '''
    
#     return code

# def main():
#     if len(sys.argv) != 2:
#         print("Usage: python busco_diagnostics.py <path_to_full_table.tsv>")
#         sys.exit(1)
    
#     busco_file = sys.argv[1]
    
#     # Diagnose the BUSCO file
#     valid_count, invalid_count = diagnose_busco_file(busco_file)
    
#     print(f"\n=== SUMMARY ===")
#     print(f"Valid entries: {valid_count}")
#     print(f"Invalid entries: {invalid_count}")
    
#     if valid_count == 0:
#         print(f"\n❌ CRITICAL: No valid BUSCO entries found!")
#         print(f"This explains why your ortholog mapping fails.")
#         print(f"\nPossible causes:")
#         print(f"1. BUSCO output format is different than expected")
#         print(f"2. All genes have N/A coordinates (incomplete annotation)")
#         print(f"3. BUSCO run failed or was incomplete")
#         print(f"\nRecommendations:")
#         print(f"1. Re-run BUSCO with proper genome annotation")
#         print(f"2. Check if you have the right BUSCO output file")
#         print(f"3. Use the fixed parser in the generated code")
#     else:
#         print(f"\n✅ Found {valid_count} valid entries")
#         print(f"The fixed parser should work with your data")
    
#     print(f"\n=== FIXED PARSER CODE ===")
#     print("Replace the enhanced_parse_busco_table function with:")
#     print(create_fixed_busco_parser())

# if __name__ == "__main__":
#     main()



#!/usr/bin/env python3
"""
Simple BUSCO file checker
"""

def check_busco_file(filename):
    print(f"\n=== Checking {filename} ===")
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        print(f"Total lines: {len(lines)}")
        
        # Show first 5 lines
        print("\nFirst 5 lines:")
        for i, line in enumerate(lines[:5]):
            parts = line.strip().split('\t')
            print(f"Line {i+1}: {len(parts)} columns")
            print(f"  Content: {repr(line.strip())}")
            if len(parts) >= 6:
                print(f"  Parsed: ID={parts[0]}, Status={parts[1]}, Seq={parts[2]}, Start={parts[3]}, End={parts[4]}, Strand={parts[5]}")
            print()
        
        # Count by status
        status_counts = {}
        complete_with_coords = 0
        
        for line in lines:
            if not line.startswith('#') and line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    status = parts[1]
                    status_counts[status] = status_counts.get(status, 0) + 1
                    
                    if status == 'Complete' and len(parts) >= 6:
                        if parts[3] != 'N/A' and parts[4] != 'N/A':
                            try:
                                start = int(parts[3])
                                end = int(parts[4])
                                if start < end:
                                    complete_with_coords += 1
                            except:
                                pass
        
        print("Status counts:")
        for status, count in status_counts.items():
            print(f"  {status}: {count}")
        
        print(f"Complete BUSCOs with valid coordinates: {complete_with_coords}")
        
        if complete_with_coords == 0:
            print("❌ NO COMPLETE BUSCOS WITH VALID COORDINATES!")
            print("This explains the parsing failures.")
        else:
            print(f"✅ Found {complete_with_coords} usable BUSCOs")
        
    except Exception as e:
        print(f"Error reading file: {e}")

# Check both BUSCO files
if __name__ == "__main__":
    check_busco_file('Bibio_marci/full_table.tsv')
    check_busco_file('Dilophus_febrilis/full_table.tsv')