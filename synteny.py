# #this is a test workflow

# # Enhanced Multi-Species Synteny with Thick Synteny Blocks
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.patches import Rectangle, Polygon
# import os
# import glob

# def create_synteny_blocks(ref_genes, target_genes, min_block_size=10):
#     """Create synteny blocks instead of individual gene lines"""
#     print(f"Creating synteny blocks (min {min_block_size} genes per block)...")
    
#     # Find common genes
#     common = ref_genes.merge(target_genes, on='busco_id', suffixes=('_ref', '_target'))
    
#     # Group by chromosome pairs
#     synteny_blocks = []
    
#     for (ref_chr, target_chr), group in common.groupby(['sequence_ref', 'sequence_target']):
#         if len(group) >= min_block_size:
#             # Sort by position in reference
#             group_sorted = group.sort_values('start_ref')
            
#             # Create blocks of consecutive genes
#             block_start = 0
#             current_block = []
            
#             for i, (_, gene) in enumerate(group_sorted.iterrows()):
#                 current_block.append(gene)
                
#                 # End block if we have enough genes or at end
#                 if len(current_block) >= min_block_size or i == len(group_sorted) - 1:
#                     if len(current_block) >= min_block_size:
#                         block_data = pd.DataFrame(current_block)
#                         synteny_blocks.append({
#                             'ref_chr': ref_chr,
#                             'target_chr': target_chr,
#                             'ref_start': block_data['start_ref'].min(),
#                             'ref_end': block_data['end_ref'].max(),
#                             'target_start': block_data['start_target'].min(),
#                             'target_end': block_data['end_target'].max(),
#                             'gene_count': len(current_block),
#                             'genes': list(block_data['busco_id'])
#                         })
#                     current_block = []
    
#     blocks_df = pd.DataFrame(synteny_blocks)
#     print(f"  Created {len(blocks_df)} synteny blocks")
#     return blocks_df

# def draw_thick_synteny_block(ax, ref_pos_start, ref_pos_end, target_pos_start, target_pos_end, 
#                             ref_y, target_y, color, gene_count):
#     """Draw thick synteny block as polygon"""
    
#     # Calculate thickness based on gene count
#     thickness = min(0.05 + (gene_count / 200), 0.3)  # Thicker for more genes
    
#     # Create polygon points for thick ribbon
#     points = np.array([
#         [ref_pos_start, ref_y - thickness/2],      # Top left
#         [ref_pos_end, ref_y - thickness/2],        # Top right  
#         [target_pos_end, target_y + thickness/2],  # Bottom right
#         [target_pos_start, target_y + thickness/2], # Bottom left
#         [ref_pos_start, ref_y - thickness/2]       # Close polygon
#     ])
    
#     # Create smooth curved ribbon using bezier-like interpolation
#     n_points = 50
#     t = np.linspace(0, 1, n_points)
    
#     # Top edge curve
#     top_x = ref_pos_start * (1-t) + target_pos_start * t
#     top_y = (ref_y - thickness/2) * (1-t)**2 + (ref_y + target_y)/2 * 2*(1-t)*t + (target_y + thickness/2) * t**2
    
#     # Bottom edge curve  
#     bottom_x = ref_pos_end * (1-t) + target_pos_end * t
#     bottom_y = (ref_y + thickness/2) * (1-t)**2 + (ref_y + target_y)/2 * 2*(1-t)*t + (target_y - thickness/2) * t**2
    
#     # Combine into polygon
#     ribbon_x = np.concatenate([top_x, bottom_x[::-1]])
#     ribbon_y = np.concatenate([top_y, bottom_y[::-1]])
    
#     # Draw filled polygon
#     polygon = Polygon(list(zip(ribbon_x, ribbon_y)), 
#                      facecolor=color, edgecolor=color, 
#                      alpha=0.7, linewidth=0.5, zorder=10)  # High zorder = in front
#     ax.add_patch(polygon)

# def create_enhanced_multi_species_synteny():
#     """Create enhanced visualization with thick synteny blocks"""
#     print("=== ENHANCED MULTI-SPECIES SYNTENY WITH THICK BLOCKS ===\n")
    
#     # [Previous data loading code here - same as before]
#     # ... load_all_species_data(), find_all_common_genes(), etc.
    
#     # For demonstration, let's focus on key species
#     key_species = ['Bibio_marci', 'Dilophus_febrilis', 'Drosophila_melanogaster', 
#                    'Aedes_aegypti', 'Calliphora_vicina']
    
#     print("Creating enhanced plot for key species...")
    
#     fig, ax = plt.subplots(figsize=(18, 10))
    
#     # Draw chromosomes (same as before but with lower zorder)
#     for i, species in enumerate(key_species):
#         y_pos = len(key_species) - i - 1
        
#         # Draw chromosome rectangles (behind synteny blocks)
#         for j in range(6):  # Assume 6 chromosomes
#             rect = Rectangle((j*60, y_pos-0.15), 50, 0.3,
#                            facecolor='lightgray', edgecolor='black', 
#                            linewidth=2, zorder=1)  # Low zorder = behind
#             ax.add_patch(rect)
            
#             # Chromosome labels
#             ax.text(j*60 + 25, y_pos + 0.35, f'Chr{j+1}',
#                    ha='center', va='center', fontsize=10, fontweight='bold')
        
#         # Species labels
#         ax.text(-20, y_pos, species.replace('_', ' '), 
#                ha='right', va='center', fontsize=12, fontstyle='italic')
    
#     # Create thick synteny blocks between species
#     reference_y = len(key_species) - 1  # Bibio at top
#     colors = plt.cm.Set3(np.linspace(0, 1, 12))
    
#     for i, target_species in enumerate(key_species[1:], 1):
#         target_y = len(key_species) - i - 1
        
#         # Simulate synteny blocks (replace with real data)
#         for block_id in range(8):  # 8 major synteny blocks
#             ref_start = np.random.uniform(block_id*40, block_id*40 + 30)
#             ref_end = ref_start + np.random.uniform(15, 25)
#             target_start = np.random.uniform(block_id*40, block_id*40 + 30)  
#             target_end = target_start + np.random.uniform(15, 25)
#             gene_count = np.random.randint(50, 200)
            
#             draw_thick_synteny_block(ax, ref_start, ref_end, target_start, target_end,
#                                    reference_y, target_y, colors[block_id], gene_count)
    
#     ax.set_xlim(-50, 400)
#     ax.set_ylim(-0.5, len(key_species) - 0.5)
#     ax.axis('off')
#     ax.set_title('Enhanced Synteny Analysis - Thick Synteny Blocks\n(High Gene Density Regions)', 
#                 fontsize=16, fontweight='bold')
    
#     plt.tight_layout()
#     plt.savefig('enhanced_thick_synteny_blocks.png', dpi=300, bbox_inches='tight')
#     plt.show()

# # Options for different visualizations:
# def visualization_options():
#     """Different ways to represent the synteny data"""
    
#     print("VISUALIZATION OPTIONS:")
#     print("="*50)
#     print("1. CURRENT: Individual gene lines (what we have)")
#     print("   - Each line = 1 BUSCO gene")
#     print("   - 4,057 thin lines")
#     print("   - Shows precise gene-to-gene relationships")
    
#     print("\n2. SYNTENY BLOCKS: Grouped gene regions")
#     print("   - Each block = 10-50 consecutive genes")
#     print("   - Thick colored ribbons")
#     print("   - Shows major conserved regions")
    
#     print("\n3. CHROMOSOME ARMS: Major genomic segments")
#     print("   - Each ribbon = entire chromosome arm")
#     print("   - Very thick connections")
#     print("   - Shows large-scale rearrangements")
    
#     print("\n4. DENSITY HEATMAP: Gene density visualization")
#     print("   - Color intensity = number of genes")
#     print("   - No individual lines")
#     print("   - Shows conservation hotspots")
    
#     print("\n5. CURVED RIBBONS: Bezier curve connections")
#     print("   - Smooth curved synteny blocks")
#     print("   - Variable thickness by gene count")
#     print("   - More aesthetically pleasing")

# if __name__ == "__main__":
#     visualization_options()
#     print("\nWhich visualization would you like me to implement?")


