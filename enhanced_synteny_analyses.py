# Enhanced Multi-Species Synteny Analysis with Synteny Blocks
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Polygon
import seaborn as sns
import networkx as nx
from collections import defaultdict, Counter
import os
import glob
import itertools

class SyntenyBlock:
    """Class to represent a synteny block"""
    def __init__(self, genes, chr1, chr2, species1, species2):
        self.genes = genes
        self.chr1 = chr1
        self.chr2 = chr2
        self.species1 = species1
        self.species2 = species2
        self.size = len(genes)
        self.block_id = f"{chr1}_{chr2}_{self.size}"
        
    def get_positions(self):
        return {
            'start1': self.genes['start_1'].min(),
            'end1': self.genes['end_1'].max(),
            'start2': self.genes['start_2'].min(), 
            'end2': self.genes['end_2'].max()
        }

def load_busco_data(file_path, species_name):
    """Load BUSCO data with error handling"""
    print(f"Loading {species_name} data...")
    
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Find header line
        header_idx = None
        for i, line in enumerate(lines):
            if line.startswith('# Busco id'):
                header_idx = i
                break
        
        if header_idx is None:
            print(f"ERROR: No header found in {species_name}!")
            return pd.DataFrame()
        
        # Parse data lines
        data_lines = []
        for line in lines[header_idx + 1:]:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                data_lines.append(parts)
        
        # Create dataframe
        df = pd.DataFrame(data_lines, columns=[
            'busco_id', 'status', 'sequence', 'start', 'end', 
            'strand', 'score', 'length'
        ])
        
        # Filter and convert
        df = df[df['status'] == 'Complete'].copy()
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        df = df.dropna(subset=['start', 'end'])
        df['species'] = species_name
        
        print(f"  Loaded {len(df)} complete genes")
        return df
        
    except Exception as e:
        print(f"ERROR loading {species_name}: {e}")
        return pd.DataFrame()

def create_synteny_blocks(species1_data, species2_data, min_block_size=10):
    """Create synteny blocks between two species"""
    print(f"Creating synteny blocks between {species1_data['species'].iloc[0]} and {species2_data['species'].iloc[0]}...")
    
    # Find common genes
    common_genes = species1_data.merge(species2_data, on='busco_id', suffixes=('_1', '_2'))
    
    # Group by chromosome pairs
    synteny_blocks = []
    for (chr1, chr2), group in common_genes.groupby(['sequence_1', 'sequence_2']):
        if len(group) >= min_block_size:
            # Sort by position in species 1
            group_sorted = group.sort_values('start_1')
            
            # Create blocks of consecutive genes (within 1MB windows)
            current_block = []
            last_pos = -1e9
            
            for _, gene in group_sorted.iterrows():
                if gene['start_1'] - last_pos > 1e6:  # 1MB gap threshold
                    if len(current_block) >= min_block_size:
                        block_genes = pd.DataFrame(current_block)
                        block = SyntenyBlock(block_genes, chr1, chr2, 
                                          species1_data['species'].iloc[0],
                                          species2_data['species'].iloc[0])
                        synteny_blocks.append(block)
                    current_block = []
                
                current_block.append(gene)
                last_pos = gene['start_1']
            
            # Add final block
            if len(current_block) >= min_block_size:
                block_genes = pd.DataFrame(current_block)
                block = SyntenyBlock(block_genes, chr1, chr2,
                                  species1_data['species'].iloc[0], 
                                  species2_data['species'].iloc[0])
                synteny_blocks.append(block)
    
    print(f"  Created {len(synteny_blocks)} synteny blocks")
    return synteny_blocks

def calculate_block_confidence(all_blocks, reference_species='Bibio_marci'):
    """Calculate confidence scores for synteny blocks"""
    print("Calculating block confidence scores...")
    
    # Group blocks by chromosome pairs
    block_support = defaultdict(list)
    for block in all_blocks:
        key = (block.chr1, block.chr2)
        block_support[key].append(block)
    
    # Calculate confidence based on species support
    confidence_scores = {}
    total_species = len(set([block.species1 for block in all_blocks] + 
                           [block.species2 for block in all_blocks]))
    
    for key, blocks in block_support.items():
        species_supporting = set([block.species1 for block in blocks] + 
                                [block.species2 for block in blocks])
        confidence = len(species_supporting) / total_species
        
        for block in blocks:
            confidence_scores[block.block_id] = {
                'confidence': confidence,
                'species_support': len(species_supporting),
                'total_genes': sum([b.size for b in blocks])
            }
    
    return confidence_scores

def create_synteny_network(all_blocks, confidence_scores, min_confidence=0.3):
    """Create network representation of synteny blocks"""
    print("Creating synteny network...")
    
    G = nx.Graph()
    
    # Add nodes for each synteny block
    for block in all_blocks:
        if confidence_scores[block.block_id]['confidence'] >= min_confidence:
            G.add_node(block.block_id, 
                      size=block.size,
                      confidence=confidence_scores[block.block_id]['confidence'],
                      species1=block.species1,
                      species2=block.species2,
                      chr1=block.chr1,
                      chr2=block.chr2)
    
    # Add edges between blocks that share chromosomes
    for block1, block2 in itertools.combinations(all_blocks, 2):
        if (confidence_scores[block1.block_id]['confidence'] >= min_confidence and
            confidence_scores[block2.block_id]['confidence'] >= min_confidence):
            
            # Connect if they share a chromosome in either species
            if (block1.chr1 == block2.chr1 or block1.chr2 == block2.chr2 or
                block1.chr1 == block2.chr2 or block1.chr2 == block2.chr1):
                
                edge_weight = (confidence_scores[block1.block_id]['confidence'] + 
                             confidence_scores[block2.block_id]['confidence']) / 2
                G.add_edge(block1.block_id, block2.block_id, weight=edge_weight)
    
    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G

def plot_synteny_network(G, output_file='synteny_network.png'):
    """Plot the synteny network"""
    print("Plotting synteny network...")
    
    if G.number_of_nodes() == 0:
        print("No nodes to plot!")
        return
    
    plt.figure(figsize=(16, 12))
    
    # Layout
    pos = nx.spring_layout(G, k=3, iterations=50)
    
    # Node properties
    node_sizes = [G.nodes[node]['size'] * 20 for node in G.nodes()]
    node_colors = [G.nodes[node]['confidence'] for node in G.nodes()]
    
    # Edge properties  
    edge_weights = [G[u][v]['weight'] * 5 for u, v in G.edges()]
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                          cmap='viridis', alpha=0.8)
    nx.draw_networkx_edges(G, pos, width=edge_weights, alpha=0.6, edge_color='gray')
    
    # Labels for largest nodes only
    large_nodes = {node: node.split('_')[-1] for node in G.nodes() 
                  if G.nodes[node]['size'] > 20}
    nx.draw_networkx_labels(G, pos, large_nodes, font_size=8)
    
    plt.colorbar(plt.cm.ScalarMappable(cmap='viridis'), 
                label='Block Confidence', shrink=0.8)
    plt.title('Synteny Block Network\n(Node size = genes, Color = confidence)', 
             fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

def create_enhanced_ribbon_plot(all_species_data, common_datasets, synteny_blocks, 
                               confidence_scores, output_file='enhanced_synteny_ribbons.png'):
    """Create enhanced ribbon plot with synteny blocks"""
    print("Creating enhanced ribbon plot...")
    
    # Select top species for visualization
    key_species = list(common_datasets.keys())[:8]  # Top 8 species
    
    fig, ax = plt.subplots(figsize=(20, len(key_species) * 1.5))
    
    # Family colors
    family_colors = {
        'Bibio_marci': '#FF6B6B', 'Dilophus_febrilis': '#FF6B6B', 'Plecia_longiforceps': '#FF6B6B',
        'Tipula_lateralis': '#4ECDC4', 'Nephrotoma_appendiculata': '#4ECDC4',
        'Chironomus_riparius': '#45B7D1',
        'Drosophila_melanogaster': '#FFA07A', 'Drosophila_simulans': '#FFA07A',
        'Aedes_aegypti': '#98D8C8', 'Anopheles_gambia': '#98D8C8', 'Culex_pipiens': '#98D8C8'
    }
    
    # Create chromosome layouts
    layouts = {}
    for i, species in enumerate(key_species):
        species_data = common_datasets[species]
        
        # Get top chromosomes
        chroms = species_data.groupby('sequence').agg({
            'start': 'min', 'end': 'max', 'busco_id': 'count'
        }).reset_index()
        chroms.columns = ['sequence', 'chr_start', 'chr_end', 'gene_count']
        chroms = chroms.sort_values('gene_count', ascending=False).head(6)
        
        # Layout
        total_genes = chroms['gene_count'].sum()
        chroms['prop_length'] = chroms['gene_count'] / total_genes * 300
        chroms['cum_start'] = chroms['prop_length'].cumsum() - chroms['prop_length']
        chroms['cum_end'] = chroms['cum_start'] + chroms['prop_length']
        chroms['chr_name'] = [f'Chr{j+1}' for j in range(len(chroms))]
        
        layouts[species] = chroms
        
        # Draw chromosomes
        y_pos = len(key_species) - i - 1
        color = family_colors.get(species, '#CCCCCC')
        
        for _, chr_data in chroms.iterrows():
            rect = Rectangle((chr_data['cum_start'], y_pos - 0.15), 
                           chr_data['prop_length'], 0.3,
                           facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.7)
            ax.add_patch(rect)
            
            # Labels
            ax.text(chr_data['cum_start'] + chr_data['prop_length']/2, y_pos + 0.35,
                   chr_data['chr_name'], ha='center', va='center', 
                   fontsize=9, fontweight='bold')
        
        # Species label
        ax.text(-30, y_pos, species.replace('_', ' '), ha='right', va='center',
               fontsize=11, fontstyle='italic', fontweight='bold')
    
    # Draw synteny block connections
    reference_species = key_species[0]
    ref_y = len(key_species) - 1
    
    high_confidence_blocks = [block for block in synteny_blocks 
                             if confidence_scores[block.block_id]['confidence'] > 0.5]
    
    print(f"Drawing {len(high_confidence_blocks)} high-confidence synteny blocks...")
    
    for block in high_confidence_blocks[:100]:  # Limit for visibility
        if (block.species1 in key_species and block.species2 in key_species and
            block.species1 != block.species2):
            
            # Get positions
            species1_idx = key_species.index(block.species1)
            species2_idx = key_species.index(block.species2)
            y1 = len(key_species) - species1_idx - 1
            y2 = len(key_species) - species2_idx - 1
            
            # Get chromosome layouts
            layout1 = layouts[block.species1]
            layout2 = layouts[block.species2]
            
            # Find chromosome positions
            chr1_data = layout1[layout1['sequence'] == block.chr1]
            chr2_data = layout2[layout2['sequence'] == block.chr2]
            
            if len(chr1_data) > 0 and len(chr2_data) > 0:
                pos1 = chr1_data['cum_start'].iloc[0] + chr1_data['prop_length'].iloc[0] / 2
                pos2 = chr2_data['cum_start'].iloc[0] + chr2_data['prop_length'].iloc[0] / 2
                
                # Draw thick ribbon based on block size
                thickness = min(0.05 + (block.size / 100), 0.2)
                confidence = confidence_scores[block.block_id]['confidence']
                alpha = 0.3 + (confidence * 0.5)
                
                # Create curved ribbon
                n_points = 20
                t = np.linspace(0, 1, n_points)
                x_curve = pos1 * (1-t) + pos2 * t
                y_curve = y1 * (1-t)**2 + (y1+y2)/2 * 2*(1-t)*t + y2 * t**2
                
                # Draw ribbon as thick line
                ax.plot(x_curve, y_curve, color='purple', linewidth=thickness*20, 
                       alpha=alpha, solid_capstyle='round')
    
    # Formatting
    ax.set_xlim(-80, 350)
    ax.set_ylim(-0.5, len(key_species) - 0.5)
    ax.axis('off')
    ax.set_title(f'Enhanced Multi-Species Synteny Analysis\n'
                f'Synteny Blocks (thickness ∝ block size, opacity ∝ confidence)',
                fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

def analyze_synteny_evolution(synteny_blocks, confidence_scores):
    """Analyze evolutionary patterns in synteny blocks"""
    print("\n=== SYNTENY EVOLUTION ANALYSIS ===")
    
    # Block size distribution
    block_sizes = [block.size for block in synteny_blocks]
    print(f"Block size statistics:")
    print(f"  Mean: {np.mean(block_sizes):.1f} genes")
    print(f"  Median: {np.median(block_sizes):.1f} genes") 
    print(f"  Max: {np.max(block_sizes)} genes")
    print(f"  Total blocks: {len(synteny_blocks)}")
    
    # Confidence distribution
    confidences = [scores['confidence'] for scores in confidence_scores.values()]
    print(f"\nConfidence statistics:")
    print(f"  Mean confidence: {np.mean(confidences):.3f}")
    print(f"  High confidence blocks (>0.7): {sum(1 for c in confidences if c > 0.7)}")
    print(f"  Medium confidence blocks (0.4-0.7): {sum(1 for c in confidences if 0.4 <= c <= 0.7)}")
    print(f"  Low confidence blocks (<0.4): {sum(1 for c in confidences if c < 0.4)}")
    
    # Family-specific patterns
    species_pairs = defaultdict(list)
    for block in synteny_blocks:
        pair = tuple(sorted([block.species1, block.species2]))
        species_pairs[pair].append(block.size)
    
    print(f"\nTop species pairs by synteny conservation:")
    sorted_pairs = sorted(species_pairs.items(), 
                         key=lambda x: sum(x[1]), reverse=True)
    for pair, sizes in sorted_pairs[:10]:
        print(f"  {pair[0]} - {pair[1]}: {len(sizes)} blocks, {sum(sizes)} total genes")

def main():
    """Enhanced multi-species synteny analysis"""
    print("=== ENHANCED MULTI-SPECIES SYNTENY ANALYSIS ===")
    print("Using synteny blocks and confidence scoring\n")
    
    # Load data
    print("1. Loading all species data...")
    all_species_data = {}
    busco_files = glob.glob('busco-data/*.tsv')
    
    for file_path in sorted(busco_files)[:10]:  # Limit to 10 species for demo
        species_name = os.path.basename(file_path).replace('.tsv', '')
        species_df = load_busco_data(file_path, species_name)
        if len(species_df) > 0:
            all_species_data[species_name] = species_df
    
    if len(all_species_data) < 2:
        print("ERROR: Need at least 2 species!")
        return
    
    # Find common genes
    print("\n2. Finding common genes...")
    species_names = list(all_species_data.keys())
    common_busco_ids = set(all_species_data[species_names[0]]['busco_id'])
    
    for species_name, species_df in all_species_data.items():
        common_busco_ids = common_busco_ids.intersection(set(species_df['busco_id']))
    
    print(f"Found {len(common_busco_ids)} genes common to all species")
    
    # Create common datasets
    common_datasets = {}
    for species_name, species_df in all_species_data.items():
        common_datasets[species_name] = species_df[
            species_df['busco_id'].isin(common_busco_ids)
        ].copy()
    
    # Create synteny blocks
    print("\n3. Creating synteny blocks...")
    all_blocks = []
    
    for i, species1 in enumerate(species_names):
        for species2 in species_names[i+1:]:
            blocks = create_synteny_blocks(
                common_datasets[species1], 
                common_datasets[species2],
                min_block_size=15
            )
            all_blocks.extend(blocks)
    
    # Calculate confidence scores
    print("\n4. Calculating confidence scores...")
    confidence_scores = calculate_block_confidence(all_blocks)
    
    # Create network
    print("\n5. Creating synteny network...")
    G = create_synteny_network(all_blocks, confidence_scores, min_confidence=0.4)
    plot_synteny_network(G)
    
    # Create enhanced ribbon plot
    print("\n6. Creating enhanced ribbon plot...")
    create_enhanced_ribbon_plot(all_species_data, common_datasets, 
                               all_blocks, confidence_scores)
    
    # Analyze evolution
    print("\n7. Analyzing synteny evolution...")
    analyze_synteny_evolution(all_blocks, confidence_scores)
    
    print(f"\n🎉 ENHANCED ANALYSIS COMPLETE!")
    print(f"📊 Created network and ribbon visualizations")
    print(f"🧬 Analyzed {len(all_blocks)} synteny blocks")
    print(f"📈 Used confidence scoring and evolutionary insights")

if __name__ == "__main__":
    main()
    