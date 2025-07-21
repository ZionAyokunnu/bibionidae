#!/usr/bin/env python3
"""
Farm-optimized Syngraph analysis with full visualization capabilities
"""

import pickle
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from collections import defaultdict
import seaborn as sns
import os
import sys

# Set matplotlib backend for cluster use
plt.switch_backend('Agg')

def load_and_analyze_syngraph(pickle_file):
    """Load and comprehensively analyze syngraph data"""
    
    print("=== Farm Syngraph Analysis ===")
    print(f"Loading data from: {pickle_file}")
    
    with open(pickle_file, 'rb') as f:
        data = pickle.load(f)
    
    print(f"Data type: {type(data)}")
    
    # Find the graph
    graph = None
    for attr in ['graph', 'G', 'network', 'syn_graph']:
        if hasattr(data, attr):
            potential_graph = getattr(data, attr)
            if hasattr(potential_graph, 'nodes') and hasattr(potential_graph, 'edges'):
                graph = potential_graph
                print(f"Found graph at attribute: {attr}")
                break
    
    if graph is None:
        print("Exploring data structure...")
        for key, value in data.__dict__.items():
            print(f"  {key}: {type(value)}")
            if hasattr(value, 'nodes') and hasattr(value, 'edges'):
                graph = value
                print(f"  -> Using this as graph!")
                break
    
    if graph is None:
        print("ERROR: Could not find graph in data structure")
        return None
    
    return graph, data

def comprehensive_graph_analysis(graph):
    """Perform comprehensive graph analysis"""
    
    print("\n=== Comprehensive Graph Analysis ===")
    
    # Basic metrics
    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()
    
    print(f"Nodes (genes): {n_nodes:,}")
    print(f"Edges (adjacencies): {n_edges:,}")
    print(f"Average degree: {(2 * n_edges) / n_nodes:.2f}")
    print(f"Graph density: {nx.density(graph):.6f}")
    
    # Degree analysis
    degrees = dict(graph.degree())
    degree_values = list(degrees.values())
    
    print(f"Degree statistics:")
    print(f"  Min degree: {min(degree_values)}")
    print(f"  Max degree: {max(degree_values)}")
    print(f"  Mean degree: {np.mean(degree_values):.2f}")
    print(f"  Median degree: {np.median(degree_values):.2f}")
    
    # Connected components
    components = list(nx.connected_components(graph))
    print(f"Connected components: {len(components)}")
    
    if len(components) > 1:
        sizes = [len(c) for c in components]
        print(f"Component sizes: {sorted(sizes, reverse=True)[:10]}")
    
    # Clustering coefficient
    clustering = nx.average_clustering(graph)
    print(f"Average clustering coefficient: {clustering:.4f}")
    
    return {
        'nodes': n_nodes,
        'edges': n_edges,
        'degrees': degrees,
        'components': components,
        'clustering': clustering
    }

def create_advanced_visualizations(graph, metrics, output_dir="syngraph_plots"):
    """Create advanced visualizations using farm resources"""
    
    print(f"\n=== Creating Advanced Visualizations ===")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Degree distribution plot
    plt.figure(figsize=(15, 10))
    
    degrees = list(metrics['degrees'].values())
    
    plt.subplot(2, 3, 1)
    plt.hist(degrees, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('Degree')
    plt.ylabel('Frequency')
    plt.title('Degree Distribution')
    plt.yscale('log')
    
    plt.subplot(2, 3, 2)
    plt.scatter(range(len(degrees)), sorted(degrees, reverse=True), alpha=0.6, s=10)
    plt.xlabel('Node Rank')
    plt.ylabel('Degree')
    plt.title('Degree Rank Plot')
    plt.yscale('log')
    
    # 2. Network layout visualization
    plt.subplot(2, 3, 3)
    
    # Sample nodes for visualization if graph is large
    if graph.number_of_nodes() > 2000:
        sample_nodes = list(graph.nodes())[:1000]
        subgraph = graph.subgraph(sample_nodes)
        print(f"Visualizing subgraph with {len(sample_nodes)} nodes")
    else:
        subgraph = graph
    
    pos = nx.spring_layout(subgraph, k=1, iterations=50)
    nx.draw(subgraph, pos, 
           node_size=20, 
           node_color='lightcoral',
           edge_color='gray',
           alpha=0.6,
           width=0.5)
    plt.title('Network Structure')
    
    # 3. Component size distribution
    plt.subplot(2, 3, 4)
    component_sizes = [len(c) for c in metrics['components']]
    plt.hist(component_sizes, bins=min(50, len(component_sizes)), alpha=0.7, color='lightgreen')
    plt.xlabel('Component Size')
    plt.ylabel('Frequency')
    plt.title('Component Size Distribution')
    plt.yscale('log')
    
    # 4. Clustering coefficient distribution
    plt.subplot(2, 3, 5)
    local_clustering = nx.clustering(graph)
    clustering_values = list(local_clustering.values())
    plt.hist(clustering_values, bins=50, alpha=0.7, color='orange')
    plt.xlabel('Local Clustering Coefficient')
    plt.ylabel('Frequency')
    plt.title('Clustering Distribution')
    
    # 5. Centrality analysis
    plt.subplot(2, 3, 6)
    if graph.number_of_nodes() <= 5000:  # Only for manageable graphs
        centrality = nx.betweenness_centrality(graph, k=1000)  # Sample for speed
        cent_values = list(centrality.values())
        plt.hist(cent_values, bins=50, alpha=0.7, color='purple')
        plt.xlabel('Betweenness Centrality')
        plt.ylabel('Frequency')
        plt.title('Centrality Distribution')
    else:
        plt.text(0.5, 0.5, 'Graph too large\nfor centrality analysis', 
                ha='center', va='center', transform=plt.gca().transAxes)
        plt.title('Centrality Analysis Skipped')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Comprehensive analysis saved to {output_dir}/comprehensive_analysis.png")

def analyze_synteny_patterns(graph, output_dir="syngraph_plots"):
    """Analyze synteny patterns in the graph"""
    
    print(f"\n=== Synteny Pattern Analysis ===")
    
    # Extract chromosome information from nodes
    chromosome_data = defaultdict(list)
    taxon_data = defaultdict(list)
    
    for node in graph.nodes():
        attrs = graph.nodes[node]
        
        # Look for chromosome/scaffold information
        for attr_name in ['chromosome', 'scaffold', 'sequence', 'seq_id']:
            if attr_name in attrs:
                chromosome_data[attrs[attr_name]].append(node)
                break
        
        # Look for taxon information
        for attr_name in ['taxon', 'species', 'organism']:
            if attr_name in attrs:
                taxon_data[attrs[attr_name]].append(node)
                break
    
    print(f"Chromosomes/scaffolds found: {len(chromosome_data)}")
    print(f"Taxa found: {len(taxon_data)}")
    
    # Visualize chromosome distribution
    if chromosome_data:
        plt.figure(figsize=(12, 8))
        
        chrom_sizes = {chrom: len(nodes) for chrom, nodes in chromosome_data.items()}
        
        plt.subplot(2, 2, 1)
        plt.bar(range(len(chrom_sizes)), sorted(chrom_sizes.values(), reverse=True))
        plt.xlabel('Chromosome Rank')
        plt.ylabel('Number of Genes')
        plt.title('Chromosome Gene Count')
        
        plt.subplot(2, 2, 2)
        plt.hist(chrom_sizes.values(), bins=min(30, len(chrom_sizes)), alpha=0.7, color='lightblue')
        plt.xlabel('Genes per Chromosome')
        plt.ylabel('Frequency')
        plt.title('Gene Count Distribution')
        
        # Top chromosomes
        top_chroms = sorted(chrom_sizes.items(), key=lambda x: x[1], reverse=True)[:10]
        
        plt.subplot(2, 2, 3)
        plt.barh(range(len(top_chroms)), [size for _, size in top_chroms])
        plt.yticks(range(len(top_chroms)), [chrom for chrom, _ in top_chroms])
        plt.xlabel('Number of Genes')
        plt.title('Top 10 Chromosomes')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/synteny_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Synteny analysis saved to {output_dir}/synteny_analysis.png")
    
    return chromosome_data, taxon_data

def main():
    """Main analysis function"""
    
    # Check for required files
    pickle_files = [
        'bibionidae_results.with_ancestors.pickle',
        'bibionidae_analysis.pickle'
    ]
    
    pickle_file = None
    for pf in pickle_files:
        if os.path.exists(pf):
            pickle_file = pf
            break
    
    if pickle_file is None:
        print("ERROR: No syngraph pickle file found!")
        print("Expected files:", pickle_files)
        sys.exit(1)
    
    # Load and analyze
    result = load_and_analyze_syngraph(pickle_file)
    if result is None:
        sys.exit(1)
    
    graph, data = result
    
    # Comprehensive analysis
    metrics = comprehensive_graph_analysis(graph)
    
    # Create visualizations
    create_advanced_visualizations(graph, metrics)
    
    # Analyze synteny
    chromosome_data, taxon_data = analyze_synteny_patterns(graph)
    
    # Save summary report
    with open('syngraph_analysis_report.txt', 'w') as f:
        f.write("SYNGRAPH ANALYSIS REPORT\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input file: {pickle_file}\n")
        f.write(f"Analysis date: {pd.Timestamp.now()}\n\n")
        
        f.write("GRAPH METRICS:\n")
        f.write(f"  Nodes (genes): {metrics['nodes']:,}\n")
        f.write(f"  Edges (adjacencies): {metrics['edges']:,}\n")
        f.write(f"  Average degree: {(2 * metrics['edges']) / metrics['nodes']:.2f}\n")
        f.write(f"  Graph density: {nx.density(graph):.6f}\n")
        f.write(f"  Connected components: {len(metrics['components'])}\n")
        f.write(f"  Average clustering: {metrics['clustering']:.4f}\n\n")
        
        f.write("SYNTENY INFORMATION:\n")
        f.write(f"  Chromosomes/scaffolds: {len(chromosome_data)}\n")
        f.write(f"  Taxa: {len(taxon_data)}\n")
        
        if chromosome_data:
            chrom_sizes = {chrom: len(nodes) for chrom, nodes in chromosome_data.items()}
            f.write(f"  Largest chromosome: {max(chrom_sizes.values())} genes\n")
            f.write(f"  Smallest chromosome: {min(chrom_sizes.values())} genes\n")
            f.write(f"  Average chromosome size: {np.mean(list(chrom_sizes.values())):.1f} genes\n")
    
    print("\n=== Analysis Complete ===")
    print("Generated files:")
    print("  - syngraph_plots/comprehensive_analysis.png")
    print("  - syngraph_plots/synteny_analysis.png")
    print("  - syngraph_analysis_report.txt")
    
    print(f"\nGraph summary:")
    print(f"  {metrics['nodes']:,} genes with {metrics['edges']:,} adjacencies")
    print(f"  {len(metrics['components'])} connected components")
    print(f"  Average clustering: {metrics['clustering']:.4f}")

if __name__ == "__main__":
    main()

