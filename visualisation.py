#!/usr/bin/env python3
"""
Comprehensive Syngraph Visualization Suite
Creates multiple types of visualizations for chromosomal rearrangement analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pickle
import json
import os
from collections import defaultdict, Counter

# Set style
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class SyngraphVisualizer:
    def __init__(self, base_dir="multi_species_analysis"):
        self.base_dir = base_dir
        self.output_dir = f"{base_dir}/visualizations"
        os.makedirs(self.output_dir, exist_ok=True)
        
    def load_syngraph_data(self, results_prefix="diptera_results_clean"):
        """Load all syngraph result files"""
        print(f"Loading syngraph data: {results_prefix}")
        
        # Load rearrangements
        rearr_file = f"{self.base_dir}/{results_prefix}.rearrangements.tsv"
        if os.path.exists(rearr_file):
            self.rearrangements = pd.read_csv(rearr_file, sep='\t')
        else:
            print(f"No rearrangements file found: {rearr_file}")
            self.rearrangements = pd.DataFrame()
        
        # Load pickle file for graph structure
        pickle_file = f"{self.base_dir}/{results_prefix}.with_ancestors.pickle"
        if os.path.exists(pickle_file):
            with open(pickle_file, 'rb') as f:
                self.syngraph_data = pickle.load(f)
        else:
            print(f"No pickle file found: {pickle_file}")
            self.syngraph_data = None
            
        print(f"✓ Loaded {len(self.rearrangements)} rearrangements")
        return self
    
    def create_rearrangement_sankey(self):
        """Create Sankey diagram of rearrangements between species"""
        print("Creating Sankey diagram...")
        
        if len(self.rearrangements) == 0:
            print("No rearrangements to visualize")
            return
        
        # Prepare data for Sankey
        sources = []
        targets = []
        values = []
        labels = set()
        
        for _, rearr in self.rearrangements.iterrows():
            parent = rearr['parent']
            child = rearr['child']
            event = rearr['event']
            mult = rearr.get('multiplicity', 1)
            
            source_label = f"{parent}_{event}"
            target_label = child
            
            labels.add(source_label)
            labels.add(target_label)
            
            sources.append(source_label)
            targets.append(target_label)
            values.append(mult)
        
        # Convert to indices
        label_list = list(labels)
        source_indices = [label_list.index(s) for s in sources]
        target_indices = [label_list.index(t) for t in targets]
        
        # Create Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            node = dict(
                pad = 15,
                thickness = 20,
                line = dict(color = "black", width = 0.5),
                label = label_list,
                color = "blue"
            ),
            link = dict(
                source = source_indices,
                target = target_indices,
                value = values,
                color = "rgba(255, 0, 255, 0.4)"
            )
        )])
        
        fig.update_layout(
            title_text="Chromosomal Rearrangement Flow (Sankey Diagram)", 
            font_size=10,
            width=1200,
            height=800
        )
        
        fig.write_html(f"{self.output_dir}/rearrangement_sankey.html")
        print(f"✓ Sankey diagram saved to: {self.output_dir}/rearrangement_sankey.html")
    
    def create_network_graph_gephi_format(self):
        """Create network graph in Gephi-compatible format"""
        print("Creating network graph for Gephi...")
        
        if self.syngraph_data is None:
            print("No syngraph data available")
            return
        
        # Extract graph structure
        try:
            if hasattr(self.syngraph_data, 'graph'):
                G = self.syngraph_data.graph
            else:
                # Try to extract graph from other attributes
                G = None
                for attr in ['G', 'network', 'syn_graph']:
                    if hasattr(self.syngraph_data, attr):
                        potential_graph = getattr(self.syngraph_data, attr)
                        if hasattr(potential_graph, 'nodes') and hasattr(potential_graph, 'edges'):
                            G = potential_graph
                            break
            
            if G is None:
                print("Could not extract graph structure")
                return
            
            # Create nodes file for Gephi
            nodes_data = []
            for node in G.nodes():
                node_attrs = G.nodes[node] if G.nodes[node] else {}
                nodes_data.append({
                    'Id': str(node),
                    'Label': str(node),
                    'Type': node_attrs.get('type', 'gene'),
                    'Chromosome': node_attrs.get('chromosome', 'unknown'),
                    'Taxon': node_attrs.get('taxon', 'unknown')
                })
            
            nodes_df = pd.DataFrame(nodes_data)
            nodes_df.to_csv(f"{self.output_dir}/syngraph_nodes.csv", index=False)
            
            # Create edges file for Gephi
            edges_data = []
            for i, (source, target) in enumerate(G.edges()):
                edge_attrs = G.edges[source, target] if G.edges[source, target] else {}
                edges_data.append({
                    'Id': i,
                    'Source': str(source),
                    'Target': str(target),
                    'Type': 'Undirected',
                    'Weight': edge_attrs.get('weight', 1),
                    'Adjacency_Type': edge_attrs.get('type', 'synteny')
                })
            
            edges_df = pd.DataFrame(edges_data)
            edges_df.to_csv(f"{self.output_dir}/syngraph_edges.csv", index=False)
            
            print(f"✓ Gephi files created:")
            print(f"  - Nodes: {self.output_dir}/syngraph_nodes.csv ({len(nodes_df)} nodes)")
            print(f"  - Edges: {self.output_dir}/syngraph_edges.csv ({len(edges_df)} edges)")
            
        except Exception as e:
            print(f"Error creating Gephi files: {e}")
    
    def create_cytoscape_network(self):
        """Create network in Cytoscape.js format"""
        print("Creating Cytoscape network...")
        
        if self.syngraph_data is None:
            print("No syngraph data available")
            return
        
        try:
            # Extract graph
            G = None
            for attr in ['graph', 'G', 'network', 'syn_graph']:
                if hasattr(self.syngraph_data, attr):
                    potential_graph = getattr(self.syngraph_data, attr)
                    if hasattr(potential_graph, 'nodes') and hasattr(potential_graph, 'edges'):
                        G = potential_graph
                        break
            
            if G is None:
                print("Could not extract graph structure")
                return
            
            # Sample large graphs for visualization
            if G.number_of_nodes() > 1000:
                print(f"Large graph ({G.number_of_nodes()} nodes), sampling 1000 nodes...")
                sample_nodes = list(G.nodes())[:1000]
                G = G.subgraph(sample_nodes)
            
            # Convert to Cytoscape format
            cytoscape_data = {
                "elements": {
                    "nodes": [],
                    "edges": []
                }
            }
            
            # Add nodes
            for node in G.nodes():
                node_attrs = G.nodes[node] if G.nodes[node] else {}
                cytoscape_data["elements"]["nodes"].append({
                    "data": {
                        "id": str(node),
                        "label": str(node),
                        "type": node_attrs.get('type', 'gene'),
                        "chromosome": node_attrs.get('chromosome', 'unknown'),
                        "taxon": node_attrs.get('taxon', 'unknown')
                    }
                })
            
            # Add edges
            for source, target in G.edges():
                edge_attrs = G.edges[source, target] if G.edges[source, target] else {}
                cytoscape_data["elements"]["edges"].append({
                    "data": {
                        "id": f"{source}-{target}",
                        "source": str(source),
                        "target": str(target),
                        "weight": edge_attrs.get('weight', 1),
                        "type": edge_attrs.get('type', 'synteny')
                    }
                })
            
            # Save JSON
            with open(f"{self.output_dir}/cytoscape_network.json", 'w') as f:
                json.dump(cytoscape_data, f, indent=2)
            
            print(f"✓ Cytoscape network saved: {self.output_dir}/cytoscape_network.json")
            print(f"  - {len(cytoscape_data['elements']['nodes'])} nodes")
            print(f"  - {len(cytoscape_data['elements']['edges'])} edges")
            
        except Exception as e:
            print(f"Error creating Cytoscape network: {e}")
    
    def create_interactive_phylogeny_plot(self):
        """Create interactive phylogenetic rearrangement plot"""
        print("Creating interactive phylogeny plot...")
        
        if len(self.rearrangements) == 0:
            print("No rearrangements to plot")
            return
        
        # Group rearrangements by branch
        branch_data = defaultdict(list)
        for _, rearr in self.rearrangements.iterrows():
            branch = f"{rearr['parent']} → {rearr['child']}"
            branch_data[branch].append({
                'event': rearr['event'],
                'multiplicity': rearr.get('multiplicity', 1)
            })
        
        # Create interactive plot
        branches = list(branch_data.keys())
        event_counts = [len(branch_data[branch]) for branch in branches]
        event_types = [', '.join([e['event'] for e in branch_data[branch]]) for branch in branches]
        
        fig = go.Figure(data=[
            go.Bar(
                x=branches,
                y=event_counts,
                text=event_types,
                textposition='auto',
                marker=dict(
                    color=event_counts,
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title="Number of Events")
                )
            )
        ])
        
        fig.update_layout(
            title="Chromosomal Rearrangements by Phylogenetic Branch",
            xaxis_title="Phylogenetic Branch",
            yaxis_title="Number of Rearrangement Events",
            xaxis_tickangle=-45,
            width=1200,
            height=600
        )
        
        fig.write_html(f"{self.output_dir}/phylogeny_rearrangements.html")
        print(f"✓ Interactive phylogeny plot saved: {self.output_dir}/phylogeny_rearrangements.html")
    
    def create_chromosome_evolution_plot(self):
        """Create chromosome number evolution visualization"""
        print("Creating chromosome evolution plot...")
        
        # This would require parsing the syngraph data for chromosome numbers
        # For now, create a placeholder with example data
        
        species = ['Bibio_marci', 'Dilophus_febrilis', 'Drosophila_melanogaster', 
                  'Aedes_aegypti', 'Calliphora_vicina']
        chromosomes = [6, 6, 4, 3, 6]  # Example data
        
        fig = go.Figure(data=[
            go.Scatter(
                x=species,
                y=chromosomes,
                mode='markers+lines',
                marker=dict(size=12, color='blue'),
                line=dict(width=2),
                name='Chromosome Count'
            )
        ])
        
        fig.update_layout(
            title="Chromosome Number Evolution Across Diptera",
            xaxis_title="Species",
            yaxis_title="Chromosome Count",
            xaxis_tickangle=-45
        )
        
        fig.write_html(f"{self.output_dir}/chromosome_evolution.html")
        print(f"✓ Chromosome evolution plot saved: {self.output_dir}/chromosome_evolution.html")
    
    def create_summary_dashboard(self):
        """Create comprehensive summary dashboard"""
        print("Creating summary dashboard...")
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Rearrangement Types', 'Events by Branch', 
                           'Event Multiplicity', 'Timeline'),
            specs=[[{"type": "pie"}, {"type": "bar"}],
                   [{"type": "histogram"}, {"type": "scatter"}]]
        )
        
        if len(self.rearrangements) > 0:
            # Plot 1: Rearrangement types pie chart
            event_counts = self.rearrangements['event'].value_counts()
            fig.add_trace(
                go.Pie(labels=event_counts.index, values=event_counts.values, name="Events"),
                row=1, col=1
            )
            
            # Plot 2: Events by branch
            branch_counts = self.rearrangements.groupby(['parent', 'child']).size()
            fig.add_trace(
                go.Bar(x=[f"{p}→{c}" for p, c in branch_counts.index], 
                      y=branch_counts.values, name="Branch Events"),
                row=1, col=2
            )
            
            # Plot 3: Multiplicity histogram
            if 'multiplicity' in self.rearrangements.columns:
                fig.add_trace(
                    go.Histogram(x=self.rearrangements['multiplicity'], name="Multiplicity"),
                    row=2, col=1
                )
        
        fig.update_layout(
            title_text="Syngraph Analysis Dashboard",
            showlegend=False,
            height=800,
            width=1200
        )
        
        fig.write_html(f"{self.output_dir}/summary_dashboard.html")
        print(f"✓ Summary dashboard saved: {self.output_dir}/summary_dashboard.html")
    
    def generate_all_visualizations(self, results_prefix="diptera_results_clean"):
        """Generate all visualization types"""
        print("=== Generating Comprehensive Visualization Suite ===")
        
        # Load data
        self.load_syngraph_data(results_prefix)
        
        # Generate all visualizations
        self.create_rearrangement_sankey()
        self.create_network_graph_gephi_format()
        self.create_cytoscape_network()
        self.create_interactive_phylogeny_plot()
        self.create_chromosome_evolution_plot()
        self.create_summary_dashboard()
        
        print(f"\n✓ All visualizations complete!")
        print(f"📁 Output directory: {self.output_dir}")
        print(f"📊 Files created:")
        
        for file in os.listdir(self.output_dir):
            file_path = os.path.join(self.output_dir, file)
            size_mb = os.path.getsize(file_path) / 1024 / 1024
            print(f"  - {file} ({size_mb:.2f} MB)")

def main():
    """Main execution function"""
    
    # Create visualizer
    viz = SyngraphVisualizer()
    
    # Generate for both result sets if available
    result_sets = ["diptera_results_clean", "diptera_results_simple"]
    
    for result_set in result_sets:
        rearr_file = f"multi_species_analysis/{result_set}.rearrangements.tsv"
        if os.path.exists(rearr_file):
            print(f"\n=== Processing {result_set} ===")
            viz.generate_all_visualizations(result_set)
        else:
            print(f"Skipping {result_set} (file not found)")
    
    print(f"\n🎉 Visualization suite complete!")
    print(f"Open the HTML files in your browser to view interactive plots")

if __name__ == "__main__":
    main()