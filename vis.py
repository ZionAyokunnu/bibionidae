# # Simple Sankey visualization
# import pandas as pd
# import plotly.graph_objects as go
# import os

# # Load data
# df = pd.read_csv('multi_species_analysis/diptera_results_clean.rearrangements.tsv', sep='\t')

# # Create output directory
# os.makedirs('multi_species_analysis/visualizations', exist_ok=True)

# # Use actual column names
# cols = df.columns.tolist()
# parent_col = cols[0]  # Usually '#parent'
# child_col = cols[1]   # Usually 'child'
# event_col = cols[2]   # Usually 'event'

# # Create Sankey
# sources = []
# targets = []
# values = []
# labels = set()

# for _, row in df.iterrows():
#     parent = str(row[parent_col])
#     child = str(row[child_col])
#     event = str(row[event_col])
    
#     source_label = f"{parent}_{event}"
#     target_label = child
    
#     labels.add(source_label)
#     labels.add(target_label)
    
#     sources.append(source_label)
#     targets.append(target_label)
#     values.append(1)

# # Convert to indices
# label_list = list(labels)
# source_indices = [label_list.index(s) for s in sources]
# target_indices = [label_list.index(t) for t in targets]

# # Create figure
# fig = go.Figure(data=[go.Sankey(
#     node=dict(
#         pad=15,
#         thickness=20,
#         line=dict(color="black", width=0.5),
#         label=label_list,
#         color="blue"
#     ),
#     link=dict(
#         source=source_indices,
#         target=target_indices,
#         value=values,
#         color="rgba(255, 0, 255, 0.4)"
#     )
# )])

# fig.update_layout(
#     title_text="Chromosomal Rearrangements (Sankey)", 
#     font_size=10,
#     width=1200,
#     height=800
# )

# fig.write_html('multi_species_analysis/visualizations/sankey.html')
# print("✓ Sankey diagram saved to multi_species_analysis/visualizations/sankey.html")


# # Create network graph for Gephi
# import pandas as pd
# import pickle
# import sys
# import os

# # Add syngraph to path
# sys.path.append('syngraph')

# # Load pickle data
# with open('multi_species_analysis/diptera_results_clean.with_ancestors.pickle', 'rb') as f:
#     data = pickle.load(f)

# # Extract graph
# G = data.graph

# # Create nodes CSV
# nodes_data = []
# for node in G.nodes():
#     node_attrs = G.nodes[node] if G.nodes[node] else {}
#     nodes_data.append({
#         'Id': str(node),
#         'Label': str(node),
#         'Type': node_attrs.get('type', 'gene'),
#         'Chromosome': node_attrs.get('chromosome', 'unknown'),
#         'Taxon': node_attrs.get('taxon', 'unknown')
#     })

# nodes_df = pd.DataFrame(nodes_data)
# nodes_df.to_csv('multi_species_analysis/visualizations/nodes.csv', index=False)

# # Create edges CSV
# edges_data = []
# for i, (source, target) in enumerate(G.edges()):
#     edge_attrs = G.edges[source, target] if G.edges[source, target] else {}
#     edges_data.append({
#         'Id': i,
#         'Source': str(source),
#         'Target': str(target),
#         'Type': 'Undirected',
#         'Weight': edge_attrs.get('weight', 1)
#     })

# edges_df = pd.DataFrame(edges_data)
# edges_df.to_csv('multi_species_analysis/visualizations/edges.csv', index=False)

# print(f"✓ Network files created:")
# print(f"  - Nodes: {len(nodes_df)} nodes")
# print(f"  - Edges: {len(edges_df)} edges")





# # Debug the pickle data structure
# import pickle
# import sys
# sys.path.append('syngraph')

# # Load and examine the data
# with open('multi_species_analysis/diptera_results_clean.with_ancestors.pickle', 'rb') as f:
#    data = pickle.load(f)

# print("Data type:", type(data))
# print("Keys:", list(data.keys()) if hasattr(data, 'keys') else 'No keys')

# # Try to find the graph in different locations
# if hasattr(data, 'graph'):
#    print("Found graph at data.graph")
#    G = data.graph
# elif 'graph' in data:
#    print("Found graph at data['graph']")
#    G = data['graph']
# else:
#    print("Available attributes:", dir(data))
#    # Try common graph attribute names
#    for attr in ['G', 'network', 'syn_graph', '_graph']:
#        if hasattr(data, attr):
#            print(f"Found potential graph at: {attr}")
#            G = getattr(data, attr)
#            break
#    else:
#        print("No graph found")
#        G = None

# if G and hasattr(G, 'nodes'):
#    print(f"Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
# else:
#    print("No valid graph object found")





# # Examine the Syngraph object structure
# import pickle
# import sys
# sys.path.append('syngraph')

# with open('multi_species_analysis/diptera_results_clean.with_ancestors.pickle', 'rb') as f:
#     data = pickle.load(f)

# print("Syngraph attributes:")
# for attr in dir(data):
#     if not attr.startswith('_'):
#         print(f"  {attr}: {type(getattr(data, attr))}")

# # Check the graph object
# G = data.graph
# print(f"\nGraph type: {type(G)}")
# print(f"Graph attributes:")
# for attr in dir(G):
#     if not attr.startswith('_'):
#         print(f"  {attr}: {type(getattr(G, attr))}")

# # Try to access graph data
# try:
#     print(f"\nGraph info:")
#     print(f"  Nodes: {len(G.nodes()) if hasattr(G, 'nodes') else 'No nodes method'}")
#     print(f"  Edges: {len(G.edges()) if hasattr(G, 'edges') else 'No edges method'}")
# except Exception as e:
#     print(f"Error accessing graph: {e}")






# # Create network graph for Gephi - corrected version
# import pandas as pd
# import pickle
# import sys
# import os

# # Add syngraph to path
# sys.path.append('syngraph')

# # Load pickle data
# with open('multi_species_analysis/diptera_results_clean.with_ancestors.pickle', 'rb') as f:
#     data = pickle.load(f)

# # The data object itself is the graph
# G = data

# print(f"Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")

# # Create nodes CSV
# nodes_data = []
# for node in G.nodes():
#     node_attrs = G.nodes[node] if G.nodes[node] else {}
#     nodes_data.append({
#         'Id': str(node),
#         'Label': str(node),
#         'Type': node_attrs.get('type', 'gene'),
#         'Chromosome': node_attrs.get('chromosome', 'unknown'),
#         'Taxon': node_attrs.get('taxon', 'unknown')
#     })

# nodes_df = pd.DataFrame(nodes_data)
# nodes_df.to_csv('multi_species_analysis/visualizations/nodes.csv', index=False)

# # Create edges CSV
# edges_data = []
# for i, (source, target) in enumerate(G.edges()):
#     edge_attrs = G.edges[source, target] if G.edges[source, target] else {}
#     edges_data.append({
#         'Id': i,
#         'Source': str(source),
#         'Target': str(target),
#         'Type': 'Undirected',
#         'Weight': edge_attrs.get('weight', 1)
#     })

# edges_df = pd.DataFrame(edges_data)
# edges_df.to_csv('multi_species_analysis/visualizations/edges.csv', index=False)

# print(f"✓ Network files created:")
# print(f"  - Nodes: {len(nodes_df)} nodes")
# print(f"  - Edges: {len(edges_df)} edges")







# Create interactive rearrangement plot
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

# Load rearrangements data
df = pd.read_csv('multi_species_analysis/diptera_results_clean.rearrangements.tsv', sep='\t')

# Create bar chart of rearrangements by event type
event_counts = df['event'].value_counts()

fig = go.Figure(data=[
    go.Bar(
        x=event_counts.index,
        y=event_counts.values,
        marker_color='steelblue',
        text=event_counts.values,
        textposition='auto'
    )
])

fig.update_layout(
    title="Chromosomal Rearrangement Events",
    xaxis_title="Event Type",
    yaxis_title="Count",
    width=800,
    height=600
)

fig.write_html('multi_species_analysis/visualizations/events_bar.html')
print("✓ Events bar chart saved")

# Create branch-specific plot
branch_data = df.groupby(['#parent', 'child', 'event']).size().reset_index(name='count')

fig2 = px.sunburst(
    branch_data,
    path=['#parent', 'child', 'event'],
    values='count',
    title="Rearrangements by Phylogenetic Branch"
)

fig2.write_html('multi_species_analysis/visualizations/branch_sunburst.html')
print("✓ Branch sunburst plot saved")

# Show multiplicity distribution
fig3 = go.Figure(data=[
    go.Histogram(
        x=df['multiplicity'],
        nbinsx=20,
        marker_color='lightcoral'
    )
])

fig3.update_layout(
    title="Distribution of Rearrangement Multiplicity",
    xaxis_title="Multiplicity",
    yaxis_title="Frequency",
    width=800,
    height=600
)

fig3.write_html('multi_species_analysis/visualizations/multiplicity_hist.html')
print("✓ Multiplicity histogram saved")

print("\nFiles created:")
print("  - multi_species_analysis/visualizations/sankey.html")
print("  - multi_species_analysis/visualizations/nodes.csv")
print("  - multi_species_analysis/visualizations/edges.csv") 
print("  - multi_species_analysis/visualizations/events_bar.html")
print("  - multi_species_analysis/visualizations/branch_sunburst.html")
print("  - multi_species_analysis/visualizations/multiplicity_hist.html")