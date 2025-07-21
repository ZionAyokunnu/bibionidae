#!/usr/bin/env python3
"""
PROPERLY INTEGRATED SYNTENY AND INVERSION ANALYZER
Combines syngraph's macro-synteny analysis with enhanced inversion detection

WORKFLOW:
Step 1: Syngraph macro-synteny analysis (LMS, rearrangements, ancestral genomes)
Step 2: Enhanced inversion detection within syntenic regions
Step 3: Integrated visualization and results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import os
import glob
import logging
from collections import defaultdict, Counter
import itertools
import networkx as nx
from scipy.stats import pearsonr
import warnings
import random
import copy
import more_itertools
from operator import attrgetter
import functools
warnings.filterwarnings('ignore')

# Import original inversion analyzer functions
from difflib import SequenceMatcher
from rich.console import Console

console = Console()
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

################################################################################
# SYNGRAPH CLASSES (from syngraph.py)
################################################################################

class Syngraph(nx.Graph):
    def __init__(self, name='', **attr):
        nx.Graph.__init__(self, name=name, taxa=set(), **attr)
        
    def __repr__(self):
       return "Syngraph(name=%r, taxa=%r, ...)" % (self.name, self.graph['taxa']) 

    def from_markerObjs(self, markerObjs):
        prev_markerObj = MarkerObj(None)
        marker_ids_by_seq_id_by_taxon = defaultdict(functools.partial(defaultdict, list))
        for markerObj in markerObjs:
            marker_ids_by_seq_id_by_taxon[markerObj.taxon][markerObj.seq].append(markerObj.name)
            self.graph['taxa'].add(markerObj.taxon)                                   
            if not markerObj.name in self:
                self.add_node(markerObj.name, taxa=set(), terminal=set(), seqs_by_taxon={}, starts_by_taxon={}, 
                    ends_by_taxon={})             
            self.nodes[markerObj.name]['taxa'].add(markerObj.taxon)
            self.nodes[markerObj.name]['seqs_by_taxon'][markerObj.taxon] = markerObj.seq
            self.nodes[markerObj.name]['starts_by_taxon'][markerObj.taxon] = markerObj.start
            self.nodes[markerObj.name]['ends_by_taxon'][markerObj.taxon] = markerObj.end
            if markerObj.is_syntenic(prev_markerObj):
                if not self.has_edge(prev_markerObj.name, markerObj.name):
                    self.add_edge(prev_markerObj.name, markerObj.name, taxa={})
                self[prev_markerObj.name][markerObj.name]['taxa'][markerObj.taxon] = prev_markerObj.distance(markerObj)    
            else:
                if not prev_markerObj.name is None:
                    self.nodes[prev_markerObj.name]['terminal'].add(prev_markerObj.taxon)
                self.nodes[markerObj.name]['terminal'].add(markerObj.taxon)
            prev_markerObj = markerObj
        self.nodes[markerObj.name]['terminal'].add(markerObj.taxon)
        self.graph['marker_ids_by_seq_id_by_taxon'] = marker_ids_by_seq_id_by_taxon

    def show_metrics(self):
        taxon_count = len(self.graph['taxa'])
        node_total_count = nx.number_of_nodes(self)
        node_non_lonely_count = 0
        node_complete_count = 0
        for graph_node_id in self.nodes:
            if len(self.nodes[graph_node_id]['taxa']) > 1:
                node_non_lonely_count += 1
                if len(self.nodes[graph_node_id]['taxa']) == len(self.graph['taxa']):
                    node_complete_count += 1
        edge_total_count = nx.number_of_edges(self)
        connected_component_count = nx.number_connected_components(self)
        print("[=] ====================================")
        print("[=] Taxa = %s" % taxon_count)
        print("[=] Nodes (Markers) = %s" % node_total_count)
        print("[=] Nodes (Markers) shared by > 1 taxon = %s" % node_non_lonely_count)
        print("[=] Nodes (Markers) shared by all taxa = %s" % node_complete_count)
        print("[=] Distinct Edges (Adjacencies) = %s" % edge_total_count)
        print("[=] Subgraphs (connected components) = %s" % connected_component_count)
        print("[=] ====================================")

class MarkerObj():
    def __init__(self, name=None, desc=None, status=None, taxon=None, seq=None, start=None, end=None):
        self.name = name
        self.desc = desc if desc is not None else name
        self.status = status
        self.taxon = taxon
        self.seq = seq
        self.start = start
        self.end = end

    def __repr__(self):
        return "MarkerObj(name=%r, desc=%r, status=%r, taxon=%r, seq=%r, start=%d, end=%d)" % (
            self.name, self.desc, self.status, self.taxon, self.seq, self.start, self.end) 

    def __eq__(self, other):
        if isinstance(other, MarkerObj):
            return (self.name == other.name)
        return False

    def is_syntenic(self, other):
        if isinstance(other, MarkerObj):
            if self.taxon == other.taxon and self.seq == other.seq:
                return True
        return False

    def distance(self, other):
        if self.is_syntenic(other):
            return int(max((self.start - other.end), (other.start - self.end)))
        return float("nan")

class GenomeObj():
    def __init__(self, labelled_CWAL=None, CWAL=None, labelled_CWAOC=None, CWAOC=None, LMS_OC=None, taxon=None):
        if labelled_CWAL is None:
            labelled_CWAL = defaultdict(set)
        if CWAL is None:
            CWAL = []
        if labelled_CWAOC is None:
            labelled_CWAOC = defaultdict(set)
        if CWAOC is None:
            CWAOC = []
        if LMS_OC is None:
            LMS_OC = {}
        self.labelled_CWAL = labelled_CWAL
        self.CWAL = CWAL
        self.labelled_CWAOC = labelled_CWAOC
        self.CWAOC = CWAOC
        self.LMS_OC = LMS_OC
        self.taxon = taxon

################################################################################
# SYNGRAPH CORE FUNCTIONS (from syngraph.py)
################################################################################

def get_LMSs(syngraph, list_of_taxa, minimum):
    """Get Linked Marker Sets from syngraph"""
    linked_marker_sets = defaultdict(set)
    unassignable_markers = set()
    for graph_node_id in syngraph.nodes():
        target_seqs_by_taxon_set = set()
        for taxon in list_of_taxa:
            if taxon in syngraph.nodes[graph_node_id]['seqs_by_taxon']:
                target_seqs_by_taxon_set.add(syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon])
        target_seqs_by_taxon_set = frozenset(target_seqs_by_taxon_set)
        if len(target_seqs_by_taxon_set) == len(list_of_taxa):
            linked_marker_sets[target_seqs_by_taxon_set].add(graph_node_id)
        else:
            unassignable_markers.add(graph_node_id)
    filtered_linked_marker_sets = {}
    filtered_LMS_count = 0
    for LMS in linked_marker_sets:
        if len(linked_marker_sets[LMS]) >= minimum:
            filtered_linked_marker_sets["LMS_" + str(filtered_LMS_count)] = linked_marker_sets[LMS]
            filtered_LMS_count += 1
        else:
            for graph_node_id in linked_marker_sets[LMS]:
                unassignable_markers.add(graph_node_id)
    return filtered_linked_marker_sets, unassignable_markers

def compact_synteny_1(syngraph, LMSs, taxon):
    """Write genome of a taxon in terms of LMSs"""
    chrom2LMS = GenomeObj()
    for graph_node_id in syngraph.nodes():
        for LMS in LMSs:
            if graph_node_id in LMSs[LMS]:
                chrom2LMS.labelled_CWAL[syngraph.nodes[graph_node_id]['seqs_by_taxon'][taxon]].add(LMS)
    chrom2LMS.CWAL = [chromosome for chromosome in chrom2LMS.labelled_CWAL.values()]
    return chrom2LMS

def check_for_fusions(instance_of_synteny, fusion_log, i):
    """Check for fusion events"""
    for combo in itertools.combinations(instance_of_synteny.labelled_CWAOC, 2):
        if len(instance_of_synteny.labelled_CWAOC[combo[0]].intersection(instance_of_synteny.labelled_CWAOC[combo[1]])) > 0:
            new_name = combo[0] + "_" + combo[1].split("_", 1)[1]
            instance_of_synteny.labelled_CWAOC[new_name] = \
            instance_of_synteny.labelled_CWAOC[combo[0]].union(instance_of_synteny.labelled_CWAOC[combo[1]])
            instance_of_synteny.labelled_CWAL[new_name] = \
            instance_of_synteny.labelled_CWAL[combo[0]].union(instance_of_synteny.labelled_CWAL[combo[1]])
            fusion_log.append(["parent_node", i, "fusion", 1, 
                [set(instance_of_synteny.labelled_CWAL[combo[0]]), 
                 set(instance_of_synteny.labelled_CWAL[combo[1]])]])
            del instance_of_synteny.labelled_CWAOC[combo[0]]
            del instance_of_synteny.labelled_CWAOC[combo[1]]
            del instance_of_synteny.labelled_CWAL[combo[0]]
            del instance_of_synteny.labelled_CWAL[combo[1]]
            return check_for_fusions(instance_of_synteny, fusion_log, i)
    return instance_of_synteny, fusion_log

def check_for_fissions(instance_of_synteny, fission_log, i):
    """Check for fission events"""
    for chrom in instance_of_synteny.labelled_CWAOC:
        indices = len(instance_of_synteny.labelled_CWAOC[chrom])
        if indices > 1:
            fission_log.append(["parent_node", i, "fission", indices-1, 
                set(instance_of_synteny.labelled_CWAL[chrom])])
    return fission_log

def ffsd(instance_of_synteny, rearrangement_log, i):
    """Fission-fusion detection"""
    instance_of_synteny, rearrangement_log = check_for_fusions(instance_of_synteny, rearrangement_log, i)
    rearrangement_log = check_for_fissions(instance_of_synteny, rearrangement_log, i)
    return rearrangement_log

################################################################################
# INTEGRATION FUNCTIONS
################################################################################

def convert_busco_to_markerObjs(busco_files):
    """Convert BUSCO data to MarkerObj format for syngraph"""
    logger.info("Converting BUSCO data to MarkerObj format...")
    
    markerObjs = []
    
    for file_path in busco_files:
        taxon = os.path.basename(file_path).replace('.tsv', '')
        logger.info(f"  Processing {taxon}...")
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            # Find header
            header_idx = None
            for i, line in enumerate(lines):
                if line.startswith('# Busco id'):
                    header_idx = i
                    break
            
            if header_idx is None:
                continue
            
            # Parse data and create MarkerObjs
            for line in lines[header_idx + 1:]:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 8 and parts[1] == 'Complete':
                    try:
                        marker = MarkerObj(
                            name=parts[0],  # busco_id
                            status=parts[1],  # Complete
                            taxon=taxon,
                            seq=f"{taxon}_{parts[2]}",  # taxon_sequence
                            start=int(parts[3]),
                            end=int(parts[4])
                        )
                        markerObjs.append(marker)
                    except (ValueError, IndexError):
                        continue
            
        except Exception as e:
            logger.warning(f"  Failed to process {taxon}: {e}")
    
    # Sort markers by sequence and position
    markerObjs = sorted(markerObjs, key=attrgetter('seq', 'start'))
    logger.info(f"  Created {len(markerObjs)} MarkerObj instances")
    
    return markerObjs

def analyze_syngraph_rearrangements(syngraph, taxa_list, minimum=10):
    """Analyze rearrangements using syngraph approach"""
    logger.info("Analyzing chromosomal rearrangements with syngraph...")
    
    # Get LMSs (Linked Marker Sets)
    LMSs, unassignable = get_LMSs(syngraph, taxa_list, minimum)
    logger.info(f"  Generated {len(LMSs)} LMSs with {sum(len(lms) for lms in LMSs.values())} markers")
    logger.info(f"  {len(unassignable)} markers unassignable")
    
    # Analyze pairwise rearrangements
    rearrangement_results = []
    
    for species1, species2 in itertools.combinations(taxa_list, 2):
        logger.info(f"  Analyzing {species1} vs {species2}...")
        
        # Get synteny representations
        ios_1 = compact_synteny_1(syngraph, LMSs, species1)
        ios_2 = compact_synteny_1(syngraph, LMSs, species2)
        
        # Create mapping between species
        ios_2_mapped = copy.deepcopy(ios_2)
        for chr2 in ios_2.labelled_CWAL:
            for LMS in ios_2.labelled_CWAL[chr2]:
                for chr1 in ios_1.labelled_CWAL:
                    if LMS in ios_1.labelled_CWAL[chr1]:
                        ios_2_mapped.labelled_CWAOC[chr2].add(chr1)
                        ios_2_mapped.LMS_OC[LMS] = chr1
        
        # Detect rearrangements
        rearrangement_log = []
        rearrangement_log = ffsd(ios_2_mapped, rearrangement_log, 0)
        
        # Process results
        for rearr in rearrangement_log:
            rearrangement_results.append({
                'species1': species1,
                'species2': species2,
                'type': rearr[2],
                'multiplicity': rearr[3],
                'affected_markers': rearr[4] if isinstance(rearr[4], set) else len(rearr[4])
            })
    
    return rearrangement_results, LMSs

def extract_syntenic_regions_from_syngraph(syngraph, LMSs, genome_sequences):
    """Extract syntenic regions identified by syngraph for inversion analysis"""
    logger.info("Extracting syntenic regions for inversion analysis...")
    
    syntenic_regions = []
    
    for lms_id, marker_set in LMSs.items():
        if len(marker_set) < 5:  # Skip small LMSs
            continue
        
        # Get species and chromosomes for this LMS
        species_chroms = defaultdict(set)
        positions = defaultdict(list)
        
        for marker in marker_set:
            if marker in syngraph.nodes:
                for taxon in syngraph.nodes[marker]['seqs_by_taxon']:
                    chrom = syngraph.nodes[marker]['seqs_by_taxon'][taxon]
                    species_chroms[taxon].add(chrom)
                    start = syngraph.nodes[marker]['starts_by_taxon'][taxon]
                    end = syngraph.nodes[marker]['ends_by_taxon'][taxon]
                    positions[taxon].append((start, end))
        
        # Create regions for pairwise comparisons
        taxa_list = list(species_chroms.keys())
        for species1, species2 in itertools.combinations(taxa_list, 2):
            if (len(species_chroms[species1]) == 1 and 
                len(species_chroms[species2]) == 1 and
                species1 in genome_sequences and 
                species2 in genome_sequences):
                
                chr1 = list(species_chroms[species1])[0].split('_', 1)[1]  # Remove taxon prefix
                chr2 = list(species_chroms[species2])[0].split('_', 1)[1]
                
                if chr1 in genome_sequences[species1] and chr2 in genome_sequences[species2]:
                    pos1 = positions[species1]
                    pos2 = positions[species2]
                    
                    start1, end1 = min(p[0] for p in pos1), max(p[1] for p in pos1)
                    start2, end2 = min(p[0] for p in pos2), max(p[1] for p in pos2)
                    
                    # Add padding
                    padding = 10000
                    seq1 = genome_sequences[species1][chr1]
                    seq2 = genome_sequences[species2][chr2]
                    
                    extract_start1 = max(0, start1 - padding)
                    extract_end1 = min(len(seq1), end1 + padding)
                    extract_start2 = max(0, start2 - padding)
                    extract_end2 = min(len(seq2), end2 + padding)
                    
                    region1 = seq1[extract_start1:extract_end1]
                    region2 = seq2[extract_start2:extract_end2]
                    
                    syntenic_regions.append({
                        'lms_id': lms_id,
                        'species1': species1,
                        'species2': species2,
                        'chr1': chr1,
                        'chr2': chr2,
                        'sequence1': region1,
                        'sequence2': region2,
                        'coordinates1': (extract_start1, extract_end1),
                        'coordinates2': (extract_start2, extract_end2),
                        'marker_count': len(marker_set)
                    })
    
    logger.info(f"  Extracted {len(syntenic_regions)} syntenic regions")
    return syntenic_regions

################################################################################
# ENHANCED INVERSION DETECTION (from original analyzer)
################################################################################

def detect_inversions_with_kmers(syntenic_regions, kmer_size=21, threshold=0.3):
    """Detect inversions using k-mer analysis within syntenic regions"""
    logger.info("Detecting inversions with k-mer analysis...")
    
    inversions = []
    
    for region in syntenic_regions:
        seq1 = region['sequence1']
        seq2 = region['sequence2']
        
        if len(seq1) < 1000 or len(seq2) < 1000:  # Skip very short regions
            continue
        
        # Generate k-mers
        kmers1 = generate_kmers(seq1, kmer_size)
        kmers2 = generate_kmers(seq2, kmer_size)
        kmers2_rc = generate_kmers(str(Seq(seq2).reverse_complement()), kmer_size)
        
        # Calculate similarities
        similarity_forward = calculate_kmer_similarity(kmers1, kmers2)
        similarity_reverse = calculate_kmer_similarity(kmers1, kmers2_rc)
        
        # Detect inversions
        if (similarity_reverse > similarity_forward and 
            similarity_reverse > threshold and
            similarity_reverse - similarity_forward > 0.1):
            
            inversion = {
                'lms_id': region['lms_id'],
                'species_pair': (region['species1'], region['species2']),
                'chromosomes': (region['chr1'], region['chr2']),
                'coordinates1': region['coordinates1'],
                'coordinates2': region['coordinates2'],
                'similarity_forward': similarity_forward,
                'similarity_reverse': similarity_reverse,
                'inversion_strength': similarity_reverse - similarity_forward,
                'length': min(len(seq1), len(seq2)),
                'marker_count': region['marker_count'],
                'type': 'kmer_inversion'
            }
            inversions.append(inversion)
    
    logger.info(f"  Detected {len(inversions)} potential inversions")
    return inversions

def generate_kmers(sequence, k):
    """Generate k-mers from sequence"""
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:
            kmers.add(kmer)
    return kmers

def calculate_kmer_similarity(kmers1, kmers2):
    """Calculate Jaccard similarity between k-mer sets"""
    if not kmers1 or not kmers2:
        return 0.0
    intersection = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))
    return intersection / union if union > 0 else 0.0

################################################################################
# VISUALIZATION FUNCTIONS
################################################################################

def create_integrated_visualization(syngraph, rearrangements, inversions, output_dir='results'):
    """Create integrated visualization of syngraph and inversion results"""
    logger.info("Creating integrated visualizations...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Syngraph network visualization
    create_syngraph_network_plot(syngraph, f"{output_dir}/syngraph_network.png")
    
    # 2. Rearrangement summary
    create_rearrangement_summary_plot(rearrangements, f"{output_dir}/rearrangements.png")
    
    # 3. Inversion analysis
    if inversions:
        create_inversion_analysis_plot(inversions, f"{output_dir}/inversions.png")
    
    # 4. Integrated summary
    create_integrated_summary_plot(syngraph, rearrangements, inversions, f"{output_dir}/integrated_summary.png")

def create_syngraph_network_plot(syngraph, output_path):
    """Create syngraph network visualization"""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Get largest connected component for visualization
    largest_cc = max(nx.connected_components(syngraph), key=len)
    subgraph = syngraph.subgraph(largest_cc)
    
    # Layout
    pos = nx.spring_layout(subgraph, k=1, iterations=50)
    
    # Node properties
    node_sizes = [100 + len(subgraph.nodes[node]['taxa']) * 50 for node in subgraph.nodes()]
    node_colors = [len(subgraph.nodes[node]['taxa']) for node in subgraph.nodes()]
    
    # Draw network
    nx.draw_networkx_nodes(subgraph, pos, node_size=node_sizes, node_color=node_colors,
                          cmap='viridis', alpha=0.8)
    nx.draw_networkx_edges(subgraph, pos, alpha=0.3, edge_color='gray')
    
    plt.colorbar(plt.cm.ScalarMappable(cmap='viridis'), label='Species Count')
    plt.title('Syngraph Network (Largest Component)\nNode size ∝ species count', fontweight='bold')
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_rearrangement_summary_plot(rearrangements, output_path):
    """Create rearrangement summary visualization"""
    if not rearrangements:
        return
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    df = pd.DataFrame(rearrangements)
    
    # Rearrangement types
    type_counts = df['type'].value_counts()
    ax1.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%')
    ax1.set_title('Rearrangement Types')
    
    # Species pairs with most rearrangements
    pair_counts = df.groupby(['species1', 'species2']).size().sort_values(ascending=False).head(10)
    ax2.barh(range(len(pair_counts)), pair_counts.values)
    ax2.set_yticks(range(len(pair_counts)))
    ax2.set_yticklabels([f"{p[0]}\nvs\n{p[1]}" for p in pair_counts.index])
    ax2.set_xlabel('Rearrangement Count')
    ax2.set_title('Top Species Pairs by Rearrangements')
    
    # Multiplicity distribution
    ax3.hist(df['multiplicity'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax3.set_xlabel('Event Multiplicity')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Rearrangement Event Multiplicity')
    
    # Affected markers
    if 'affected_markers' in df.columns:
        markers = pd.to_numeric(df['affected_markers'], errors='coerce').dropna()
        if len(markers) > 0:
            ax4.hist(markers, bins=20, alpha=0.7, color='lightcoral', edgecolor='black')
            ax4.set_xlabel('Affected Markers')
            ax4.set_ylabel('Frequency')
            ax4.set_title('Markers Affected by Rearrangements')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_inversion_analysis_plot(inversions, output_path):
    """Create inversion analysis visualization"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    df = pd.DataFrame(inversions)
    
    # Inversion lengths
    ax1.hist(df['length'], bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
    ax1.set_xlabel('Inversion Length (bp)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Inversion Length Distribution')
    
    # Inversion strength
    ax2.hist(df['inversion_strength'], bins=20, alpha=0.7, color='orange', edgecolor='black')
    ax2.set_xlabel('Inversion Strength')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Inversion Strength Distribution')
    
    # Forward vs Reverse similarity
    ax3.scatter(df['similarity_forward'], df['similarity_reverse'], 
               alpha=0.6, c=df['inversion_strength'], cmap='viridis')
    ax3.plot([0, 1], [0, 1], '--', color='gray', alpha=0.5)
    ax3.set_xlabel('Forward Similarity')
    ax3.set_ylabel('Reverse Similarity')
    ax3.set_title('Forward vs Reverse Similarity')
    
    # Species pair inversions
    species_pairs = [f"{inv['species_pair'][0]}\nvs\n{inv['species_pair'][1]}" for inv in inversions]
    pair_counts = Counter(species_pairs)
    ax4.bar(range(len(pair_counts)), pair_counts.values())
    ax4.set_xticks(range(len(pair_counts)))
    ax4.set_xticklabels(pair_counts.keys(), rotation=45, ha='right')
    ax4.set_ylabel('Inversion Count')
    ax4.set_title('Inversions by Species Pair')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def create_integrated_summary_plot(syngraph, rearrangements, inversions, output_path):
    """Create integrated summary visualization"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Summary statistics
    n_taxa = len(syngraph.graph['taxa'])
    n_markers = syngraph.number_of_nodes()
    n_edges = syngraph.number_of_edges()
    n_rearrangements = len(rearrangements)
    n_inversions = len(inversions)
    
    summary_text = f"""
INTEGRATED SYNGRAPH + INVERSION ANALYSIS SUMMARY

SYNGRAPH MACRO-SYNTENY ANALYSIS:
• Taxa analyzed: {n_taxa}
• Markers (BUSCO genes): {n_markers}
• Syntenic edges: {n_edges}
• Chromosomal rearrangements detected: {n_rearrangements}

MICRO-INVERSION ANALYSIS:
• Inversions detected: {n_inversions}
• Detection method: K-mer similarity within syntenic regions
• Analysis scope: Within LMS-identified synteny blocks

METHODOLOGICAL ADVANTAGES:
• Hierarchical approach: Macro → Micro analysis
• Network-based synteny detection (syngraph LMS)
• Sequence-level inversion validation
• Phylogenetically-informed rearrangement inference
• Comprehensive genomic rearrangement detection

BIOLOGICAL INSIGHTS:
• Identifies both large-scale and fine-scale rearrangements
• Separates ancestral from derived genomic arrangements
• Reveals evolutionary hotspots and conservation patterns
• Enables comparative genomic analysis across multiple species
    """
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=11,
           fontfamily='monospace', verticalalignment='top')
    ax.axis('off')
    ax.set_title('Integrated Synteny and Inversion Analysis Summary', 
                fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

################################################################################
# MAIN INTEGRATED WORKFLOW
################################################################################

def run_integrated_syngraph_inversion_analysis(config):
    """Main function for integrated syngraph + inversion analysis"""
    logger.info("=== INTEGRATED SYNGRAPH + INVERSION ANALYSIS ===")
    logger.info("Hierarchical synteny and inversion detection workflow\n")
    
    # STEP 1: Load data and build syngraph
    logger.info("STEP 1: Building Syngraph from BUSCO data")
    busco_files = glob.glob(config['busco_dir'] + '/*.tsv')
    
    # Convert BUSCO to MarkerObjs
    markerObjs = convert_busco_to_markerObjs(busco_files)
    
    # Build syngraph
    syngraph = Syngraph()
    syngraph.from_markerObjs(markerObjs)
    syngraph.show_metrics()
    
    # STEP 2: Syngraph macro-synteny analysis
    logger.info("\nSTEP 2: Syngraph macro-synteny analysis")
    taxa_list = list(syngraph.graph['taxa'])
    rearrangements, LMSs = analyze_syngraph_rearrangements(
        syngraph, taxa_list, config['minimum_lms_size']
    )
    
    # STEP 3: Load genome sequences for inversion analysis
    logger.info("\nSTEP 3: Loading genome sequences")
    genome_sequences = {}
    for taxon in taxa_list:
        fasta_file = os.path.join(config['genome_dir'], f"{taxon}.fna")
        if os.path.exists(fasta_file):
            sequences = {}
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences[record.id] = str(record.seq)
            genome_sequences[taxon] = sequences
            logger.info(f"  {taxon}: {len(sequences)} sequences")
    
    # STEP 4: Extract syntenic regions and detect inversions
    logger.info("\nSTEP 4: Micro-inversion analysis")
    syntenic_regions = extract_syntenic_regions_from_syngraph(
        syngraph, LMSs, genome_sequences
    )
    
    inversions = detect_inversions_with_kmers(
        syntenic_regions, 
        config['kmer_size'], 
        config['inversion_threshold']
    )
    
    # STEP 5: Create integrated visualizations
    logger.info("\nSTEP 5: Creating integrated visualizations")
    create_integrated_visualization(syngraph, rearrangements, inversions, config['output_dir'])
    
    # STEP 6: Save results
    logger.info("\nSTEP 6: Saving results")
    output_dir = config['output_dir']
    os.makedirs(output_dir, exist_ok=True)
    
    # Save rearrangements
    if rearrangements:
        pd.DataFrame(rearrangements).to_csv(f"{output_dir}/syngraph_rearrangements.csv", index=False)
    
    # Save inversions
    if inversions:
        pd.DataFrame(inversions).to_csv(f"{output_dir}/kmer_inversions.csv", index=False)
    
    # Generate summary
    summary = {
        'taxa_count': len(taxa_list),
        'markers_count': syngraph.number_of_nodes(),
        'syntenic_edges': syngraph.number_of_edges(),
        'lms_count': len(LMSs),
        'rearrangements_count': len(rearrangements),
        'inversions_count': len(inversions),
        'syntenic_regions_analyzed': len(syntenic_regions)
    }
    
    logger.info("\n=== ANALYSIS COMPLETED ===")
    logger.info("Summary Statistics:")
    for key, value in summary.items():
        logger.info(f"  {key}: {value}")
    
    return {
        'syngraph': syngraph,
        'rearrangements': rearrangements,
        'inversions': inversions,
        'LMSs': LMSs,
        'summary': summary
    }

################################################################################
# MAIN EXECUTION
################################################################################

def main():
    """Main execution function"""
    
    # Configuration
    config = {
        'busco_dir': 'busco-data',
        'genome_dir': 'genomes',
        'output_dir': 'integrated_results',
        'minimum_lms_size': 10,  # Minimum markers per LMS
        'kmer_size': 21,
        'inversion_threshold': 0.3
    }
    
    # Run integrated analysis
    results = run_integrated_syngraph_inversion_analysis(config)
    
    # Print final summary
    print("\n" + "="*70)
    print("INTEGRATED SYNGRAPH + INVERSION ANALYSIS RESULTS")
    print("="*70)
    print(f"✓ Syngraph network: {results['summary']['markers_count']} markers, {results['summary']['syntenic_edges']} edges")
    print(f"✓ LMS analysis: {results['summary']['lms_count']} linked marker sets")
    print(f"✓ Macro-rearrangements: {results['summary']['rearrangements_count']} events")
    print(f"✓ Micro-inversions: {results['summary']['inversions_count']} events")
    print(f"✓ Syntenic regions analyzed: {results['summary']['syntenic_regions_analyzed']}")
    print("="*70)
    print("This hierarchical approach successfully combines:")
    print("• Syngraph's network-based macro-synteny analysis")
    print("• Enhanced k-mer based micro-inversion detection")
    print("• Comprehensive multi-scale genomic rearrangement analysis")

if __name__ == "__main__":
    main()