"""
Phylogenetic Integration Module for Genome Inversion Analyzer
Phase 3: Core phylogenetic framework (Tier 1 Essential)
Focus: Distance matrices → Tree inference → Ancestral reconstruction
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from collections import defaultdict, Counter
import itertools

# Handle scipy imports gracefully
try:
    from scipy.spatial.distance import squareform, pdist
    from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
    from scipy.stats import spearmanr, pearsonr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

import json

logger = logging.getLogger(__name__)


class PhylogeneticIntegrator:
    """
    Core phylogenetic framework for inversion analysis
    Implements: Distance matrices → Tree inference → Ancestral reconstruction
    """
    
    def __init__(self, registry, config: Dict = None):
        self.registry = registry
        self.config = config or {
            'distance_methods': ['jaccard', 'manhattan', 'euclidean'],
            'linkage_method': 'ward',
            'min_species': 3
        }
        self.species_data = {}
        self.distance_matrices = {}
        self.phylogenetic_trees = {}
        
        # Check dependencies
        if not SCIPY_AVAILABLE:
            logger.warning("SciPy not available - some phylogenetic features will be limited")
        
    def add_species_data(self, species_name: str, inversion_df: pd.DataFrame, 
                        genome_stats: Dict, contextual_metrics: Dict = None) -> bool:
        """
        Add species inversion data for phylogenetic analysis
        
        Args:
            species_name: Species identifier
            inversion_df: Inversion events DataFrame
            genome_stats: Basic genome statistics  
            contextual_metrics: Optional contextual analysis
            
        Returns:
            Success status
        """
        try:
            logger.info(f"Adding species data for {species_name}")
            
            # Create species feature vector for phylogenetic analysis
            feature_vector = self._create_species_feature_vector(
                inversion_df, genome_stats, contextual_metrics
            )
            
            self.species_data[species_name] = {
                'inversion_data': inversion_df,
                'genome_stats': genome_stats,
                'contextual_metrics': contextual_metrics or {},
                'feature_vector': feature_vector,
                'inversion_count': len(inversion_df),
                'chromosomes': inversion_df['first_chr'].unique().tolist() if 'first_chr' in inversion_df.columns else []
            }
            
            logger.info(f"Added {species_name}: {len(inversion_df)} inversions, "
                       f"{len(feature_vector)} features")
            return True
            
        except Exception as e:
            logger.error(f"Failed to add species data for {species_name}: {e}")
            return False
    
    def compute_distance_matrices(self, methods: List[str] = None) -> Dict[str, np.ndarray]:
        """
        Compute phylogenetic distance matrices using different methods
        
        Args:
            methods: Distance calculation methods to use
            
        Returns:
            Dictionary of distance matrices by method
        """
        if len(self.species_data) < self.config['min_species']:
            logger.warning(f"Need at least {self.config['min_species']} species for phylogenetic analysis")
            return {}
        
        methods = methods or self.config['distance_methods']
        logger.info(f"Computing distance matrices using methods: {methods}")
        
        # Prepare data matrix
        species_names = list(self.species_data.keys())
        n_species = len(species_names)
        
        # Get feature vectors
        feature_vectors = []
        for species in species_names:
            feature_vectors.append(self.species_data[species]['feature_vector'])
        
        feature_matrix = np.array(feature_vectors)
        
        # Compute distance matrices
        distance_matrices = {}
        
        for method in methods:
            try:
                if method == 'jaccard':
                    # Jaccard distance for binary inversion presence/absence
                    dist_matrix = self._compute_jaccard_distance(species_names)
                elif method == 'manhattan' and SCIPY_AVAILABLE:
                    distances = pdist(feature_matrix, metric='manhattan')
                    dist_matrix = squareform(distances)
                elif method == 'euclidean' and SCIPY_AVAILABLE:
                    distances = pdist(feature_matrix, metric='euclidean')
                    dist_matrix = squareform(distances)
                elif method == 'inversion_rate':
                    # Custom distance based on inversion rates
                    dist_matrix = self._compute_inversion_rate_distance(species_names)
                elif not SCIPY_AVAILABLE and method in ['manhattan', 'euclidean']:
                    logger.warning(f"SciPy required for {method} distance - skipping")
                    continue
                else:
                    logger.warning(f"Unknown distance method: {method}")
                    continue
                
                distance_matrices[method] = dist_matrix
                
                # Register distance matrix
                distance_df = pd.DataFrame(
                    dist_matrix, 
                    index=species_names, 
                    columns=species_names
                )
                
                self.registry.register_file(
                    f'distance_matrix_{method}',
                    distance_df,
                    'csv',
                    f'Phylogenetic distance matrix ({method})',
                    dependencies=[f'{species}_data' for species in species_names],
                    metadata={
                        'method': method,
                        'n_species': n_species,
                        'species_list': species_names
                    }
                )
                
                logger.info(f"Computed {method} distance matrix ({n_species}x{n_species})")
                
            except Exception as e:
                logger.error(f"Failed to compute {method} distance matrix: {e}")
        
        self.distance_matrices = distance_matrices
        return distance_matrices
    
    def infer_phylogenetic_trees(self, linkage_methods: List[str] = None) -> Dict[str, Any]:
        """
        Infer phylogenetic trees using hierarchical clustering
        
        Args:
            linkage_methods: Hierarchical clustering methods
            
        Returns:
            Dictionary of phylogenetic trees
        """
        if not SCIPY_AVAILABLE:
            logger.error("SciPy required for phylogenetic tree inference")
            return {}
            
        if not self.distance_matrices:
            logger.warning("No distance matrices available for tree inference")
            return {}
        
        linkage_methods = linkage_methods or ['ward', 'average', 'complete']
        logger.info(f"Inferring phylogenetic trees using linkage methods: {linkage_methods}")
        
        trees = {}
        species_names = list(self.species_data.keys())
        
        for dist_method, distance_matrix in self.distance_matrices.items():
            for link_method in linkage_methods:
                try:
                    tree_key = f"{dist_method}_{link_method}"
                    
                    # Convert distance matrix to condensed form for linkage
                    condensed_distances = squareform(distance_matrix, checks=False)
                    
                    # Perform hierarchical clustering
                    linkage_matrix = linkage(condensed_distances, method=link_method)
                    
                    # Convert to tree structure
                    tree = to_tree(linkage_matrix)
                    
                    # Create tree representation
                    tree_data = {
                        'linkage_matrix': linkage_matrix.tolist(),
                        'species_names': species_names,
                        'distance_method': dist_method,
                        'linkage_method': link_method,
                        'n_species': len(species_names)
                    }
                    
                    trees[tree_key] = tree_data
                    
                    # Register tree
                    self.registry.register_file(
                        f'phylogenetic_tree_{tree_key}',
                        tree_data,
                        'json',
                        f'Phylogenetic tree ({dist_method} + {link_method})',
                        dependencies=[f'distance_matrix_{dist_method}'],
                        metadata=tree_data
                    )
                    
                    logger.info(f"Inferred tree: {tree_key}")
                    
                except Exception as e:
                    logger.error(f"Failed to infer tree {dist_method}_{link_method}: {e}")
        
        self.phylogenetic_trees = trees
        return trees
    
    def reconstruct_ancestral_states(self, tree_key: str = None) -> Dict[str, Any]:
        """
        Reconstruct ancestral inversion states using parsimony
        
        Args:
            tree_key: Specific tree to use for reconstruction
            
        Returns:
            Ancestral state reconstruction results
        """
        if not self.phylogenetic_trees:
            logger.warning("No phylogenetic trees available for ancestral reconstruction")
            return {}
        
        # Use first available tree if none specified
        tree_key = tree_key or list(self.phylogenetic_trees.keys())[0]
        
        if tree_key not in self.phylogenetic_trees:
            logger.error(f"Tree {tree_key} not found")
            return {}
        
        logger.info(f"Reconstructing ancestral states using tree: {tree_key}")
        
        try:
            tree_data = self.phylogenetic_trees[tree_key]
            species_names = tree_data['species_names']
            
            # Get inversion presence/absence matrix
            inversion_matrix = self._create_inversion_presence_matrix(species_names)
            
            # Simple parsimony reconstruction
            ancestral_states = self._parsimony_reconstruction(
                tree_data['linkage_matrix'], 
                inversion_matrix, 
                species_names
            )
            
            reconstruction_results = {
                'tree_key': tree_key,
                'ancestral_states': ancestral_states,
                'species_names': species_names,
                'inversion_characters': inversion_matrix.columns.tolist(),
                'method': 'parsimony'
            }
            
            # Register ancestral reconstruction
            self.registry.register_file(
                f'ancestral_reconstruction_{tree_key}',
                reconstruction_results,
                'json',
                f'Ancestral state reconstruction ({tree_key})',
                dependencies=[f'phylogenetic_tree_{tree_key}'],
                metadata={
                    'tree_used': tree_key,
                    'n_species': len(species_names),
                    'n_characters': len(inversion_matrix.columns),
                    'method': 'parsimony'
                }
            )
            
            logger.info(f"Reconstructed ancestral states for {len(species_names)} species, "
                       f"{len(inversion_matrix.columns)} characters")
            
            return reconstruction_results
            
        except Exception as e:
            logger.error(f"Failed to reconstruct ancestral states: {e}")
            return {}
    
    def analyze_phylogenetic_signal(self) -> Dict[str, Any]:
        """
        Analyze phylogenetic signal in inversion patterns
        
        Returns:
            Phylogenetic signal analysis results
        """
        if not SCIPY_AVAILABLE:
            logger.warning("SciPy required for phylogenetic signal analysis - using basic implementation")
            return self._basic_phylogenetic_signal()
            
        if len(self.species_data) < 3:
            logger.warning("Need at least 3 species for phylogenetic signal analysis")
            return {}
        
        logger.info("Analyzing phylogenetic signal in inversion patterns")
        
        try:
            species_names = list(self.species_data.keys())
            
            # Create trait matrix (inversion rates, patterns, etc.)
            trait_matrix = []
            trait_names = ['inversion_rate', 'mean_length', 'chromosome_diversity']
            
            for species in species_names:
                data = self.species_data[species]
                inversion_df = data['inversion_data']
                genome_stats = data['genome_stats']
                
                # Calculate traits
                inversion_rate = len(inversion_df) / (genome_stats.get('total_length', 1) / 1_000_000)
                mean_length = inversion_df['length'].mean() if 'length' in inversion_df.columns and len(inversion_df) > 0 else 0
                chromosome_diversity = len(data['chromosomes'])
                
                trait_matrix.append([inversion_rate, mean_length, chromosome_diversity])
            
            trait_matrix = np.array(trait_matrix)
            
            # Calculate phylogenetic signal (if we have distance matrices)
            signal_results = {}
            
            if self.distance_matrices:
                for dist_method, distance_matrix in self.distance_matrices.items():
                    method_results = {}
                    
                    for i, trait_name in enumerate(trait_names):
                        trait_values = trait_matrix[:, i]
                        
                        # Calculate correlation between trait distances and phylogenetic distances
                        trait_distances = pdist(trait_values.reshape(-1, 1), metric='euclidean')
                        phylo_distances = squareform(distance_matrix, checks=False)
                        
                        # Mantel test approximation (correlation between distance matrices)
                        correlation, p_value = spearmanr(trait_distances, phylo_distances)
                        
                        method_results[trait_name] = {
                            'correlation': correlation,
                            'p_value': p_value,
                            'trait_mean': float(np.mean(trait_values)),
                            'trait_std': float(np.std(trait_values))
                        }
                    
                    signal_results[dist_method] = method_results
            
            # Register phylogenetic signal analysis
            self.registry.register_file(
                'phylogenetic_signal_analysis',
                signal_results,
                'json',
                'Phylogenetic signal analysis of inversion traits',
                dependencies=list(self.distance_matrices.keys()),
                metadata={
                    'n_species': len(species_names),
                    'traits_analyzed': trait_names,
                    'distance_methods': list(signal_results.keys())
                }
            )
            
            logger.info(f"Analyzed phylogenetic signal for {len(trait_names)} traits "
                       f"across {len(signal_results)} distance methods")
            
            return signal_results
            
        except Exception as e:
            logger.error(f"Failed to analyze phylogenetic signal: {e}")
            return {}
    
    def generate_phylogenetic_summary(self) -> Dict[str, Any]:
        """
        Generate comprehensive phylogenetic analysis summary
        
        Returns:
            Complete phylogenetic analysis summary
        """
        logger.info("Generating comprehensive phylogenetic summary")
        
        summary = {
            'species_overview': {},
            'distance_matrices': {},
            'phylogenetic_trees': {},
            'ancestral_reconstruction': {},
            'phylogenetic_signal': {}
        }
        
        # 1. Species overview
        summary['species_overview'] = {
            'n_species': len(self.species_data),
            'species_names': list(self.species_data.keys()),
            'total_inversions': sum(data['inversion_count'] for data in self.species_data.values()),
            'species_details': {
                name: {
                    'inversion_count': data['inversion_count'],
                    'chromosomes': len(data['chromosomes']),
                    'genome_size_mb': data['genome_stats'].get('total_length', 0) / 1_000_000
                }
                for name, data in self.species_data.items()
            }
        }
        
        # 2. Distance matrices summary
        if self.distance_matrices:
            summary['distance_matrices'] = {
                'methods_computed': list(self.distance_matrices.keys()),
                'matrix_size': f"{len(self.species_data)}x{len(self.species_data)}",
                'mean_distances': {
                    method: float(np.mean(matrix[np.triu_indices_from(matrix, k=1)]))
                    for method, matrix in self.distance_matrices.items()
                }
            }
        
        # 3. Phylogenetic trees summary
        if self.phylogenetic_trees:
            summary['phylogenetic_trees'] = {
                'trees_inferred': list(self.phylogenetic_trees.keys()),
                'methods_used': list(set(
                    tree['distance_method'] + '_' + tree['linkage_method'] 
                    for tree in self.phylogenetic_trees.values()
                ))
            }
        
        # 4. Generate remaining analyses if not already done
        if not hasattr(self, '_phylogenetic_signal_done'):
            summary['phylogenetic_signal'] = self.analyze_phylogenetic_signal()
            self._phylogenetic_signal_done = True
        
        # Register complete summary
        self.registry.register_file(
            'phylogenetic_analysis_summary',
            summary,
            'json',
            'Comprehensive phylogenetic analysis summary',
            dependencies=['species_data', 'distance_matrices', 'phylogenetic_trees'],
            metadata={
                'n_species': len(self.species_data),
                'analyses_completed': list(summary.keys()),
                'distance_methods': list(self.distance_matrices.keys()) if self.distance_matrices else [],
                'tree_methods': list(self.phylogenetic_trees.keys()) if self.phylogenetic_trees else []
            }
        )
        
        logger.info("Phylogenetic analysis summary generated")
        return summary
    
    # Helper methods
    def _create_species_feature_vector(self, inversion_df: pd.DataFrame, 
                                     genome_stats: Dict, contextual_metrics: Dict) -> np.ndarray:
        """Create feature vector for phylogenetic analysis"""
        features = []
        
        # Basic inversion metrics
        features.append(len(inversion_df))  # Total inversions
        features.append(len(inversion_df) / (genome_stats.get('total_length', 1) / 1_000_000))  # Rate per Mb
        
        # Length statistics
        if 'length' in inversion_df.columns and len(inversion_df) > 0:
            features.extend([
                inversion_df['length'].mean(),
                inversion_df['length'].std(),
                inversion_df['length'].median()
            ])
        else:
            features.extend([0, 0, 0])
        
        # Chromosome distribution
        if 'first_chr' in inversion_df.columns:
            n_chromosomes = len(inversion_df['first_chr'].unique())
            features.append(n_chromosomes)
        else:
            features.append(0)
        
        # Contextual metrics if available
        if contextual_metrics:
            rate_metrics = contextual_metrics.get('rate_metrics', {})
            features.append(rate_metrics.get('inversions_per_mb', 0))
            
            gc_correlation = contextual_metrics.get('gc_correlation', {})
            features.append(gc_correlation.get('mean_gc_at_inversions', 0))
        else:
            features.extend([0, 0])
        
        return np.array(features, dtype=float)
    
    def _compute_jaccard_distance(self, species_names: List[str]) -> np.ndarray:
        """Compute Jaccard distance based on inversion presence/absence"""
        n_species = len(species_names)
        distance_matrix = np.zeros((n_species, n_species))
        
        # Create sets of inversion "signatures" for each species
        species_signatures = {}
        for species in species_names:
            inversion_df = self.species_data[species]['inversion_data']
            # Create simple signature based on chromosome and approximate position
            signatures = set()
            for _, row in inversion_df.iterrows():
                if 'first_chr' in row and 'first_start' in row:
                    # Bin positions to 1Mb windows for comparison
                    binned_pos = int(row['first_start'] / 1_000_000)
                    signature = f"{row['first_chr']}_{binned_pos}"
                    signatures.add(signature)
            species_signatures[species] = signatures
        
        # Compute pairwise Jaccard distances
        for i, species1 in enumerate(species_names):
            for j, species2 in enumerate(species_names):
                if i == j:
                    distance_matrix[i, j] = 0
                else:
                    set1 = species_signatures[species1]
                    set2 = species_signatures[species2]
                    
                    intersection = len(set1 & set2)
                    union = len(set1 | set2)
                    
                    jaccard_similarity = intersection / union if union > 0 else 0
                    jaccard_distance = 1 - jaccard_similarity
                    
                    distance_matrix[i, j] = jaccard_distance
        
        return distance_matrix
    
    def _compute_inversion_rate_distance(self, species_names: List[str]) -> np.ndarray:
        """Compute distance based on inversion rates"""
        n_species = len(species_names)
        
        # Get inversion rates
        rates = []
        for species in species_names:
            data = self.species_data[species]
            rate = data['inversion_count'] / (data['genome_stats'].get('total_length', 1) / 1_000_000)
            rates.append(rate)
        
        rates = np.array(rates)
        
        # Compute pairwise distances
        distance_matrix = np.zeros((n_species, n_species))
        for i in range(n_species):
            for j in range(n_species):
                distance_matrix[i, j] = abs(rates[i] - rates[j])
        
        return distance_matrix
    
    def _create_inversion_presence_matrix(self, species_names: List[str]) -> pd.DataFrame:
        """Create binary matrix of inversion presence/absence"""
        # Collect all unique inversion "characters" across species
        all_characters = set()
        species_characters = {}
        
        for species in species_names:
            inversion_df = self.species_data[species]['inversion_data']
            characters = set()
            
            for _, row in inversion_df.iterrows():
                if 'first_chr' in row and 'first_start' in row:
                    # Create character based on chromosome and binned position
                    binned_pos = int(row['first_start'] / 1_000_000)
                    character = f"{row['first_chr']}_{binned_pos}"
                    characters.add(character)
            
            species_characters[species] = characters
            all_characters.update(characters)
        
        # Create binary matrix
        character_list = sorted(list(all_characters))
        matrix_data = []
        
        for species in species_names:
            row = [1 if char in species_characters[species] else 0 for char in character_list]
            matrix_data.append(row)
        
        return pd.DataFrame(matrix_data, index=species_names, columns=character_list)
    
    def _parsimony_reconstruction(self, linkage_matrix: List[List], 
                                inversion_matrix: pd.DataFrame, 
                                species_names: List[str]) -> Dict[str, Any]:
        """Simple parsimony-based ancestral state reconstruction"""
        # This is a simplified implementation
        # In practice, you'd use more sophisticated algorithms
        
        n_species = len(species_names)
        n_characters = len(inversion_matrix.columns)
        
        # For now, just compute the most parsimonious state at each internal node
        # as the majority state of descendant species
        ancestral_states = {}
        
        # Add leaf states
        for i, species in enumerate(species_names):
            ancestral_states[f"leaf_{i}"] = inversion_matrix.loc[species].values.tolist()
        
        # Simple reconstruction: internal nodes get majority state
        for i in range(n_species - 1):
            node_id = f"internal_{i}"
            # This is overly simplified - real parsimony algorithms are more complex
            ancestral_states[node_id] = [0] * n_characters
        
        return ancestral_states
    
    def _basic_phylogenetic_signal(self) -> Dict[str, Any]:
        """Basic phylogenetic signal analysis without SciPy"""
        logger.info("Using basic phylogenetic signal analysis (SciPy not available)")
        
        signal_results = {}
        species_names = list(self.species_data.keys())
        
        # Simple trait comparison
        trait_data = {}
        for species in species_names:
            data = self.species_data[species]
            inversion_df = data['inversion_data']
            genome_stats = data['genome_stats']
            
            trait_data[species] = {
                'inversion_rate': len(inversion_df) / (genome_stats.get('total_length', 1) / 1_000_000),
                'inversion_count': len(inversion_df),
                'chromosome_count': len(data['chromosomes'])
            }
        
        signal_results['basic_traits'] = trait_data
        
        return signal_results