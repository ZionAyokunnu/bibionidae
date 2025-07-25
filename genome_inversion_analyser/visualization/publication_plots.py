"""
Publication-quality plots for genome inversion analysis
Integrates with external tools like synteny_plotter and provides tree annotation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import subprocess
import tempfile
import shutil
import logging
import json
from typing import Dict, List, Optional, Tuple, Any

logger = logging.getLogger(__name__)

try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    logger.warning("ete3 not available - tree manipulation will be limited")


class PublicationPlotGenerator:
    """
    Creates publication-quality plots using external tools and annotations
    """
    
    def __init__(self, registry, config: Dict):
        self.registry = registry
        self.config = config
        self.pub_config = config.get('publication_config', {})
        
    def create_publication_suite(self, all_results: Dict, species_stats: Dict, output_dir: Path):
        """Create complete publication visualization suite"""
        
        if not self.pub_config.get('publication_suite', {}).get('enabled', False):
            logger.info("Publication suite disabled in config")
            return
            
        logger.info("🎨 Creating publication-quality visualization suite...")
        
        pub_dir = output_dir / 'publication_plots'
        pub_dir.mkdir(exist_ok=True)
        
        results = {}
        
        # 1. Curved synteny plots
        if self.pub_config.get('synteny_visualization', {}).get('enabled', False):
            logger.info("  📊 Starting synteny visualization...")
            synteny_results = self._create_synteny_plots(all_results, pub_dir)
            results['synteny_plots'] = synteny_results
            
        # 2. Annotated phylogenetic tree
        if self.pub_config.get('tree_annotation', {}).get('enabled', False):
            logger.info("  🌳 Starting tree annotation...")
            tree_results = self._create_annotated_phylogeny(all_results, species_stats, pub_dir)
            results['annotated_trees'] = tree_results
            
        # 3. Register all outputs
        for plot_type, plot_data in results.items():
            if isinstance(plot_data, dict):
                for plot_name, plot_path in plot_data.items():
                    self.registry.register_file(
                        f'publication_{plot_type}_{plot_name}',
                        plot_path,
                        'visualization',
                        f'Publication-quality {plot_type}: {plot_name}'
                    )
        
        logger.info(f"✅ Publication suite completed: {len(results)} plot types created")
        return results
    
    def _create_synteny_plots(self, all_results: Dict, output_dir: Path) -> Dict:
        """Create publication-quality synteny plots using synteny_plotter"""
        
        logger.info("  📊 Creating curved synteny plots...")
        
        synteny_dir = output_dir / 'synteny_curves'
        synteny_dir.mkdir(exist_ok=True)
        
        synteny_plotter_path = self.pub_config.get('external_tools', {}).get('synteny_plotter')
        
        if not synteny_plotter_path or not Path(synteny_plotter_path).exists():
            logger.warning(f"Synteny plotter not found: {synteny_plotter_path}")
            logger.info("  Creating fallback matplotlib synteny plots...")
            return self._create_fallback_synteny_plots(all_results, synteny_dir)
        
        results = {}
        
        for pair_name, pair_data in all_results.items():
            if 'full_results' not in pair_data:
                continue
                
            try:
                species1, species2 = pair_data['species_pair']
                logger.info(f"    • Creating synteny plot: {pair_name}")
                
                # Convert data to synteny_plotter format
                temp_files = self._prepare_synteny_plotter_input(
                    pair_data['ortholog_df'], 
                    pair_data['inversion_df'],
                    species1, species2,
                    synteny_dir
                )
                
                # Call synteny_plotter
                output_plot = synteny_dir / f'curved_synteny_{pair_name}.png'
                
                success = self._run_synteny_plotter(temp_files, output_plot, species1, species2)
                
                if success:
                    results[pair_name] = output_plot
                    logger.info(f"      ✅ Created: {output_plot}")
                else:
                    logger.warning(f"      ❌ Failed: {pair_name}, creating fallback")
                    # Create fallback plot
                    fallback_plot = self._create_single_fallback_plot(
                        pair_data['ortholog_df'], pair_data['inversion_df'], 
                        species1, species2, synteny_dir
                    )
                    if fallback_plot:
                        results[pair_name] = fallback_plot
                
            except Exception as e:
                logger.error(f"      ❌ Error creating synteny plot for {pair_name}: {e}")
        
        logger.info(f"  ✅ Created {len(results)} synteny plots")
        return results
    
    def _prepare_synteny_plotter_input(self, ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame,
                                     species1: str, species2: str, output_dir: Path) -> Dict:
        """Convert ortholog data to synteny_plotter BED format"""
        
        temp_dir = output_dir / f'temp_{species1}_vs_{species2}'
        temp_dir.mkdir(exist_ok=True)
        
        # Create BED files for each species
        bed1_data = []
        bed2_data = []
        links_data = []
        
        for _, row in ortholog_df.iterrows():
            # BED format: chr start end name score strand
            bed1_data.append([
                row['first_chromosome'],
                int(row['first_start']),
                int(row['first_end']),
                row['busco_id'],
                int(row['similarity'] * 1000),  # Convert to integer score
                row.get('first_strand', '+')
            ])
            
            bed2_data.append([
                row['second_chromosome'],
                int(row['second_start']),
                int(row['second_end']),
                row['busco_id'],
                int(row['similarity'] * 1000),
                row.get('second_strand', '+')
            ])
            
            # Links format: chr1 start1 end1 chr2 start2 end2 similarity
            links_data.append([
                row['first_chromosome'],
                int(row['first_start']),
                int(row['first_end']),
                row['second_chromosome'],
                int(row['second_start']),
                int(row['second_end']),
                row['similarity']
            ])
        
        # Write files
        bed1_file = temp_dir / f'{species1}.bed'
        bed2_file = temp_dir / f'{species2}.bed'
        links_file = temp_dir / f'{species1}_vs_{species2}_links.tsv'
        
        # Write BED files
        with open(bed1_file, 'w') as f:
            for row in bed1_data:
                f.write('\t'.join(map(str, row)) + '\n')
                
        with open(bed2_file, 'w') as f:
            for row in bed2_data:
                f.write('\t'.join(map(str, row)) + '\n')
        
        # Write links file
        with open(links_file, 'w') as f:
            f.write('chr1\tstart1\tend1\tchr2\tstart2\tend2\tsimilarity\n')
            for row in links_data:
                f.write('\t'.join(map(str, row)) + '\n')
        
        return {
            'bed1': bed1_file,
            'bed2': bed2_file,
            'links': links_file,
            'temp_dir': temp_dir
        }
    
    def _run_synteny_plotter(self, temp_files: Dict, output_plot: Path, 
                           species1: str, species2: str) -> bool:
        """Run synteny_plotter with the prepared input files"""
        
        synteny_config = self.pub_config.get('synteny_visualization', {})
        options = synteny_config.get('options', {})
        
        cmd = [
            'python', self.pub_config['external_tools']['synteny_plotter'],
            '--bed1', str(temp_files['bed1']),
            '--bed2', str(temp_files['bed2']),
            '--links', str(temp_files['links']),
            '--output', str(output_plot),
            '--species1', species1,
            '--species2', species2
        ]
        
        # Add optional parameters
        if options.get('curved', True):
            cmd.append('--curved')
        if options.get('color_by'):
            cmd.extend(['--color_by', options['color_by']])
        if options.get('figure_width'):
            cmd.extend(['--figure_width', str(options['figure_width'])])
        if options.get('figure_height'):
            cmd.extend(['--figure_height', str(options['figure_height'])])
        if options.get('dpi'):
            cmd.extend(['--dpi', str(options['dpi'])])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Clean up temp files
            if temp_files['temp_dir'].exists():
                shutil.rmtree(temp_files['temp_dir'])
            
            return output_plot.exists()
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Synteny plotter failed: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Error running synteny plotter: {e}")
            return False
    
    def _create_fallback_synteny_plots(self, all_results: Dict, output_dir: Path) -> Dict:
        """Create matplotlib-based synteny plots as fallback"""
        
        logger.info("    Creating matplotlib fallback synteny plots...")
        
        results = {}
        
        for pair_name, pair_data in all_results.items():
            if 'full_results' not in pair_data:
                continue
                
            try:
                species1, species2 = pair_data['species_pair']
                plot_file = self._create_single_fallback_plot(
                    pair_data['ortholog_df'], 
                    pair_data['inversion_df'],
                    species1, species2, 
                    output_dir
                )
                
                if plot_file:
                    results[pair_name] = plot_file
                    logger.info(f"    ✅ Fallback plot: {pair_name}")
                    
            except Exception as e:
                logger.error(f"    ❌ Fallback plot failed for {pair_name}: {e}")
        
        return results
    
    def _create_single_fallback_plot(self, ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame,
                                species1: str, species2: str, output_dir: Path) -> Optional[Path]:
        """Create a single matplotlib synteny plot with curved connections"""
        
        if ortholog_df.empty:
            return None
        
        try:
            fig, ax = plt.subplots(figsize=(14, 10))
            
            # Fix column name mapping
            ortholog_df_mapped = ortholog_df.rename(columns={
                'first_chr': 'first_chromosome',
                'second_chr': 'second_chromosome'
            })
            
            # Get chromosome information
            chr1_list = sorted(ortholog_df_mapped['first_chromosome'].unique())
            chr2_list = sorted(ortholog_df_mapped['second_chromosome'].unique())
            
            # Create chromosome positions
            chr1_positions = {chr_name: i for i, chr_name in enumerate(chr1_list)}
            chr2_positions = {chr_name: i for i, chr_name in enumerate(chr2_list)}
            
            y1_base = 0.7  # Species 1 y-position
            y2_base = 0.3  # Species 2 y-position
            
            # Plot chromosome backbones
            for i, chr_name in enumerate(chr1_list):
                ax.plot([0, 1], [y1_base + i*0.05, y1_base + i*0.05], 'k-', linewidth=3, alpha=0.7)
                ax.text(-0.05, y1_base + i*0.05, chr_name, ha='right', va='center', fontsize=8)
            
            for i, chr_name in enumerate(chr2_list):
                ax.plot([0, 1], [y2_base - i*0.05, y2_base - i*0.05], 'k-', linewidth=3, alpha=0.7)
                ax.text(-0.05, y2_base - i*0.05, chr_name, ha='right', va='center', fontsize=8)
            
            # Plot ortholog connections with curves
            for _, row in ortholog_df_mapped.iterrows():
                chr1_idx = chr1_positions[row['first_chromosome']]
                chr2_idx = chr2_positions[row['second_chromosome']]
                
                y1 = y1_base + chr1_idx * 0.05
                y2 = y2_base - chr2_idx * 0.05
                
                # Normalized positions along chromosomes (improved calculation)
                x1 = 0.1 + 0.8 * (row['first_start'] % 10000000) / 10000000  # Simple normalization
                x2 = 0.1 + 0.8 * (row['second_start'] % 10000000) / 10000000
                
                # Color by similarity
                color = plt.cm.viridis(row['similarity'])
                alpha = 0.6
                
                # Check if this is an inversion
                is_inversion = False
                if not inversion_df.empty:
                    for _, inv_row in inversion_df.iterrows():
                        if (row['busco_id'] in str(inv_row.get('genes', '')) or 
                            abs(row['first_start'] - inv_row.get('first_start', 0)) < 10000):
                            is_inversion = True
                            break
                
                if is_inversion:
                    color = 'red'
                    alpha = 0.8
                
                # Create curved connection
                x_curve = np.linspace(x1, x2, 50)
                y_curve = []
                
                for x in x_curve:
                    # Bezier-like curve
                    t = (x - x1) / (x2 - x1) if x2 != x1 else 0
                    y = y1 * (1 - t) + y2 * t + 0.1 * np.sin(np.pi * t)  # Add curve
                    y_curve.append(y)
                
                ax.plot(x_curve, y_curve, color=color, alpha=alpha, linewidth=1)
            
            # Styling
            ax.set_xlim(-0.2, 1.1)
            ax.set_ylim(min(y2_base - len(chr2_list)*0.05, y1_base - 0.1), 
                    max(y1_base + len(chr1_list)*0.05, y2_base + 0.1))
            
            ax.set_title(f'Synteny Plot: {species1} vs {species2}', fontsize=16, fontweight='bold')
            ax.text(0.5, y1_base + len(chr1_list)*0.05 + 0.05, species1, 
                ha='center', va='bottom', fontsize=14, fontweight='bold')
            ax.text(0.5, y2_base - len(chr2_list)*0.05 - 0.05, species2, 
                ha='center', va='top', fontsize=14, fontweight='bold')
            
            ax.axis('off')
            
            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label('Similarity', rotation=270, labelpad=15)
            
            plt.tight_layout()
            
            # Save plot
            plot_file = output_dir / f'matplotlib_synteny_{species1}_vs_{species2}.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            return plot_file
            
        except Exception as e:
            logger.error(f"Fallback plot creation failed: {e}")
            plt.close()
            return None
    
    def _create_annotated_phylogeny(self, all_results: Dict, species_stats: Dict, 
                                  output_dir: Path) -> Dict:
        """Create annotated phylogenetic tree from Diptera tree"""
        
        logger.info("  🌳 Creating annotated phylogenetic tree...")
        
        tree_dir = output_dir / 'annotated_trees'
        tree_dir.mkdir(exist_ok=True)
        
        tree_config = self.pub_config.get('tree_annotation', {})
        source_tree_path = tree_config.get('source_tree_path')
        
        if not source_tree_path or not Path(source_tree_path).exists():
            logger.error(f"Source tree not found: {source_tree_path}")
            return {}
        
        try:
            # Load and prune tree
            target_species = list(species_stats.keys())
            pruned_tree = self._prune_diptera_tree(source_tree_path, target_species)
            
            if not pruned_tree:
                logger.error("Failed to prune tree")
                return {}
            
            # Calculate inversion annotations
            inversion_annotations = self._calculate_inversion_annotations(all_results, species_stats)
            
            # Annotate tree with inversion data
            annotated_tree = self._annotate_tree_with_inversions(pruned_tree, inversion_annotations)
            
            # Create visualizations
            results = {}
            
            # Save annotated tree
            tree_file = tree_dir / 'annotated_bibionidae_tree.newick'
            if ETE3_AVAILABLE:
                annotated_tree.write(outfile=str(tree_file))
                results['newick'] = tree_file
            
            # Create tree plot
            plot_file = tree_dir / 'annotated_tree_plot.png'
            self._create_tree_plot(annotated_tree, plot_file, inversion_annotations)
            results['plot'] = plot_file
            
            # Create tree heatmap
            heatmap_file = tree_dir / 'tree_inversion_heatmap.png'
            self._create_tree_heatmap(inversion_annotations, heatmap_file)
            results['heatmap'] = heatmap_file
            
            logger.info(f"  ✅ Created annotated tree: {len(results)} outputs")
            return results
            
        except Exception as e:
            logger.error(f"Tree annotation failed: {e}")
            import traceback
            traceback.print_exc()
            return {}
    

        # def _create_single_fallback_plot(self, ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame,
        #                         species1: str, species2: str, output_dir: Path) -> Optional[Path]:
        # """Create a single matplotlib synteny plot with curved connections"""
        
        # if ortholog_df.empty:
        #     return None
        
        # try:
        #     fig, ax = plt.subplots(figsize=(14, 10))
            
        #     # Use correct column names directly (no renaming needed)
        #     # Get chromosome information using the actual column names
        #     chr1_list = sorted(ortholog_df['first_chr'].unique())
        #     chr2_list = sorted(ortholog_df['second_chr'].unique())
            
        #     # Create chromosome positions
        #     chr1_positions = {chr_name: i for i, chr_name in enumerate(chr1_list)}
        #     chr2_positions = {chr_name: i for i, chr_name in enumerate(chr2_list)}
            
        #     y1_base = 0.7  # Species 1 y-position
        #     y2_base = 0.3  # Species 2 y-position
            
        #     # Plot chromosome backbones
        #     for i, chr_name in enumerate(chr1_list):
        #         ax.plot([0, 1], [y1_base + i*0.05, y1_base + i*0.05], 'k-', linewidth=3, alpha=0.7)
        #         # Truncate long chromosome names
        #         display_name = str(chr_name)[:10] + "..." if len(str(chr_name)) > 10 else str(chr_name)
        #         ax.text(-0.05, y1_base + i*0.05, display_name, ha='right', va='center', fontsize=8)
            
        #     for i, chr_name in enumerate(chr2_list):
        #         ax.plot([0, 1], [y2_base - i*0.05, y2_base - i*0.05], 'k-', linewidth=3, alpha=0.7)
        #         # Truncate long chromosome names
        #         display_name = str(chr_name)[:10] + "..." if len(str(chr_name)) > 10 else str(chr_name)
        #         ax.text(-0.05, y2_base - i*0.05, display_name, ha='right', va='center', fontsize=8)
            
        #     # Plot ortholog connections with curves
        #     for _, row in ortholog_df.iterrows():
        #         try:
        #             chr1_idx = chr1_positions[row['first_chr']]
        #             chr2_idx = chr2_positions[row['second_chr']]
                    
        #             y1 = y1_base + chr1_idx * 0.05
        #             y2 = y2_base - chr2_idx * 0.05
                    
        #             # Normalized positions along chromosomes (simple approach)
        #             x1 = 0.1 + 0.8 * ((row['first_start'] % 10000000) / 10000000)
        #             x2 = 0.1 + 0.8 * ((row['second_start'] % 10000000) / 10000000)
                    
        #             # Color by similarity
        #             color = plt.cm.viridis(row['similarity'])
        #             alpha = 0.6
                    
        #             # Check if this is an inversion (simple check)
        #             is_inversion = False
        #             if not inversion_df.empty and 'genes' in inversion_df.columns:
        #                 for _, inv_row in inversion_df.iterrows():
        #                     if row['busco_id'] in str(inv_row.get('genes', '')):
        #                         is_inversion = True
        #                         break
                    
        #             if is_inversion:
        #                 color = 'red'
        #                 alpha = 0.8
                    
        #             # Create curved connection
        #             x_curve = np.linspace(x1, x2, 50)
        #             y_curve = []
                    
        #             for x in x_curve:
        #                 # Bezier-like curve
        #                 t = (x - x1) / (x2 - x1) if x2 != x1 else 0
        #                 y = y1 * (1 - t) + y2 * t + 0.1 * np.sin(np.pi * t)  # Add curve
        #                 y_curve.append(y)
                    
        #             ax.plot(x_curve, y_curve, color=color, alpha=alpha, linewidth=1)
                    
        #         except (KeyError, ValueError) as e:
        #             # Skip problematic rows
        #             continue
            
        #     # Styling
        #     ax.set_xlim(-0.2, 1.1)
        #     ax.set_ylim(min(y2_base - len(chr2_list)*0.05, y1_base - 0.1), 
        #             max(y1_base + len(chr1_list)*0.05, y2_base + 0.1))
            
        #     ax.set_title(f'Synteny Plot: {species1} vs {species2}', fontsize=16, fontweight='bold')
        #     ax.text(0.5, y1_base + len(chr1_list)*0.05 + 0.05, species1, 
        #         ha='center', va='bottom', fontsize=14, fontweight='bold')
        #     ax.text(0.5, y2_base - len(chr2_list)*0.05 - 0.05, species2, 
        #         ha='center', va='top', fontsize=14, fontweight='bold')
            
        #     ax.axis('off')
            
        #     # Add colorbar
        #     sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
        #     sm.set_array([])
        #     cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        #     cbar.set_label('Similarity', rotation=270, labelpad=15)
            
        #     plt.tight_layout()
            
        #     # Save plot
        #     plot_file = output_dir / f'matplotlib_synteny_{species1}_vs_{species2}.png'
        #     plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        #     plt.close()
            
        #     return plot_file
            
        # except Exception as e:
        #     logger.error(f"Fallback plot creation failed: {e}")
        #     import traceback
        #     traceback.print_exc()
        #     plt.close()
        #     return None
    
    def _prune_diptera_tree(self, tree_path: str, target_species: List[str]) -> Optional['Tree']:
        """Prune large Diptera tree to target species"""
        
        if not ETE3_AVAILABLE:
            logger.warning("ete3 not available - creating simple tree structure")
            return self._create_simple_tree(target_species)
        
        try:
            # Load tree
            tree = Tree(tree_path)
            
            # Find target species in tree
            found_species = {}
            for species in target_species:
                # Try different name formats
                possible_names = [
                    species,
                    species.replace('_', ' '),
                    ' '.join(species.split('_')[:2])  # Genus species format
                ]
                
                for name in possible_names:
                    nodes = tree.search_nodes(name=name)
                    if nodes:
                        found_species[species] = nodes[0]
                        break
            
            if len(found_species) < 2:
                logger.warning(f"Not enough species found in tree: {found_species.keys()}")
                logger.info("Creating simple tree structure")
                return self._create_simple_tree(target_species)
            
            logger.info(f"Found {len(found_species)} species in tree: {list(found_species.keys())}")
            
            # Get common ancestor and prune
            common_ancestor = tree.get_common_ancestor(list(found_species.values()))
            pruned_tree = common_ancestor.copy()
            
            # Keep only target species
            target_nodes = [found_species[sp] for sp in found_species.keys()]
            pruned_tree.prune(target_nodes, preserve_branch_length=True)
            
            return pruned_tree
            
        except Exception as e:
            logger.error(f"Tree pruning failed: {e}")
            return self._create_simple_tree(target_species)
    
    def _create_simple_tree(self, species_list: List[str]) -> Optional['Tree']:
        """Create simple tree structure if ete3 not available or pruning fails"""
        
        if not ETE3_AVAILABLE:
            logger.info("Creating matplotlib-based tree visualization")
            return None
        
        try:
            # Create simple balanced tree
            if len(species_list) == 4:
                tree_str = f"(({species_list[0]}:0.1,{species_list[1]}:0.1):0.1,({species_list[2]}:0.1,{species_list[3]}:0.1):0.1);"
            elif len(species_list) == 3:
                tree_str = f"(({species_list[0]}:0.1,{species_list[1]}:0.1):0.1,{species_list[2]}:0.1);"
            elif len(species_list) == 2:
                tree_str = f"({species_list[0]}:0.1,{species_list[1]}:0.1);"
            else:
                # For more species, create a simple star topology
                species_with_branch = [f"{sp}:0.1" for sp in species_list]
                tree_str = f"({','.join(species_with_branch)});"
            
            tree = Tree(tree_str)
            return tree
            
        except Exception as e:
            logger.error(f"Simple tree creation failed: {e}")
            return None
    
    def _calculate_inversion_annotations(self, all_results: Dict, species_stats: Dict) -> Dict:
        """Calculate inversion metrics for tree annotation"""
        
        inversion_data = {}
        
        for species_name, stats in species_stats.items():
            genome_size = stats['quality']['metrics'].get('total_length', 1000000)
            
            # Collect all inversions involving this species
            species_inversions = []
            for pair_name, pair_data in all_results.items():
                if 'full_results' in pair_data and species_name in pair_data['species_pair']:
                    inv_df = pair_data['inversion_df']
                    if not inv_df.empty:
                        species_inversions.extend(inv_df.to_dict('records'))
            
            # Calculate metrics
            total_inversions = len(species_inversions)
            rate_per_mb = total_inversions / (genome_size / 1_000_000)
            
            inversion_data[species_name] = {
                'total_inversions': total_inversions,
                'rate_per_mb': rate_per_mb,
                'genome_size': genome_size,
                'quality_score': stats['quality']['quality_score']
            }
        
        # Calculate normalized scores
        if inversion_data:
            max_rate = max(data['rate_per_mb'] for data in inversion_data.values())
            for species_name in inversion_data:
                if max_rate > 0:
                    inversion_data[species_name]['normalized_score'] = inversion_data[species_name]['rate_per_mb'] / max_rate
                else:
                    inversion_data[species_name]['normalized_score'] = 0
        
        return inversion_data
    
    def _annotate_tree_with_inversions(self, tree: 'Tree', inversion_data: Dict) -> 'Tree':
        """Annotate tree nodes with inversion data"""
        
        if not ETE3_AVAILABLE or not tree:
            return tree
        
        try:
            # Annotate leaf nodes
            for leaf in tree.get_leaves():
                species_name = leaf.name
                if species_name in inversion_data:
                    data = inversion_data[species_name]
                    leaf.add_features(
                        inversion_count=data['total_inversions'],
                        inversion_rate=data['rate_per_mb'],
                        normalized_score=data['normalized_score'],
                        genome_size=data['genome_size']
                    )
            
            # Annotate internal nodes (ancestral state reconstruction)
            for node in tree.traverse():
                if not node.is_leaf():
                    # Calculate average metrics for this clade
                    leaves = node.get_leaves()
                    if leaves:
                        avg_rate = np.mean([getattr(leaf, 'inversion_rate', 0) for leaf in leaves])
                        total_count = sum([getattr(leaf, 'inversion_count', 0) for leaf in leaves])
                        avg_normalized = np.mean([getattr(leaf, 'normalized_score', 0) for leaf in leaves])
                        
                        node.add_features(
                            inversion_count=total_count,
                            inversion_rate=avg_rate,
                            normalized_score=avg_normalized
                        )
            
            return tree
            
        except Exception as e:
            logger.error(f"Tree annotation failed: {e}")
            return tree
    
    def _create_tree_plot(self, tree: 'Tree', output_file: Path, inversion_data: Dict):
        """Create publication-quality tree plot with inversion annotations"""
        
        try:
            if ETE3_AVAILABLE and tree:
                # Use ete3 for tree plotting
                self._create_ete3_tree_plot(tree, output_file, inversion_data)
            else:
                # Fallback to matplotlib
                self._create_matplotlib_tree_plot(inversion_data, output_file)
                
        except Exception as e:
            logger.error(f"Tree plot creation failed: {e}")
            # Create simple fallback
            self._create_matplotlib_tree_plot(inversion_data, output_file)
    
    def _create_ete3_tree_plot(self, tree: 'Tree', output_file: Path, inversion_data: Dict):
        """Create tree plot using ete3"""
        
        from ete3 import TreeStyle, NodeStyle, faces, AttrFace
        
        # Create tree style
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        ts.show_branch_support = True
        ts.mode = "r"  # rectangular
        ts.scale = 200
        
        # Color nodes based on inversion rate
        for node in tree.traverse():
            ns = NodeStyle()
            
            if node.is_leaf():
                # Color leaf nodes by inversion rate
                if hasattr(node, 'normalized_score'):
                    score = node.normalized_score
                    color_intensity = int(255 * (1 - score))  # Higher score = redder
                    ns.bgcolor = f"rgb(255,{color_intensity},{color_intensity})"
                    ns.size = 10
                else:
                    ns.bgcolor = "lightblue"
                    ns.size = 8
                
                # Add inversion count as text
                if hasattr(node, 'inversion_count'):
                    count_face = AttrFace("inversion_count", fsize=8)
                    node.add_face(count_face, column=1, position="branch-right")
            else:
                # Internal nodes
                ns.size = 5
                ns.shape = "square"
                if hasattr(node, 'normalized_score'):
                    score = node.normalized_score
                    color_intensity = int(255 * (1 - score))
                    ns.bgcolor = f"rgb(255,{color_intensity},{color_intensity})"
                else:
                    ns.bgcolor = "gray"
            
            node.set_style(ns)
        
        # Render tree
        tree.render(str(output_file), tree_style=ts, dpi=300)
        logger.info(f"    ✅ ETE3 tree plot created: {output_file}")
    
    def _create_matplotlib_tree_plot(self, inversion_data: Dict, output_file: Path):
        """Create simple tree plot using matplotlib"""
        
        try:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            species_names = list(inversion_data.keys())
            n_species = len(species_names)
            
            if n_species < 2:
                ax.text(0.5, 0.5, "Not enough species for tree", ha='center', va='center')
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()
                return
            
            # Create simple dendrogram-like structure
            from scipy.cluster.hierarchy import dendrogram
            from scipy.spatial.distance import pdist, squareform
            
            # Create distance matrix based on inversion rates
            rates = [inversion_data[sp]['rate_per_mb'] for sp in species_names]
            rate_matrix = np.array(rates).reshape(-1, 1)
            
            if n_species > 1:
                distances = pdist(rate_matrix, metric='euclidean')
                
                # Simple linkage
                linkage_matrix = []
                for i in range(n_species - 1):
                    # Simple clustering
                    linkage_matrix.append([i, i+1, distances[0] if len(distances) > 0 else 1.0, 2])
                
                linkage_matrix = np.array(linkage_matrix)
                
                # Create dendrogram
                dendro = dendrogram(linkage_matrix, labels=species_names, ax=ax)
                
                # Color-code by inversion rates
                for i, species in enumerate(species_names):
                    rate = inversion_data[species]['rate_per_mb']
                    normalized_rate = rate / max(rates) if max(rates) > 0 else 0
                    color = plt.cm.Reds(normalized_rate)
                    
                    # Add colored annotation
                    ax.text(dendro['icoord'][i][1], -0.1, f"{rate:.2f} inv/Mb", 
                           rotation=45, ha='right', va='top', color=color, fontweight='bold')
            
            ax.set_title('Bibionidae Phylogeny with Inversion Rates', fontsize=14, fontweight='bold')
            ax.set_ylabel('Distance')
            ax.set_xlabel('Species')
            
            # Add colorbar
            sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=0, vmax=max(rates)))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.6)
            cbar.set_label('Inversion Rate (per Mb)', rotation=270, labelpad=15)
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"    ✅ Matplotlib tree plot created: {output_file}")
            
        except Exception as e:
            logger.error(f"Matplotlib tree plot failed: {e}")
            # Create very simple plot
            fig, ax = plt.subplots(figsize=(10, 6))
            species_names = list(inversion_data.keys())
            rates = [inversion_data[sp]['rate_per_mb'] for sp in species_names]
            
            bars = ax.bar(species_names, rates, color=plt.cm.Reds([r/max(rates) if max(rates) > 0 else 0 for r in rates]))
            ax.set_title('Inversion Rates by Species', fontsize=14, fontweight='bold')
            ax.set_ylabel('Inversions per Mb')
            ax.set_xlabel('Species')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
    
    def _create_tree_heatmap(self, inversion_data: Dict, output_file: Path):
        """Create heatmap showing inversion metrics across species"""
        
        try:
            # Prepare data for heatmap
            species_names = list(inversion_data.keys())
            metrics = ['total_inversions', 'rate_per_mb', 'normalized_score', 'genome_size']
            
            heatmap_data = []
            for metric in metrics:
                row = [inversion_data[sp].get(metric, 0) for sp in species_names]
                # Normalize each metric to 0-1 scale
                max_val = max(row) if max(row) > 0 else 1
                normalized_row = [val / max_val for val in row]
                heatmap_data.append(normalized_row)
            
            heatmap_df = pd.DataFrame(heatmap_data, 
                                    index=['Total Inversions', 'Rate per Mb', 'Normalized Score', 'Genome Size'],
                                    columns=species_names)
            
            # Create heatmap
            fig, ax = plt.subplots(figsize=(10, 6))
            
            sns.heatmap(heatmap_df, annot=True, fmt='.2f', cmap='RdYlBu_r', 
                       cbar_kws={'label': 'Normalized Value'}, ax=ax)
            
            ax.set_title('Inversion Metrics Heatmap', fontsize=14, fontweight='bold')
            ax.set_xlabel('Species')
            ax.set_ylabel('Metrics')
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"    ✅ Tree heatmap created: {output_file}")
            
        except Exception as e:
            logger.error(f"Tree heatmap creation failed: {e}")


# Integration function for existing visualization system
def create_publication_plots(ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame, 
                           synteny_df: pd.DataFrame, registry, config: Dict, 
                           output_dir: Path, species_stats: Dict = None):
    """
    Main integration function for publication plots
    Can be called from existing visualization system
    """
    
    logger.info("🎨 Creating publication-quality plots...")
    
    try:
        # Create publication config if not exists
        pub_config = config.get('publication_config', {
            'publication_suite': {'enabled': True},
            'synteny_visualization': {'enabled': True, 'options': {}},
            'tree_annotation': {'enabled': True}
        })
        
        # Update config
        updated_config = config.copy()
        updated_config['publication_config'] = pub_config
        
        # Create plot generator
        plot_generator = PublicationPlotGenerator(registry, updated_config)
        
        # Prepare all_results format (if coming from single analysis)
        if not isinstance(ortholog_df, dict):
            # Single pair analysis - convert to multi-species format
            all_results = {
                'single_pair': {
                    'species_pair': ('Species_A', 'Species_B'),
                    'ortholog_df': ortholog_df,
                    'inversion_df': inversion_df,
                    'synteny_df': synteny_df,
                    'full_results': {
                        'ortholog_df': ortholog_df,
                        'inversion_df': inversion_df,
                        'synteny_df': synteny_df
                    }
                }
            }
            
            # Create dummy species stats if not provided
            if species_stats is None:
                species_stats = {
                    'Species_A': {
                        'quality': {
                            'metrics': {'total_length': 100000000},
                            'quality_score': 0.8
                        }
                    },
                    'Species_B': {
                        'quality': {
                            'metrics': {'total_length': 100000000},
                            'quality_score': 0.8
                        }
                    }
                }
        else:
            # Multi-species format
            all_results = ortholog_df  # ortholog_df is actually all_results in this case
        
        # Create publication suite
        results = plot_generator.create_publication_suite(all_results, species_stats, output_dir)
        
        logger.info("✅ Publication plots completed")
        return results
        
    except Exception as e:
        logger.error(f"Publication plots failed: {e}")
        import traceback
        traceback.print_exc()
        return {}


# Convenience functions for individual plot types
def create_curved_synteny_plot(ortholog_df: pd.DataFrame, inversion_df: pd.DataFrame,
                             species1: str, species2: str, output_dir: Path,
                             registry=None, config: Dict = None) -> Optional[Path]:
    """Create a single curved synteny plot"""
    
    if config is None:
        config = {'publication_config': {'synteny_visualization': {'enabled': True}}}
    
    plot_generator = PublicationPlotGenerator(registry, config)
    
    result = plot_generator._create_single_fallback_plot(
        ortholog_df, inversion_df, species1, species2, output_dir
    )
    
    return result


def create_annotated_phylogeny(all_results: Dict, species_stats: Dict, 
                             tree_path: str, output_dir: Path,
                             registry=None) -> Dict:
    """Create annotated phylogenetic tree"""
    
    config = {
        'publication_config': {
            'tree_annotation': {
                'enabled': True,
                'source_tree_path': tree_path
            }
        }
    }
    
    plot_generator = PublicationPlotGenerator(registry, config)
    
    return plot_generator._create_annotated_phylogeny(all_results, species_stats, output_dir)