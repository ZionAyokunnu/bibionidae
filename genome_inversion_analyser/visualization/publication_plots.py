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
import matplotlib.cm as cm
from matplotlib.colors import Normalize

logger = logging.getLogger(__name__)

try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    logger.warning("ete3 not available - tree manipulation will be limited")

try:
    from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace
    ETE3_AVAILABLE = True
    logger.info("✅ ete3 fully available")
except ImportError as e:
    logger.warning(f"⚠️ ete3 import issue: {e}")
    try:
        from ete3 import Tree
        ETE3_AVAILABLE = True
        logger.info("✅ ete3 Tree available (limited features)")
    except ImportError:
        ETE3_AVAILABLE = False
        logger.warning("❌ ete3 not available - tree manipulation will be limited")


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
        
        logger.info("🎨 Creating publication-quality visualization suite...")
        
        pub_dir = output_dir / 'publication_plots'
        pub_dir.mkdir(exist_ok=True)
        
        results = {}
        
        # Check if publication suite is enabled
        if not self.pub_config.get('publication_suite', {}).get('enabled', False):
            logger.info("Publication suite disabled in config")
            return results  # Return empty results if disabled
        
        # 1. Curved synteny plots
        if self.pub_config.get('synteny_visualization', {}).get('enabled', False):
            logger.info("  📊 Starting synteny visualization...")
            synteny_results = self._create_synteny_plots(all_results, pub_dir)
            results['synteny_plots'] = synteny_results
            
        # 2. Annotated phylogenetic tree
        if self.pub_config.get('tree_annotation', {}).get('enabled', False):
            logger.info("  🌳 Starting tree annotation...")
            tree_results = self.create_annotated_phylogeny(all_results, species_stats, pub_dir)
            results['annotated_trees'] = tree_results
            
        # 3. BUSCO phylogenetic tree
        if self.pub_config.get('busco_phylogeny', {}).get('enabled', False):
            logger.info("  🧬 Starting BUSCO phylogenetic tree...")
            busco_tree_results = self.create_busco_phylogenetic_tree(all_results, species_stats, pub_dir)
            results['busco_phylogenetic_tree'] = busco_tree_results
            
        # 4. Register all outputs
        for plot_type, plot_data in results.items():
            if isinstance(plot_data, dict):
                for plot_name, plot_path in plot_data.items():
                    if plot_path and Path(plot_path).exists():
                        self.registry.register_file(
                            f'publication_{plot_type}_{plot_name}',
                            plot_path,
                            file_type='plot',
                            description=f'Publication-quality {plot_type}: {plot_name}'
                        )

        logger.info(f"✅ Publication suite completed: {len(results)} plot types created")
        return results
    
    def _create_synteny_plots(self, all_results: Dict, output_dir: Path) -> Dict:
        """Create publication-quality synteny plots using synteny_plotter"""
        
        logger.info("  📊 Creating curved synteny plots...")
        
        synteny_dir = output_dir / 'synteny_curves'
        synteny_dir.mkdir(exist_ok=True)
        
        # CORRECTED: Check for R script, not Python script
        synteny_plotter_dir = self.pub_config.get('external_tools', {}).get('synteny_plotter')
        
        if synteny_plotter_dir:
            synteny_plotter_path = Path(synteny_plotter_dir).resolve() / 'scripts' / 'generate_synteny_plot.R'
            logger.info(f"DEBUG: Looking for R script at: {synteny_plotter_path.absolute()}")
        else:
            synteny_plotter_path = None
        
        if not synteny_plotter_path or not synteny_plotter_path.exists():
            logger.warning(f"Synteny plotter R script not found: {synteny_plotter_path}")
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
        """Convert ortholog data to synteny_plotter BUSCO TSV format"""
        
        # Get the synteny_plotter directory for relative paths
        synteny_plotter_dir = Path(self.pub_config['external_tools']['synteny_plotter'])
        
        # Create temp directory INSIDE synteny_plotter for relative paths
        temp_dir = synteny_plotter_dir / f'temp_{species1}_vs_{species2}'
        temp_dir.mkdir(exist_ok=True)
        
        # Create BUSCO TSV files in the EXACT format the R script expects
        busco1_data = []
        busco2_data = []
        
        for _, row in ortholog_df.iterrows():
            # Format for reference species (R prefix) - EXACT format for R script
            busco1_data.append([
                row['busco_id'],           # busco
                'Complete',                # status  
                row['first_chr'],          # chrR (will be renamed by R script)
                int(row['first_start']),   # Rstart
                int(row['first_end']),     # Rend
                row.get('first_strand', '+')  # Rstrand
            ])
            
            # Format for query species (Q prefix) - EXACT format for R script
            busco2_data.append([
                row['busco_id'],           # busco
                'Complete',                # status
                row['second_chr'],         # chrQ (will be renamed by R script)
                int(row['second_start']),  # Qstart
                int(row['second_end']),    # Qend
                row.get('second_strand', '+')  # Qstrand
            ])
        
        # Write BUSCO TSV files - NO HEADERS, just data (R script adds its own headers)
        busco1_file = temp_dir / f'{species1}.tsv'
        busco2_file = temp_dir / f'{species2}.tsv'
        
        with open(busco1_file, 'w') as f:
            for row in busco1_data:
                f.write('\t'.join(map(str, row)) + '\n')

        with open(busco2_file, 'w') as f:
            for row in busco2_data:
                f.write('\t'.join(map(str, row)) + '\n')
        
        # Create chromosome info files (chr, length, order, direction, invert)
        chrom1_info = {}
        chrom2_info = {}
        
        for _, row in ortholog_df.iterrows():
            chr1 = row['first_chr']
            chr2 = row['second_chr']
            
            if chr1 not in chrom1_info:
                max_end = ortholog_df[ortholog_df['first_chr'] == chr1]['first_end'].max()
                chrom1_info[chr1] = {
                    'chr': chr1,
                    'length': max_end,
                    'order': len(chrom1_info) + 1,
                    'direction': 1,
                    'invert': 'FALSE'  # R script needs this column
                }
                
            if chr2 not in chrom2_info:
                max_end = ortholog_df[ortholog_df['second_chr'] == chr2]['second_end'].max()
                chrom2_info[chr2] = {
                    'chr': chr2,
                    'length': max_end,
                    'order': len(chrom2_info) + 1,
                    'direction': 1,
                    'invert': 'FALSE'  # R script needs this column
                }
        
        # Write chromosome info files WITH headers
        chrom1_file = temp_dir / f'{species1}_info.tsv'
        chrom2_file = temp_dir / f'{species2}_info.tsv'
        
        pd.DataFrame(list(chrom1_info.values())).to_csv(chrom1_file, sep='\t', index=False, header=True)
        pd.DataFrame(list(chrom2_info.values())).to_csv(chrom2_file, sep='\t', index=False, header=True)
        
        return {
            'busco1': busco1_file,
            'busco2': busco2_file, 
            'chrom1': chrom1_file,
            'chrom2': chrom2_file,
            'temp_dir': temp_dir
        }
    def _run_synteny_plotter(self, temp_files: Dict, output_plot: Path, 
                        species1: str, species2: str) -> bool:
        """Run Charlotte's R-based synteny_plotter with the correct argument format"""
        
        # Get the path to the synteny_plotter directory
        synteny_plotter_dir = Path(self.pub_config['external_tools']['synteny_plotter'])
        r_script_path = synteny_plotter_dir / 'scripts' / 'generate_synteny_plot.R'
        
        if not r_script_path.exists():
            logger.error(f"R script not found: {r_script_path}")
            return False
        
        # Build R command with CORRECT argument format:
        # -busco1 + -chrom1 for reference species
        # -busco_list + -chrom_list for query species

        cmd = [
            'Rscript', 'scripts/generate_synteny_plot.R',
            '-busco1', f"{temp_files['busco1'].relative_to(synteny_plotter_dir)}",  # Reference species
            '-chrom1', f"{temp_files['chrom1'].relative_to(synteny_plotter_dir)}",  # Reference chr info
            '-busco_list', f"{temp_files['busco2'].relative_to(synteny_plotter_dir)}",  # Query species
            '-chrom_list', f"{temp_files['chrom2'].relative_to(synteny_plotter_dir)}",  # Query chr info
            '-o', str(output_plot.with_suffix('').name)  # Output prefix
        ]
        
        # Add optional parameters
        synteny_config = self.pub_config.get('synteny_visualization', {})
        options = synteny_config.get('options', {})
        
        if options.get('filter_threshold'):
            cmd.extend(['-f', str(options['filter_threshold'])])
        if options.get('gap'):
            cmd.extend(['-g', str(options['gap'])])
        if options.get('alpha'):
            cmd.extend(['-alpha', str(options['alpha'])])
        
        try:
            logger.info(f"Running synteny_plotter: {' '.join(cmd)}")
            logger.info(f"Working directory: {synteny_plotter_dir}")
            
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=str(synteny_plotter_dir)  # CRITICAL: Run from synteny_plotter directory
            )
            
            logger.info(f"Synteny plotter stdout: {result.stdout}")
            
            # Clean up temp files
            if temp_files['temp_dir'].exists():
                shutil.rmtree(temp_files['temp_dir'])
            
            # R script creates output in its working directory, move to expected location
            # Check for various possible output files
            possible_outputs = []
            for ext in ['.png', '.pdf', '.svg', '']:
                possible_outputs.append(synteny_plotter_dir / f"{output_plot.with_suffix('').name}{ext}")
            
            output_created = False
            for possible_output in possible_outputs:
                if possible_output.exists():
                    # Move to expected location
                    final_output = output_plot.with_suffix(possible_output.suffix or '.png')
                    shutil.move(str(possible_output), str(final_output))
                    logger.info(f"Moved output: {possible_output} -> {final_output}")
                    output_created = True
                    break
            
            if not output_created:
                logger.warning("Synteny plotter completed but no output file found")
                # List what files were created
                created_files = list(synteny_plotter_dir.glob(f"{output_plot.with_suffix('').name}*"))
                if created_files:
                    logger.info(f"Files created: {[f.name for f in created_files]}")
            
            return output_created
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Synteny plotter failed: {e.stderr}")
            logger.error(f"Command: {' '.join(cmd)}")
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
            print(f"    No ortholog data for {species1} vs {species2}")
            return None
        
        try:
            print(f"    Creating matplotlib plot for {species1} vs {species2}")
            fig, ax = plt.subplots(figsize=(14, 10))
            
            # Use correct column names directly
            chr1_list = sorted(ortholog_df['first_chr'].unique())
            chr2_list = sorted(ortholog_df['second_chr'].unique())
            
            print(f"    Chromosomes: {len(chr1_list)} vs {len(chr2_list)}")
            
            # Create chromosome positions
            chr1_positions = {chr_name: i for i, chr_name in enumerate(chr1_list)}
            chr2_positions = {chr_name: i for i, chr_name in enumerate(chr2_list)}
            
            y1_base = 0.7
            y2_base = 0.3
            
            # Plot chromosome backbones
            for i, chr_name in enumerate(chr1_list):
                ax.plot([0, 1], [y1_base + i*0.05, y1_base + i*0.05], 'k-', linewidth=3, alpha=0.7)
                display_name = str(chr_name)[:10] + "..." if len(str(chr_name)) > 10 else str(chr_name)
                ax.text(-0.05, y1_base + i*0.05, display_name, ha='right', va='center', fontsize=8)
            
            for i, chr_name in enumerate(chr2_list):
                ax.plot([0, 1], [y2_base - i*0.05, y2_base - i*0.05], 'k-', linewidth=3, alpha=0.7)
                display_name = str(chr_name)[:10] + "..." if len(str(chr_name)) > 10 else str(chr_name)
                ax.text(-0.05, y2_base - i*0.05, display_name, ha='right', va='center', fontsize=8)
            
            # Plot connections (limit to avoid overcrowding)
            connections_plotted = 0
            max_connections = min(200, len(ortholog_df))  # Limit for readability
            
            for _, row in ortholog_df.head(max_connections).iterrows():
                try:
                    chr1_idx = chr1_positions[row['first_chr']]
                    chr2_idx = chr2_positions[row['second_chr']]
                    
                    y1 = y1_base + chr1_idx * 0.05
                    y2 = y2_base - chr2_idx * 0.05
                    
                    # Simple position mapping
                    x1 = 0.1 + 0.8 * (connections_plotted / max_connections)
                    x2 = 0.1 + 0.8 * (connections_plotted / max_connections)
                    
                    # Color by similarity
                    color = cm.get_cmap('viridis')(row['similarity'])
                    alpha = 0.6
                    
                    # Create curved connection
                    x_curve = np.linspace(x1, x2, 20)
                    y_curve = []
                    
                    for x in x_curve:
                        t = (x - x1) / (x2 - x1) if x2 != x1 else 0
                        y = y1 * (1 - t) + y2 * t + 0.1 * np.sin(np.pi * t)
                        y_curve.append(y)
                    
                    ax.plot(x_curve, y_curve, color=color, alpha=alpha, linewidth=1)
                    connections_plotted += 1
                    
                except (KeyError, ValueError) as e:
                    continue
            
            print(f"    Plotted {connections_plotted} connections")
            
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
            sm = cm.ScalarMappable(cmap='viridis', norm=Normalize(vmin=0, vmax=1))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
            cbar.set_label('Similarity', rotation=270, labelpad=15)
            
            plt.tight_layout()
            
            # Save plot
            plot_file = output_dir / f'matplotlib_synteny_{species1}_vs_{species2}.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"    ✅ Saved: {plot_file}")
            return plot_file
            
        except Exception as e:
            print(f"    ❌ Plot creation failed: {e}")
            import traceback
            traceback.print_exc()
            if 'fig' in locals():
                plt.close()
            return None
    
    def create_annotated_phylogeny(self, all_results: Dict, species_stats: Dict, 
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
        """Create tree plot using ete3 with better error handling"""
        
        try:
            from ete3 import TreeStyle, NodeStyle, AttrFace
        except ImportError as e:
            logger.error(f"Cannot import ete3 components: {e}")
            return False

        try:
            # Create tree style
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_length = True
            ts.show_branch_support = True
            ts.mode = "r"  # rectangular
            ts.scale = 200
            
            # Configure tree style
            ts.branch_vertical_margin = 10
            ts.scale = 120
            
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
            return True
            
        except Exception as e:
            logger.error(f"ETE3 tree plot creation failed: {e}")
            return False
    
    def _create_matplotlib_tree_plot(self, inversion_data: Dict, output_file: Path):
        """Create simple tree plot using matplotlib with proper linkage"""
        
        try:
            from scipy.cluster.hierarchy import dendrogram, linkage
            from scipy.spatial.distance import pdist
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            species_names = list(inversion_data.keys())
            n_species = len(species_names)
            
            if n_species < 2:
                ax.text(0.5, 0.5, "Not enough species for tree", ha='center', va='center')
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()
                return
            
            # Create distance matrix based on inversion rates
            rates = [inversion_data[sp]['rate_per_mb'] for sp in species_names]
            
            if n_species == 2:
                # Special case for 2 species
                distance = abs(rates[0] - rates[1]) if rates[0] != rates[1] else 0.1
                linkage_matrix = np.array([[0, 1, distance, 2]])
            else:
                # For 3+ species, create proper distance matrix
                rate_matrix = np.array(rates).reshape(-1, 1)
                distances = pdist(rate_matrix, metric='euclidean')
                linkage_matrix = linkage(distances, method='average')  # Use 'average' instead of 'ward'
            
            # Create dendrogram
            dendro = dendrogram(linkage_matrix, labels=species_names, ax=ax, orientation='top')
            
            # Color-code by inversion rates
            max_rate = max(rates) if rates else 1
            for i, species in enumerate(species_names):
                rate = inversion_data[species]['rate_per_mb']
                normalized_rate = rate / max_rate if max_rate > 0 else 0
                color = cm.get_cmap('Reds')(normalized_rate)
                
                # Find the corresponding x position in the dendrogram
                if i < len(dendro['icoord']):
                    x_pos = (dendro['icoord'][i][1] + dendro['icoord'][i][2]) / 2
                    ax.text(x_pos, -0.1, f"{rate:.2f}", 
                        rotation=45, ha='right', va='top', color=color, fontweight='bold')
            
            ax.set_title('Bibionidae Phylogeny with Inversion Rates', fontsize=14, fontweight='bold')
            ax.set_ylabel('Distance')
            ax.set_xlabel('Species')
            
            # Add colorbar
            sm = cm.ScalarMappable(cmap='Reds', norm=Normalize(vmin=0, vmax=max_rate))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, shrink=0.6)
            cbar.set_label('Inversion Rate (per Mb)', rotation=270, labelpad=15)
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f"    ✅ Matplotlib tree plot created: {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"Matplotlib tree plot failed: {e}")
            # Create simple bar plot fallback
            fig, ax = plt.subplots(figsize=(10, 6))
            species_names = list(inversion_data.keys())
            rates = [inversion_data[sp]['rate_per_mb'] for sp in species_names]
            
            max_rate = max(rates) if rates else 1
            colors = [cm.get_cmap('Reds')(r/max_rate) for r in rates]
            bars = ax.bar(species_names, rates, color=colors)
            
            ax.set_title('Inversion Rates by Species', fontsize=14, fontweight='bold')
            ax.set_ylabel('Inversions per Mb')
            ax.set_xlabel('Species')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            return True
    
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
        
    def create_busco_phylogenetic_tree(self, all_results: Dict, species_stats: Dict, output_dir: Path) -> Dict:
        """Create phylogenetic tree from existing Newick file with configurable node annotations"""
        
        logger.info("  🌳 Creating BUSCO phylogenetic tree from existing file...")
        
        tree_config = self.pub_config.get('tree_annotation', {})
        source_tree = tree_config.get('source_tree_path', 'diptera_clean_20species.newick')
        
        if not Path(source_tree).exists():
            logger.warning(f"Source tree file not found: {source_tree}")
            return {}
        
        try:
            from ete3 import Tree
            
            # Load the tree
            tree = Tree(source_tree)
            logger.info(f"    • Loaded tree with {len(tree)} nodes")
            
            # Get species names
            species_names = list(species_stats.keys())
            logger.info(f"    • Target species: {species_names}")
            
            # Prune tree to target species (if enabled)
            if tree_config.get('prune_to_target_species', True):
                tree = self.prune_tree_to_species(tree, species_names)
                if tree is None:
                    logger.warning("Tree pruning failed")
                    return {}
            
            # Calculate metrics for node annotation
            node_metrics = self.calculate_tree_node_metrics(all_results, species_stats)
            
            # Annotate tree with metrics
            annotated_tree = self.annotate_tree_nodes(tree, node_metrics)
            
            # Create tree plots
            tree_dir = output_dir / 'busco_phylogenetic_trees'
            tree_dir.mkdir(exist_ok=True)
            
            results = {}
            
            # Save annotated tree
            tree_file = tree_dir / 'busco_phylogenetic_tree.newick'
            annotated_tree.write(format=1, outfile=str(tree_file))
            results['newick'] = tree_file
            
            # Create tree plot
            plot_file = tree_dir / 'busco_phylogenetic_tree.png'
            self.create_annotated_tree_plot(annotated_tree, plot_file, node_metrics)
            results['plot'] = plot_file
            
            logger.info(f"    ✅ BUSCO phylogenetic tree created: {len(results)} outputs")
            return results
            
        except ImportError:
            logger.error("ete3 required for phylogenetic tree manipulation")
            return {}
        except Exception as e:
            logger.error(f"BUSCO phylogenetic tree creation failed: {e}")
            return {}

    def calculate_tree_node_metrics(self, all_results: Dict, species_stats: Dict) -> Dict:
        """Calculate metrics for tree node annotation from config"""
        
        tree_config = self.pub_config.get('tree_annotation', {})
        metrics_config = tree_config.get('annotation_metrics', {})
        
        node_metrics = {}
        
        for species_name, stats in species_stats.items():
            metrics = {}
            
            # Default metric: inversion rate per MB
            if metrics_config.get('inversion_rate_per_mb', True):
                genome_size = stats['quality']['metrics'].get('total_length', 1000000)
                species_inversions = sum(1 for pair_results in all_results.values() 
                                    if 'full_results' in pair_results and 
                                    species_name in pair_results['species_pair'] and
                                    len(pair_results['inversion_df']) > 0)
                metrics['inversion_rate_per_mb'] = species_inversions / (genome_size / 1_000_000)
            
            # Inversion count
            if metrics_config.get('inversion_count', True):
                species_inversions = []
                for pair_results in all_results.values():
                    if 'full_results' in pair_results and species_name in pair_results['species_pair']:
                        species_inversions.extend(pair_results['inversion_df'].to_dict('records'))
                metrics['inversion_count'] = len(species_inversions)
            
            # Normalized inversion score
            if metrics_config.get('normalized_inversion_score', True):
                if 'inversion_rate_per_mb' in metrics:
                    all_rates = [node_metrics.get(sp, {}).get('inversion_rate_per_mb', 0) 
                            for sp in species_stats.keys()]
                    max_rate = max(all_rates + [metrics['inversion_rate_per_mb']])
                    metrics['normalized_score'] = metrics['inversion_rate_per_mb'] / max_rate if max_rate > 0 else 0
            
            # Assembly quality
            if metrics_config.get('assembly_quality', False):
                metrics['quality_score'] = stats['quality'].get('quality_score', 0.5)
            
            # Genome size
            if metrics_config.get('genome_size', False):
                metrics['genome_size_mb'] = stats['quality']['metrics'].get('total_length', 0) / 1_000_000
            
            node_metrics[species_name] = metrics
        
        return node_metrics

    def prune_tree_to_species(self, tree: 'Tree', target_species: List[str]) -> Optional['Tree']:
        """Prune tree to target species"""
        try:
            # Get all leaf names
            leaf_names = [leaf.name for leaf in tree.get_leaves()]
            
            # Find matches (allowing partial matches)
            species_in_tree = []
            for species in target_species:
                matches = [leaf for leaf in leaf_names if species in leaf or leaf in species]
                if matches:
                    species_in_tree.append(matches[0])  # Take first match
                else:
                    logger.warning(f"Species {species} not found in tree")
            
            if len(species_in_tree) < 2:
                logger.error(f"Too few species found in tree: {species_in_tree}")
                return None
            
            # Prune tree
            tree.prune(species_in_tree)
            logger.info(f"    • Pruned tree to {len(species_in_tree)} species: {species_in_tree}")
            
            return tree
            
        except Exception as e:
            logger.error(f"Tree pruning failed: {e}")
            return None

    def annotate_tree_nodes(self, tree: 'Tree', node_metrics: Dict) -> 'Tree':
        """Annotate tree nodes with calculated metrics"""
        
        for leaf in tree.get_leaves():
            leaf_name = leaf.name
            
            # Find matching species (allowing partial matches)
            matching_species = None
            for species in node_metrics.keys():
                if species in leaf_name or leaf_name in species:
                    matching_species = species
                    break
            
            if matching_species and matching_species in node_metrics:
                metrics = node_metrics[matching_species]
                
                # Add metrics as node features
                for metric_name, metric_value in metrics.items():
                    setattr(leaf, metric_name, metric_value)
        
        return tree

    # def create_annotated_tree_plot(self, tree: 'Tree', output_file: Path, node_metrics: Dict): # Uses ETE3 GUI components
    #     """Create tree plot with node annotations"""
        
    #     try:
    #         from ete3 import TreeStyle, NodeStyle, TextFace
            
    #         ts = TreeStyle()
    #         ts.show_leaf_name = True
    #         ts.show_branch_length = True
    #         ts.mode = "r"  # rectangular
    #         ts.scale = 120
            
    #         # Color nodes by primary metric (inversion rate)
    #         primary_metric = 'inversion_rate_per_mb'
    #         max_value = max(metrics.get(primary_metric, 0) for metrics in node_metrics.values())
            
    #         for node in tree.traverse():
    #             ns = NodeStyle()
                
    #             if node.is_leaf() and hasattr(node, primary_metric):
    #                 # Color by metric value
    #                 metric_value = getattr(node, primary_metric, 0)
    #                 normalized_value = metric_value / max_value if max_value > 0 else 0
                    
    #                 # Color scale: blue (low) to red (high)
    #                 color_intensity = int(255 * normalized_value)
    #                 ns.bgcolor = f"rgb({color_intensity},0,{255-color_intensity})"
    #                 ns.size = 10
                    
    #                 # Add metric as text
    #                 metric_text = f"{metric_value:.2f}"
    #                 metric_face = TextFace(metric_text, fsize=8)
    #                 node.add_face(metric_face, column=1, position="branch-right")
                
    #             node.set_style(ns)
            
    #         # Render tree
    #         tree.render(str(output_file), tree_style=ts, dpi=300)
    #         logger.info(f"    ✅ Tree plot created: {output_file}")
            
    #     except Exception as e:
    #         logger.error(f"Tree plot creation failed: {e}")
            
    def create_annotated_tree_plot(self, tree, output_file, node_metrics):
        """Create tree plot using ete3 render without GUI components"""
        try:
            # Use basic ete3 render - NO GUI components
            tree.render(str(output_file), w=800, h=600, dpi=300)
            logger.info(f" • BUSCO tree rendered: {output_file}")
            
        except Exception as e:
            logger.error(f"ete3 render failed: {e}, falling back to matplotlib")
            # Fallback to matplotlib
            import matplotlib.pyplot as plt
            import numpy as np
            
            fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            
            leaves = tree.get_leaves()
            species_names = [leaf.name.replace('_', ' ') for leaf in leaves]
            
            for i, name in enumerate(species_names):
                ax.text(1.0, i, name, ha='left', va='center', fontsize=11, weight='bold')
                ax.plot([0.8, 1.0], [i, i], 'k-', linewidth=2)
                
                # Add inversion count if available
                inversions = getattr(leaves[i], 'inversions', 0)
                ax.text(1.5, i, f"({inversions} inv)", ha='left', va='center', 
                    fontsize=9, color='red')
            
            ax.set_xlim(0, 2.0)
            ax.set_ylim(-0.5, len(leaves) - 0.5)
            ax.set_title('BUSCO Phylogenetic Tree', fontsize=14, weight='bold')
            ax.axis('off')
            
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            logger.info(f" • Matplotlib fallback tree saved: {output_file}")

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
        
        return plot_generator.create_annotated_phylogeny(all_results, species_stats, output_dir)
