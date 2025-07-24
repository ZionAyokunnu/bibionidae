#!/usr/bin/env python3
"""
Multi-Species Genome Inversion Analysis - Complete Pipeline
Usage: python run_multi_species_analysis.py --config species_config.json
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, List

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_example_config(output_file: str = "species_config.json"):
    """Create an example configuration file"""
    
    example_config = {
        "output_directory": "./multi_species_analysis_output",
        "species_data": [
            {
                "name": "Species_A",
                "fasta": "/path/to/species_A_genome.fasta",
                "busco": "/path/to/species_A_busco_results.tsv",
                "metadata": {
                    "common_name": "Species A",
                    "lineage": "diptera_odb10",
                    "assembly_date": "2024-01-01"
                }
            },
            {
                "name": "Species_B", 
                "fasta": "/path/to/species_B_genome.fasta",
                "busco": "/path/to/species_B_busco_results.tsv",
                "metadata": {
                    "common_name": "Species B",
                    "lineage": "diptera_odb10",
                    "assembly_date": "2024-01-15"
                }
            },
            {
                "name": "Species_C",
                "fasta": "/path/to/species_C_genome.fasta", 
                "busco": "/path/to/species_C_busco_results.tsv",
                "metadata": {
                    "common_name": "Species C",
                    "lineage": "diptera_odb10",
                    "assembly_date": "2024-02-01"
                }
            }
        ],
        "analysis_settings": {
            "alignment_strategy": "hybrid_enhanced",
            "similarity_threshold": 0.85,
            "confidence_threshold": 0.8,
            "min_block_size": 3,
            "enable_contextual_analysis": True,
            "enable_phylogenetic_analysis": True,
            "create_syri_visualizations": True
        },
        "phylogenetic_settings": {
            "distance_methods": ["jaccard", "manhattan", "euclidean", "inversion_rate"],
            "linkage_methods": ["ward", "average", "complete"],
            "min_species_for_trees": 3
        },
        "visualization_settings": {
            "create_standard_plots": True,
            "create_syri_style_plots": True,
            "plot_formats": ["png", "pdf"],
            "dpi": 300
        },
        "optional_data": {
            "functional_annotations": {
                "enabled": False,
                "gff_files": {
                    "Species_A": "/path/to/species_A_genes.gff3",
                    "Species_B": "/path/to/species_B_genes.gff3"
                }
            },
            "repeat_elements": {
                "enabled": False,
                "repeat_files": {
                    "Species_A": "/path/to/species_A_repeats.bed",
                    "Species_B": "/path/to/species_B_repeats.bed"
                }
            }
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(example_config, f, indent=2)
    
    print(f"Example configuration created: {output_file}")
    print("Edit this file with your actual file paths and run:")
    print(f"python run_multi_species_analysis.py --config {output_file}")


def validate_config(config: Dict) -> bool:
    """Validate the configuration file"""
    
    required_keys = ['output_directory', 'species_data']
    
    for key in required_keys:
        if key not in config:
            logger.error(f"Missing required configuration key: {key}")
            return False
    
    # Validate species data
    species_data = config['species_data']
    if not isinstance(species_data, list) or len(species_data) < 2:
        logger.error("Need at least 2 species for multi-species analysis")
        return False
    
    # Check each species has required fields
    for i, species in enumerate(species_data):
        required_species_keys = ['name', 'fasta', 'busco']
        for key in required_species_keys:
            if key not in species:
                logger.error(f"Species {i}: missing required key '{key}'")
                return False
        
        # Check if files exist
        fasta_path = Path(species['fasta'])
        busco_path = Path(species['busco'])
        
        if not fasta_path.exists():
            logger.error(f"FASTA file not found: {fasta_path}")
            return False
        
        if not busco_path.exists():
            logger.error(f"BUSCO file not found: {busco_path}")
            return False
    
    logger.info("Configuration validation passed")
    return True


def run_multi_species_analysis(config_file: str):
    """Run the complete multi-species analysis pipeline"""
    
    logger.info("🧬 STARTING MULTI-SPECIES GENOME INVERSION ANALYSIS")
    logger.info("=" * 80)
    
    try:
        # Load configuration
        logger.info(f"Loading configuration from: {config_file}")
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Validate configuration
        if not validate_config(config):
            logger.error("Configuration validation failed")
            return False
        
        # Initialize pipeline
        output_dir = config['output_directory']
        pipeline = MultiSpeciesPipeline(output_dir, config.get('analysis_settings', {}))
        
        logger.info(f"Pipeline initialized: {output_dir}")
        
        # Add all species
        logger.info("Adding species to pipeline...")
        species_data = config['species_data']
        
        for species in species_data:
            success = pipeline.add_species(
                species['name'],
                species['fasta'],
                species['busco'],
                species.get('metadata', {})
            )
            
            if success:
                logger.info(f"✅ Added species: {species['name']}")
            else:
                logger.error(f"❌ Failed to add species: {species['name']}")
                return False
        
        # Add optional functional annotations
        if config.get('optional_data', {}).get('functional_annotations', {}).get('enabled', False):
            logger.info("Adding functional annotations...")
            gff_files = config['optional_data']['functional_annotations'].get('gff_files', {})
            
            for species_name, gff_file in gff_files.items():
                if Path(gff_file).exists():
                    pipeline.context_analyzer.add_functional_annotations(gff_file, 'genes')
                    logger.info(f"Added functional annotations for {species_name}")
        
        # Add optional repeat data
        if config.get('optional_data', {}).get('repeat_elements', {}).get('enabled', False):
            logger.info("Adding repeat element data...")
            repeat_files = config['optional_data']['repeat_elements'].get('repeat_files', {})
            
            for species_name, repeat_file in repeat_files.items():
                if Path(repeat_file).exists():
                    pipeline.context_analyzer.add_repeat_data(repeat_file, 'bed')
                    logger.info(f"Added repeat data for {species_name}")
        
        # Run complete pipeline
        logger.info("🚀 RUNNING COMPLETE MULTI-SPECIES PIPELINE")
        logger.info("=" * 80)
        
        pipeline_results = pipeline.run_complete_pipeline(
            include_functional_annotations=config.get('optional_data', {}).get('functional_annotations', {}).get('enabled', False),
            include_repeat_analysis=config.get('optional_data', {}).get('repeat_elements', {}).get('enabled', False),
            create_visualizations=config.get('visualization_settings', {}).get('create_syri_style_plots', True)
        )
        
        # Print comprehensive summary
        print_pipeline_summary(pipeline_results, output_dir)
        
        logger.info("🎉 MULTI-SPECIES ANALYSIS COMPLETED SUCCESSFULLY!")
        return True
        
    except Exception as e:
        logger.error(f"Multi-species analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def print_pipeline_summary(results: Dict, output_dir: str):
    """Print comprehensive pipeline summary"""
    
    print("\n" + "=" * 80)
    print("🧬 MULTI-SPECIES GENOME INVERSION ANALYSIS - COMPLETE RESULTS")
    print("=" * 80)
    
    # Species processing summary
    species_summary = results.get('summary', {}).get('species_processing_summary', {})
    print(f"\n📊 SPECIES PROCESSING:")
    print(f"  • Successful species: {species_summary.get('successful_species', 0)}")
    print(f"  • Failed species: {species_summary.get('failed_species', 0)}")
    print(f"  • Total genes analyzed: {species_summary.get('total_genes_analyzed', 0):,}")
    
    # Pairwise analysis summary
    pairwise_summary = results.get('summary', {}).get('pairwise_analysis_summary', {})
    print(f"\n🔄 PAIRWISE COMPARISONS:")
    print(f"  • Successful pairs: {pairwise_summary.get('successful_pairs', 0)}")
    print(f"  • Failed pairs: {pairwise_summary.get('failed_pairs', 0)}")
    print(f"  • Total inversions detected: {pairwise_summary.get('total_inversions_detected', 0):,}")
    
    # Phylogenetic analysis summary
    phylo_summary = results.get('summary', {}).get('phylogenetic_summary', {})
    print(f"\n🌳 PHYLOGENETIC ANALYSIS:")
    print(f"  • Trees constructed: {phylo_summary.get('trees_constructed', 0)}")
    print(f"  • Distance methods: {', '.join(phylo_summary.get('distance_methods_used', []))}")
    print(f"  • Phylogenetic signal detected: {phylo_summary.get('phylogenetic_signal_detected', False)}")
    
    # Visualization summary
    viz_summary = results.get('summary', {}).get('visualization_summary', {})
    print(f"\n📈 VISUALIZATIONS:")
    print(f"  • Visualizations created: {viz_summary.get('visualizations_created', 0)}")
    print(f"  • Failed visualizations: {viz_summary.get('failed_visualizations', 0)}")
    
    # Output structure
    print(f"\n📁 OUTPUT STRUCTURE:")
    print(f"  {output_dir}/")
    print(f"  ├── data/                    # Individual species data")
    print(f"  ├── exports/                 # Standardized format exports")
    print(f"  │   ├── bed/                 # Genome browser files")
    print(f"  │   ├── json/                # Machine-readable data")
    print(f"  │   └── phylogenetic/        # Phylogenetic analysis files")
    print(f"  ├── visualizations/          # All plots and visualizations")
    print(f"  │   ├── Species_A_vs_B/      # Pairwise comparisons")
    print(f"  │   ├── Species_A_vs_C/")
    print(f"  │   └── phylogenetic_trees/  # Phylogenetic visualizations")
    print(f"  ├── reports/                 # Generated reports")
    print(f"  ├── registry.json            # Complete file registry")
    print(f"  └── file_manifest.json       # File manifest")
    
    # Key files for downstream analysis
    print(f"\n🎯 KEY OUTPUT FILES:")
    print(f"  📊 file_manifest.json        → Complete analysis overview")
    print(f"  🧬 exports/bed/*              → Load into genome browsers (IGV, UCSC)")
    print(f"  🌳 exports/phylogenetic/*     → Phylogenetic analysis results")
    print(f"  📈 visualizations/*/          → SyRI-style comparative plots")
    print(f"  📋 reports/*                  → Publication-ready summaries")
    
    # Next steps
    print(f"\n🚀 NEXT STEPS:")
    print(f"  1. 📖 Check reports/ for comprehensive analysis summaries")
    print(f"  2. 🌐 Load BED files into genome browsers for visual exploration")
    print(f"  3. 🌳 Use phylogenetic results for evolutionary analysis")
    print(f"  4. 📊 Explore visualizations/ for publication-quality plots")
    print(f"  5. 🔄 Use registry.json for complete data provenance")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Multi-Species Genome Inversion Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create example configuration
  python run_multi_species_analysis.py --create-config
  
  # Run analysis with configuration file
  python run_multi_species_analysis.py --config species_config.json
  
  # Run with custom output directory
  python run_multi_species_analysis.py --config species_config.json --output ./my_analysis
        """
    )
    
    parser.add_argument(
        '--config', 
        type=str,
        help='Configuration file (JSON format)'
    )
    
    parser.add_argument(
        '--create-config',
        action='store_true',
        help='Create an example configuration file'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        help='Override output directory from config file'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create example config
    if args.create_config:
        config_file = "species_config.json"
        create_example_config(config_file)
        return
    
    # Validate arguments
    if not args.config:
        parser.error("Must specify --config file or use --create-config")
    
    if not Path(args.config).exists():
        parser.error(f"Configuration file not found: {args.config}")
    
    # Run analysis
    success = run_multi_species_analysis(args.config)
    
    if success:
        print("\n🎉 Multi-species analysis completed successfully!")
        sys.exit(0)
    else:
        print("\n❌ Multi-species analysis failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()