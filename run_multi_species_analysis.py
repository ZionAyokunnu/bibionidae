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
        "output_directory": "./bibionidae_multi_species_analysis",
        "species_data": [
            {
                "name": "Bibio_marci",
                "fasta": "bibionidae/your_bibio_marci_genome.fasta",
                "busco": "bibionidae/busco-tables/your_bibio_marci_busco.tsv",
                "metadata": {
                    "common_name": "Bibio marci",
                    "lineage": "diptera_odb10",
                    "family": "Bibionidae"
                }
            },
            {
                "name": "Species_B", 
                "fasta": "bibionidae/your_species_B_genome.fasta",
                "busco": "bibionidae/busco-tables/your_species_B_busco.tsv",
                "metadata": {
                    "common_name": "Species B",
                    "lineage": "diptera_odb10",
                    "family": "Bibionidae"
                }
            },
            {
                "name": "Species_C",
                "fasta": "bibionidae/your_species_C_genome.fasta",
                "busco": "bibionidae/busco-tables/your_species_C_busco.tsv",
                "metadata": {
                    "common_name": "Species C",
                    "lineage": "diptera_odb10",
                    "family": "Bibionidae"
                }
            },
            {
                "name": "Species_D",
                "fasta": "bibionidae/your_species_D_genome.fasta",
                "busco": "bibionidae/busco-tables/your_species_D_busco.tsv",
                "metadata": {
                    "common_name": "Species D",
                    "lineage": "diptera_odb10",
                    "family": "Bibionidae"
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
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(example_config, f, indent=2)
    
    print(f"✅ Example configuration created: {output_file}")
    print("📝 Edit this file with your actual file paths and run:")
    print(f"python run_multi_species_analysis.py --config {output_file}")
    print("")
    print("📂 Your file structure should look like:")
    print("  bibionidae/")
    print("  ├── busco-tables/")
    print("  │   ├── species1_busco.tsv")
    print("  │   └── species2_busco.tsv")
    print("  ├── species1_genome.fasta")
    print("  └── species2_genome.fasta")


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
        
        logger.info("✅ Configuration validated successfully")
        logger.info(f"📊 Species to analyze: {len(config['species_data'])}")
        
        # For now, just show what would be analyzed
        print("\n📋 ANALYSIS PLAN:")
        for i, species in enumerate(config['species_data'], 1):
            print(f"  {i}. {species['name']}")
            print(f"     FASTA: {species['fasta']}")
            print(f"     BUSCO: {species['busco']}")
        
        print(f"\n📂 Output directory: {config['output_directory']}")
        print("\n⚠️  FULL PIPELINE NOT YET IMPLEMENTED")
        print("🔧 Working on pipeline import fixes...")
        
        return True
        
    except Exception as e:
        logger.error(f"Multi-species analysis failed: {e}")
        return False


def print_pipeline_summary(results: Dict, output_dir: str):
    """Print comprehensive pipeline summary"""
    print("📊 Pipeline summary would go here")


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
        print("\n🎉 Configuration validated and analysis plan created!")
        sys.exit(0)
    else:
        print("\n❌ Configuration validation failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()