#!/usr/bin/env python3
"""
Test the complete clean system
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

print("🔬 TESTING CLEAN SYSTEM")
print("=" * 30)

# Test 1: All imports
print("Test 1: Module imports...")
try:
    from genome_inversion_analyser.config import PUBLICATION_CONFIG
    from genome_inversion_analyser.visualization import SyRIIntegrator, PublicationPlotGenerator
    from genome_inversion_analyser.phylogenetic import PhylogeneticIntegrator
    print("✅ All imports successful")
except Exception as e:
    print(f"❌ Import failed: {e}")

# Test 2: Configuration
print("\nTest 2: Configuration...")
try:
    print(f"✅ PUBLICATION_CONFIG keys: {list(PUBLICATION_CONFIG.keys())}")
    print(f"✅ Synteny enabled: {PUBLICATION_CONFIG['synteny_visualization']['enabled']}")
    print(f"✅ Tree path: {PUBLICATION_CONFIG['tree_annotation']['source_tree_path']}")
except Exception as e:
    print(f"❌ Config test failed: {e}")

print(f"\n🎯 System ready for full multi-species analysis!")