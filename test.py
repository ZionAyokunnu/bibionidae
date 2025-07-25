#!/usr/bin/env python3
"""
CORRECTED TEST - runs from the right directory
"""

import subprocess
import pandas as pd
from pathlib import Path
import os

print("🎨 Testing Charlotte's R-based synteny_plotter (CORRECTED)...")

# Check if the R script exists
synteny_dir = Path("genome_inversion_analyser/external_tools/synteny_plotter")
r_script_path = synteny_dir / "scripts" / "generate_synteny_plot.R"

if not r_script_path.exists():
    print(f"❌ R script not found: {r_script_path}")
    exit(1)

print(f"✅ R script found: {r_script_path}")

# Check if R is available
try:
    r_result = subprocess.run(['Rscript', '--version'], capture_output=True, text=True)
    print(f"✅ Rscript available")
except FileNotFoundError:
    print(f"❌ Rscript not found. Please install R first.")
    exit(1)

# Test with the provided example data first
print(f"\n📊 Testing with provided example data...")

test_data_dir = synteny_dir / "test_data"
if not test_data_dir.exists():
    print(f"❌ Test data directory not found: {test_data_dir}")
    exit(1)

# Check what's in test_data
print(f"📁 Test data contents:")
try:
    for item in test_data_dir.iterdir():
        if item.is_file():
            print(f"  • {item.name}")
except Exception as e:
    print(f"❌ Cannot list test_data contents: {e}")

# Create output directory
output_dir = Path("synteny_test_output")
output_dir.mkdir(exist_ok=True)

# CORRECTED: Use relative paths from synteny_plotter directory
cmd = [
    'Rscript', 'scripts/generate_synteny_plot.R',
    '-busco1', 'test_data/Melitaea_cinxia.tsv',
    '-chrom1', 'test_data/Melitaea_cinxia_info.tsv',
    '-busco_list', 'test_data/Vanessa_cardui.tsv', 
    '-chrom_list', 'test_data/Vanessa_cardui_info.tsv',
    '-o', str(output_dir.absolute() / 'example_test')  # Absolute path for output
]

print(f"\n🔧 Running CORRECTED command:")
print(f"Working directory: {synteny_dir.absolute()}")
print(f"Command: {' '.join(cmd)}")

try:
    # CRITICAL: Run from synteny_plotter directory with relative paths
    result = subprocess.run(
        cmd, 
        capture_output=True, 
        text=True, 
        timeout=60, 
        cwd=str(synteny_dir.absolute())  # This is the key fix!
    )
    
    print(f"\nReturn code: {result.returncode}")
    print(f"STDOUT:\n{result.stdout}")
    
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")
    
    if result.returncode == 0:
        print("✅ Example test succeeded!")
        
        # Check for output files
        output_files = list(output_dir.glob('example_test*'))
        if output_files:
            print(f"✅ Output files created:")
            for f in output_files:
                print(f"  - {f} ({f.stat().st_size} bytes)")
        else:
            print(f"⚠️ Script succeeded but no output files found in {output_dir}")
            
            # Check if files were created in the synteny_plotter directory instead
            synteny_output_files = list(synteny_dir.glob('example_test*'))
            if synteny_output_files:
                print(f"✅ Output files found in synteny_plotter directory:")
                for f in synteny_output_files:
                    print(f"  - {f} ({f.stat().st_size} bytes)")
    else:
        print("❌ Example test failed")
        
        # Check R dependencies if the error suggests missing packages
        if "there is no package called" in result.stderr:
            print(f"\n📦 Missing R packages detected!")
            print(f"Install missing packages in R:")
            print(f"  install.packages(c('dplyr', 'argparse', 'ggplot2'))")
        
except subprocess.TimeoutExpired:
    print(f"❌ Example test timed out")
except Exception as e:
    print(f"❌ Error running example test: {e}")

print(f"\n🎯 Test completed!")

# Manual test instructions
print(f"\n📝 MANUAL TEST (if needed):")
print(f"cd {synteny_dir.absolute()}")
print(f"Rscript scripts/generate_synteny_plot.R \\")
print(f"  -busco1 test_data/Melitaea_cinxia.tsv \\")
print(f"  -chrom1 test_data/Melitaea_cinxia_info.tsv \\")
print(f"  -busco_list test_data/Vanessa_cardui.tsv \\")
print(f"  -chrom_list test_data/Vanessa_cardui_info.tsv \\")
print(f"  -o test_manual")

print(f"\n💡 Key insight: R script must be run from synteny_plotter directory!")