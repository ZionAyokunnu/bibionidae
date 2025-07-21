#!/usr/bin/env python3
"""Simple analysis of syngraph results"""

import pandas as pd
import os

def check_files():
    """Check what files exist"""
    files_to_check = [
        "conservative_results.rearrangements.tsv",
        "sensitive_results.rearrangements.tsv", 
        "conservative_final.tsv",
        "sensitive_final.tsv"
    ]
    
    print("=== FILE CHECK ===")
    for file in files_to_check:
        if os.path.exists(file):
            size = os.path.getsize(file)
            print(f"✓ {file} exists ({size} bytes)")
        else:
            print(f"✗ {file} missing")
    print()

def analyze_rearrangements():
    """Analyze rearrangement results"""
    print("=== REARRANGEMENT ANALYSIS ===")
    
    for analysis in ["conservative", "sensitive"]:
        filename = f"{analysis}_results.rearrangements.tsv"
        print(f"\n{analysis.upper()} RESULTS:")
        
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            print(f"Total lines: {len(lines)}")
            print("Content:")
            for i, line in enumerate(lines):
                print(f"  Line {i+1}: {line.strip()}")
            
            # Count actual rearrangements (exclude header)
            rearrangements = len(lines) - 1 if lines[0].startswith('#') else len(lines)
            print(f"Actual rearrangements: {rearrangements}")
        else:
            print(f"File {filename} not found!")

def main():
    check_files()
    analyze_rearrangements()
    
    print("\n=== INTERPRETATION ===")
    print("If both analyses show 0 rearrangements, this means:")
    print("1. The two species have very similar chromosome organization")
    print("2. No detectable fusions or fissions occurred")
    print("3. This suggests stable genome structure between species")
    print("\nThis is actually a meaningful biological result!")

if __name__ == "__main__":
    main()