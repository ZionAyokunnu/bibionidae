#!/usr/bin/env python3
import os

def format_busco_file(input_path, output_path):
    """Format BUSCO file for syngraph"""
    genes = []
    
    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                busco_id, status, sequence, start, end = parts[:5]
                
                if status == 'Complete' and start != 'N/A' and end != 'N/A':
                    genes.append(f"{busco_id}\t{sequence}\t{start}\t{end}")
    
    with open(output_path, 'w') as f:
        for gene in genes:
            f.write(gene + '\n')
    
    return len(genes)

# Process files
os.makedirs('syngraph_input', exist_ok=True)

files = {
    'Bibio_marci': 'Bibio_marci/full_table.tsv',
    'Dilophus_febrilis': 'Dilophus_febrilis/full_table.tsv'
}

for taxon, path in files.items():
    if os.path.exists(path):
        output = f'syngraph_input/{taxon}.tsv'
        count = format_busco_file(path, output)
        print(f"✓ {taxon}: {count} genes")
    else:
        print(f"✗ {taxon}: file not found at {path}")
