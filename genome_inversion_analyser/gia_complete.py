#!/usr/bin/env python3
import sys, os
# add parent directory to Python path so the package root is on sys.path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# change working directory to project root so relative file paths resolve
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from genome_inversion_analyser import analyze_genomes

if __name__ == "__main__":
    results = analyze_genomes(
        first_fasta='../GCA_910594885.2_idBibMarc1.2_genomic.fna',
        second_fasta='../GCA_958336335.1_idDilFebr1.1_genomic.fna',
        first_busco='../Bibio_marci/full_table.tsv',
        second_busco='../Dilophus_febrilis/full_table.tsv',
        output_dir='../v4/complete_results'
    )
    print("Done!", results)