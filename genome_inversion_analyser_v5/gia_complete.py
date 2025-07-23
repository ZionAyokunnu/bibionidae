#!/usr/bin/env python3
from genome_inversion_analyser_v5 import analyze_genomes

if __name__ == "__main__":
    results = analyze_genomes(
        first_fasta='GCA_910594885.2_idBibMarc1.2_genomic.fna',
        second_fasta='GCA_958336335.1_idDilFebr1.1_genomic.fna',
        first_busco='Bibio_marci/full_table.tsv',
        second_busco='Dilophus_febrilis/full_table.tsv',
        output_dir='v4/complete_results'
    )
    print("Done!", results)