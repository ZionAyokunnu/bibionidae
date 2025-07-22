# from Bio import SeqIO

# # Read both files
# original_seqs = list(SeqIO.parse('GCA_910594885.2_idBibMarc1.2_genomic.fna', 'fasta'))
# inverted_seqs = list(SeqIO.parse('Bibio_marci_inverted.fasta', 'fasta'))

# print(f"Original file: {len(original_seqs)} sequences")
# print(f"Inverted file: {len(inverted_seqs)} sequences")

# # Check chromosome 1 specifically (where the inversion happened)
# orig_chr1 = str(original_seqs[0].seq)
# inv_chr1 = str(inverted_seqs[0].seq)

# print(f"Original chr1 length: {len(orig_chr1):,}")
# print(f"Inverted chr1 length: {len(inv_chr1):,}")
# print(f"Sequences identical: {orig_chr1 == inv_chr1}")

# # Check specific inversion region (bp 72,960,001 to 81,165,000)
# start_bp = 72960000  # 0-based
# end_bp = 81165000

# orig_segment = orig_chr1[start_bp:end_bp]
# inv_segment = inv_chr1[start_bp:end_bp]

# print(f"Inversion region length: {len(orig_segment):,}")
# print(f"Inversion region identical: {orig_segment == inv_segment}")
# print(f"Inversion region is reversed: {orig_segment == inv_segment[::-1]}")


import pandas as pd
from Bio import SeqIO

def process_fasta_to_markers(fasta_path, marker_size=5000):
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    markers_list = []
    
    for seq_record in sequences:
        chromosome = seq_record.id
        chr_size = len(seq_record.seq)
        sequence = str(seq_record.seq)
        
        marker_count = chr_size // marker_size
        
        if marker_count > 0:
            for marker_idx in range(1, min(11, marker_count + 1)):  # Just first 10 markers
                start_bp = (marker_idx - 1) * marker_size + 1
                end_bp = marker_idx * marker_size
                marker_id = f"{chromosome}_{marker_idx:06d}"
                
                # Get actual sequence content
                marker_seq = sequence[start_bp-1:end_bp]  # 0-based indexing
                
                markers_list.append({
                    'marker_id': marker_id,
                    'start_bp': start_bp,
                    'end_bp': end_bp,
                    'sequence_start': marker_seq[:50] + "..." if len(marker_seq) > 50 else marker_seq
                })
                
    return pd.DataFrame(markers_list)

# Check both files
print("=== ORIGINAL MARKERS ===")
orig_markers = process_fasta_to_markers('GCA_910594885.2_idBibMarc1.2_genomic.fna')
print(orig_markers.head())

print("\n=== INVERTED MARKERS ===")
inv_markers = process_fasta_to_markers('Bibio_marci_inverted.fasta')
print(inv_markers.head())

print("\n=== COMPARISON ===")
for i in range(5):
    orig_seq = orig_markers.iloc[i]['sequence_start']
    inv_seq = inv_markers.iloc[i]['sequence_start']
    print(f"Marker {i+1}: {'SAME' if orig_seq == inv_seq else 'DIFFERENT'}")