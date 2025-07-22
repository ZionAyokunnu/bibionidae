from Bio import SeqIO

# Check markers around the inversion (14593-16233)
def check_inversion_region():
    # Read both files
    original_seqs = list(SeqIO.parse('GCA_910594885.2_idBibMarc1.2_genomic.fna', 'fasta'))
    inverted_seqs = list(SeqIO.parse('Bibio_marci_inverted.fasta', 'fasta'))
    
    orig_chr1 = str(original_seqs[0].seq)
    inv_chr1 = str(inverted_seqs[0].seq)
    
    print(f"Original chr1 length: {len(orig_chr1):,}")
    print(f"Inverted chr1 length: {len(inv_chr1):,}")
    
    # Check marker regions around the inversion
    marker_size = 5000
    
    print("\n=== CHECKING MARKERS AROUND INVERSION (14593-16233) ===")
    
    # Check marker 14590 (before inversion)
    start_bp = (14590 - 1) * marker_size
    end_bp = 14590 * marker_size
    orig_seq = orig_chr1[start_bp:end_bp]
    inv_seq = inv_chr1[start_bp:end_bp]
    print(f"Marker 14590 (before): {'SAME' if orig_seq == inv_seq else 'DIFFERENT'}")
    
    # Check marker 14593 (start of inversion)
    start_bp = (14593 - 1) * marker_size
    end_bp = 14593 * marker_size
    orig_seq = orig_chr1[start_bp:end_bp]
    inv_seq = inv_chr1[start_bp:end_bp]
    print(f"Marker 14593 (inversion start): {'SAME' if orig_seq == inv_seq else 'DIFFERENT'}")
    
    # Check marker 16233 (end of inversion)
    start_bp = (16233 - 1) * marker_size
    end_bp = 16233 * marker_size
    orig_seq = orig_chr1[start_bp:end_bp] if end_bp <= len(orig_chr1) else orig_chr1[start_bp:]
    inv_seq = inv_chr1[start_bp:end_bp] if end_bp <= len(inv_chr1) else inv_chr1[start_bp:]
    print(f"Marker 16233 (inversion end): {'SAME' if orig_seq == inv_seq else 'DIFFERENT'}")
    
    # Check marker 16240 (after inversion)
    start_bp = (16240 - 1) * marker_size
    end_bp = 16240 * marker_size
    orig_seq = orig_chr1[start_bp:end_bp] if end_bp <= len(orig_chr1) else orig_chr1[start_bp:]
    inv_seq = inv_chr1[start_bp:end_bp] if end_bp <= len(inv_chr1) else inv_chr1[start_bp:]
    print(f"Marker 16240 (after): {'SAME' if orig_seq == inv_seq else 'DIFFERENT'}")

check_inversion_region()