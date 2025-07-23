"""
BUSCO processing module for the Genome Inversion Analyzer
Handles parsing, filtering, and sequence extraction from BUSCO results
"""

import pandas as pd
import logging
from Bio import SeqIO
from Bio.Seq import Seq

logger = logging.getLogger(__name__)


def enhanced_parse_busco_table(busco_path, config):
    """Enhanced BUSCO table parsing that correctly handles negative strand genes"""
    logger.info(f"Parsing BUSCO table: {busco_path}")
    
    with open(busco_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]
    
    busco_data = []
    parsing_errors = []
    duplicate_entries = []
    
    seen_entries = set()
    
    for line_num, line in enumerate(lines, 1):
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                try:
                    entry_key = (parts[0], parts[2], parts[3], parts[4])  # busco_id, sequence, start, end
                    
                    # Check for duplicates
                    if config.get('enable_duplicate_handling', False):
                        if entry_key in seen_entries:
                            duplicate_entries.append(line_num)
                            continue
                        seen_entries.add(entry_key)
                    
                    # Parse coordinates - handle both positive and negative strand
                    start_str = parts[3]
                    end_str = parts[4]
                    
                    if start_str != 'N/A' and end_str != 'N/A':
                        coord1 = int(start_str)
                        coord2 = int(end_str)
                        
                        # For genomic coordinates, always use min as start, max as end
                        gene_start = min(coord1, coord2)
                        gene_end = max(coord1, coord2)
                        
                        # Determine actual strand if not provided correctly
                        strand = parts[5] if len(parts) > 5 else '+'
                        if coord1 > coord2 and strand == '+':
                            strand = '-'  # Correct strand based on coordinates
                        elif coord1 < coord2 and strand == '-':
                            strand = '+'  # Correct strand based on coordinates
                        
                        entry = {
                            'busco_id': parts[0],
                            'status': parts[1],
                            'sequence': parts[2],
                            'gene_start': gene_start,
                            'gene_end': gene_end,
                            'strand': strand,
                            'score': float(parts[6]) if len(parts) > 6 and parts[6] != 'N/A' else None,
                            'length': int(parts[7]) if len(parts) > 7 and parts[7] != 'N/A' else (gene_end - gene_start),
                            'line_number': line_num,
                            'original_start': coord1,  # Keep original for debugging
                            'original_end': coord2
                        }
                        
                        busco_data.append(entry)
                    else:
                        parsing_errors.append(f"Line {line_num}: N/A coordinates")
                        
                except (ValueError, IndexError) as e:
                    parsing_errors.append(f"Line {line_num}: {e}")
            else:
                parsing_errors.append(f"Line {line_num}: Too few columns ({len(parts)})")
    
    busco_df = pd.DataFrame(busco_data)
    
    # Detect and handle paralogs if enabled
    if config.get('enable_paralog_detection', False):
        busco_df = detect_and_annotate_paralogs(busco_df, config)
    
    # Report parsing issues only if debug enabled
    total_lines = len(lines)
    success_rate = len(busco_data) / total_lines * 100 if total_lines > 0 else 0
    
    logger.info(f"  Parsing results:")
    logger.info(f"    Total data lines: {total_lines}")
    logger.info(f"    Successfully parsed: {len(busco_data)} ({success_rate:.1f}%)")
    logger.info(f"    Parsing errors: {len(parsing_errors)}")
    
    if len(parsing_errors) > 0 and config.get('enable_debug_output', False):
        logger.info(f"    First few errors:")
        for error in parsing_errors[:3]:
            logger.info(f"      {error}")
    
    if duplicate_entries and config.get('enable_debug_output', False):
        logger.warning(f"  {len(duplicate_entries)} duplicate entries removed")
    
    logger.info(f"  Found {len(busco_df)} valid BUSCO entries")
    return busco_df


def detect_and_annotate_paralogs(busco_df, config):
    """Enhanced paralog detection with detailed analysis"""
    logger.info("  Detecting and analyzing paralogs...")
    
    # Count occurrences of each BUSCO ID
    busco_counts = busco_df['busco_id'].value_counts()
    paralogs = busco_counts[busco_counts > 1]
    
    # Add paralog annotations
    busco_df['is_paralog'] = busco_df['busco_id'].isin(paralogs.index)
    busco_df['paralog_count'] = busco_df['busco_id'].map(busco_counts)
    
    # Initialize rank and cluster columns
    busco_df['paralog_rank'] = 1
    busco_df['paralog_cluster'] = busco_df['busco_id'] + '_singleton'
    
    # Process each unique BUSCO ID that has paralogs
    for busco_id in paralogs.index:
        paralog_group = busco_df[busco_df['busco_id'] == busco_id].copy()
        
        # Multi-criteria ranking
        ranking_criteria = []
        
        # Primary: score (if available)
        if 'score' in paralog_group.columns and paralog_group['score'].notna().any():
            ranking_criteria.append(paralog_group['score'].fillna(0))
        
        # Secondary: length (if available)
        if 'length' in paralog_group.columns and paralog_group['length'].notna().any():
            ranking_criteria.append(paralog_group['length'].fillna(0))
        
        # Tertiary: genomic position (earlier positions get higher rank)
        if paralog_group['gene_start'].notna().any():
            ranking_criteria.append(-paralog_group['gene_start'].fillna(float('inf')))
        
        if ranking_criteria:
            # Combine criteria with weights
            weights = [0.5, 0.3, 0.2][:len(ranking_criteria)]
            combined_score = sum(w * criteria for w, criteria in zip(weights, ranking_criteria))
            ranks = combined_score.rank(ascending=False, method='first')
        else:
            # If no ranking criteria available, rank by dataframe order
            ranks = pd.Series(range(1, len(paralog_group) + 1), index=paralog_group.index)
        
        # Update the main dataframe with ranks for this BUSCO ID
        busco_df.loc[paralog_group.index, 'paralog_rank'] = ranks
        busco_df.loc[paralog_group.index, 'paralog_cluster'] = f"{busco_id}_cluster"
    
    logger.info(f"    Found {len(paralogs)} BUSCO IDs with paralogs")
    logger.info(f"    Total paralogous genes: {len(busco_df[busco_df['is_paralog']])}")
    logger.info(f"    Max paralogs for single BUSCO: {busco_counts.max()}")
    
    return busco_df


def enhanced_filter_busco_genes(busco_df, config, quality_info=None):
    """Enhanced BUSCO filtering with adaptive parameters and detailed reporting"""
    logger.info("Enhanced BUSCO gene filtering...")
    
    initial_count = len(busco_df)
    filtering_stats = {'initial': initial_count}
    
    # Get adaptive parameters
    if quality_info and config.get('enable_adaptive_thresholds', False):
        filter_params = quality_info['adjustments']
    else:
        filter_params = {
            'min_busco_length': config['base_min_busco_length'],
            'similarity_threshold': config['base_similarity_threshold']
        }
    
    # Filter by status
    status_filter = config.get('busco_status_filter', ['Complete'])
    filtered_df = busco_df[busco_df['status'].isin(status_filter)].copy()
    filtering_stats['status_filter'] = len(filtered_df)
    
    # Filter by length if available and enabled
    if 'length' in filtered_df.columns and filter_params.get('min_busco_length'):
        length_before = len(filtered_df)
        min_length = filter_params['min_busco_length']
        filtered_df = filtered_df[filtered_df['length'] >= min_length]
        filtering_stats['length_filter'] = len(filtered_df)
        logger.info(f"  Length filter (>={min_length}bp): {length_before} -> {len(filtered_df)}")
    
    # Remove genes with missing coordinates
    coord_before = len(filtered_df)
    filtered_df = filtered_df.dropna(subset=['gene_start', 'gene_end'])
    filtering_stats['coordinate_filter'] = len(filtered_df)
    
    # Paralog handling
    if config.get('enable_paralog_detection', False):
        if config.get('enable_paralog_ortholog_mapping', False):
            # Keep all paralogs for downstream analysis
            pass
        else:
            # Keep only best paralog per BUSCO ID
            paralog_before = len(filtered_df)
            filtered_df = filtered_df.sort_values('paralog_rank').groupby('busco_id').first().reset_index()
            filtering_stats['paralog_filter'] = len(filtered_df)
            logger.info(f"  Paralog filter (best per BUSCO): {paralog_before} -> {len(filtered_df)}")
    
    final_count = len(filtered_df)
    exclusion_rate = (initial_count - final_count) / initial_count if initial_count > 0 else 0
    
    # Warning about excessive exclusions
    if config.get('enable_exclusion_warnings', False):
        if exclusion_rate > 0.7:
            logger.warning(f"Very high exclusion rate: {exclusion_rate:.1%} of BUSCOs excluded")
            logger.warning("Consider relaxing filtering parameters or checking assembly quality")
        elif exclusion_rate > 0.5:
            logger.warning(f"High exclusion rate: {exclusion_rate:.1%} of BUSCOs excluded")
    
    logger.info(f"  Filtering summary: {initial_count} -> {final_count} ({exclusion_rate:.1%} excluded)")
    
    # Add filtering statistics to dataframe
    if config.get('enable_debug_output', False):
        filtered_df['filtering_stats'] = str(filtering_stats)
    
    return filtered_df


def validate_gene_boundaries(gene_info, genome_seqs, config):
    """Enhanced gene boundary validation with comprehensive checks"""
    seq_id = gene_info['sequence']
    start = int(gene_info['gene_start']) - 1  # Convert to 0-based
    end = int(gene_info['gene_end'])
    strand = gene_info['strand']
    busco_id = gene_info['busco_id']
    
    validation_info = {
        'valid': False,
        'warnings': [],
        'errors': []
    }
    
    # Check if sequence exists in genome
    if seq_id not in genome_seqs:
        validation_info['errors'].append(f"Sequence {seq_id} not found in genome")
        return None, validation_info
    
    full_seq = genome_seqs[seq_id]
    chromosome_length = len(full_seq)
    
    # Check chromosome bounds
    if config.get('enable_chromosome_bounds_check', False):
        if start < 0:
            validation_info['warnings'].append(f"Gene start < 0, adjusted to 0")
            start = 0
        
        if end > chromosome_length:
            validation_info['warnings'].append(f"Gene end > chromosome length, adjusted")
            end = chromosome_length
        
        if start >= chromosome_length:
            validation_info['errors'].append(f"Gene start >= chromosome length")
            return None, validation_info
    
    # Validate coordinate order
    if start >= end:
        validation_info['errors'].append(f"Invalid coordinates: start >= end")
        return None, validation_info
    
    # Extract sequence
    gene_seq = full_seq[start:end]
    original_length = len(gene_seq)
    
    # Strand validation and processing
    if config.get('enable_strand_validation', False):
        if strand == '-':
            gene_seq = str(Seq(gene_seq).reverse_complement())
        
        # Check for valid translation if enabled
        if config.get('enable_translation_check', False):
            validation_info.update(validate_translation(gene_seq, busco_id))
    else:
        # Simple strand processing
        if strand == '-':
            gene_seq = str(Seq(gene_seq).reverse_complement())
    
    # Final validation
    if len(gene_seq) < 10:  # Minimum reasonable gene length
        validation_info['warnings'].append(f"Very short gene sequence ({len(gene_seq)} bp)")
    
    validation_info['valid'] = len(validation_info['errors']) == 0
    
    extracted_info = {
        'sequence': gene_seq,
        'validated_start': start + 1,  # Convert back to 1-based
        'validated_end': end,
        'validated_strand': strand,
        'length': len(gene_seq),
        'original_length': original_length,
        'chromosome_length': chromosome_length
    }
    
    return extracted_info, validation_info


def validate_translation(gene_seq, busco_id):
    """Validate gene translation for CDS quality"""
    validation_info = {'translation_warnings': [], 'translation_errors': []}
    
    try:
        # Check if sequence length is multiple of 3
        if len(gene_seq) % 3 != 0:
            validation_info['translation_warnings'].append("Sequence length not multiple of 3")
        
        # Attempt translation
        protein_seq = str(Seq(gene_seq).translate())
        
        # Check for premature stop codons
        internal_stops = protein_seq[:-1].count('*')
        if internal_stops > 0:
            validation_info['translation_warnings'].append(f"{internal_stops} internal stop codons")
        
        # Check for start codon (if long enough)
        if len(gene_seq) >= 3:
            start_codon = gene_seq[:3].upper()
            if start_codon not in ['ATG', 'GTG', 'TTG']:
                validation_info['translation_warnings'].append(f"Non-standard start codon: {start_codon}")
        
        # Check for reasonable protein length
        if len(protein_seq) < 50:
            validation_info['translation_warnings'].append(f"Short protein sequence ({len(protein_seq)} aa)")
        
    except Exception as e:
        validation_info['translation_errors'].append(f"Translation failed: {e}")
    
    return validation_info


def extract_enhanced_busco_sequences(busco_df, fasta_path, config):
    """Extract BUSCO sequences with comprehensive validation and error handling"""
    logger.info(f"Extracting BUSCO sequences from {fasta_path} with enhanced validation...")
    
    # Load genome sequences
    genome_seqs = {}
    sequence_lengths = {}
    
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            genome_seqs[record.id] = str(record.seq)
            sequence_lengths[record.id] = len(record.seq)
    except Exception as e:
        logger.error(f"Failed to load genome sequences: {e}")
        return pd.DataFrame()
    
    logger.info(f"  Loaded {len(genome_seqs)} sequences from genome")
    
    busco_with_seqs = []
    extraction_stats = {
        'total_attempted': len(busco_df),
        'successful': 0,
        'failed': 0,
        'warnings': 0
    }
    
    # Process each gene
    for idx, gene in busco_df.iterrows():
        try:
            extracted_info, validation_info = validate_gene_boundaries(gene, genome_seqs, config)
            
            if extracted_info and validation_info['valid']:
                # Create enhanced gene record
                gene_record = {
                    'busco_id': gene['busco_id'],
                    'sequence_id': gene['sequence'],
                    'gene_start': extracted_info['validated_start'],
                    'gene_end': extracted_info['validated_end'],
                    'strand': extracted_info['validated_strand'],
                    'gene_sequence': extracted_info['sequence'],
                    'gene_length': extracted_info['length'],
                    'status': gene['status'],
                    'is_paralog': gene.get('is_paralog', False),
                    'paralog_rank': gene.get('paralog_rank', 1),
                    'paralog_count': gene.get('paralog_count', 1),
                    'chromosome_length': extracted_info['chromosome_length']
                }
                
                # Add validation information if debug enabled
                if config.get('enable_debug_output', False):
                    gene_record.update({
                        'validation_warnings': ';'.join(validation_info.get('warnings', [])),
                        'translation_warnings': ';'.join(validation_info.get('translation_warnings', [])),
                        'extraction_method': 'enhanced'
                    })
                
                busco_with_seqs.append(gene_record)
                extraction_stats['successful'] += 1
                
                if validation_info.get('warnings'):
                    extraction_stats['warnings'] += 1
                    
            else:
                extraction_stats['failed'] += 1
                if config.get('enable_debug_output', False) and extraction_stats['failed'] <= 10:
                    errors = '; '.join(validation_info.get('errors', ['Unknown error']))
                    logger.warning(f"  Failed to extract {gene['busco_id']}: {errors}")
                    
        except Exception as e:
            extraction_stats['failed'] += 1
            if config.get('enable_debug_output', False) and extraction_stats['failed'] <= 10:
                logger.warning(f"  Extraction error for {gene['busco_id']}: {e}")
    
    # Report extraction statistics
    success_rate = extraction_stats['successful'] / extraction_stats['total_attempted'] * 100
    logger.info(f"  Extraction completed: {extraction_stats['successful']}/{extraction_stats['total_attempted']} ({success_rate:.1f}% success)")
    
    if extraction_stats['warnings'] > 0:
        logger.info(f"  {extraction_stats['warnings']} genes extracted with warnings")
    
    if extraction_stats['failed'] > 10 and config.get('enable_debug_output', False):
        logger.warning(f"  ... and {extraction_stats['failed'] - 10} more extraction failures")
    
    busco_seq_df = pd.DataFrame(busco_with_seqs)
    
    # Add extraction statistics to config for downstream use
    if config.get('enable_debug_output', False):
        busco_seq_df.attrs['extraction_stats'] = extraction_stats
    
    return busco_seq_df