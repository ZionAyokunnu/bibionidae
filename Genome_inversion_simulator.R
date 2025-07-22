################################################################################
# MODULE 1: GENOME INVERSION SIMULATOR
# Purpose: Load FASTA, simulate inversions, create dotplot, save inverted FASTA
################################################################################

# Load required libraries
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(stringr)

################################################################################
# CONFIGURATION PARAMETERS
################################################################################
CONFIG <- list(
  # Input/Output paths
  input_fasta_path = "GCA_910594885.2_idBibMarc1.2_genomic.fna",
  inverted_fasta_output = "Bibio_marci_inverted.fasta",
  dotplot_output = "manual_inversion_dotplot.png",
  
  # Simulation parameters
  simulate_inversions = TRUE,
  num_inversions = 1,
  marker_size_bp = 5000,  # Size of each marker block in base pairs
  
  # Optional features
  use_busco = FALSE,
  save_intermediate_states = TRUE
)

################################################################################
# CORE FUNCTIONS
################################################################################

#' Load and process FASTA file into genome markers
#' @param fasta_path Path to input FASTA file
#' @param marker_size Size of each marker block in base pairs
#' @return tibble with genome markers and positions
load_and_process_genome <- function(fasta_path, marker_size = 1000) {
  cat("Loading genome from:", fasta_path, "\n")
  
  # Read FASTA sequences
  seqs <- readDNAStringSet(fasta_path)
  
  # Create scaffold/chromosome information
  genome_data <- tibble(
    chromosome = names(seqs),
    chromosome_size_bp = as.numeric(width(seqs)),
    sequence = as.character(seqs)
  )
  
  cat("Loaded", nrow(genome_data), "scaffolds/chromosomes\n")
  cat("Total genome size:", sum(genome_data$chromosome_size_bp), "bp\n")
  
  # Create markers (1 marker per specified bp block)
  genome_df <- genome_data %>%
    mutate(marker_count = floor(chromosome_size_bp / marker_size)) %>%
    filter(marker_count > 0) %>%  # Remove very small scaffolds
    uncount(marker_count, .id = "marker_index") %>%
    group_by(chromosome) %>%
    mutate(
      start_bp = (marker_index - 1) * marker_size + 1,
      end_bp = marker_index * marker_size,
      marker_id = paste0(chromosome, "_", sprintf("%06d", marker_index)),
      # Create linear genome position for plotting
      scaffold_start = min(start_bp),
      scaffold_end = max(end_bp)
    ) %>%
    ungroup() %>%
    arrange(chromosome, marker_index)
  
  # Add cumulative positions across all chromosomes
  genome_df <- genome_df %>%
    group_by(chromosome) %>%
    mutate(
      chr_cumulative_start = cumsum(c(0, rep(marker_size, n()-1))),
      chr_cumulative_end = chr_cumulative_start + marker_size
    ) %>%
    ungroup()
  
  # Add global linear positions
  chr_offsets <- genome_df %>%
    group_by(chromosome) %>%
    summarise(
      chr_max = max(chr_cumulative_end),
      .groups = "drop"
    ) %>%
    mutate(
      chr_offset = cumsum(c(0, head(chr_max, -1)))
    )
  
  genome_df <- genome_df %>%
    left_join(chr_offsets, by = "chromosome") %>%
    mutate(
      linear_start = chr_cumulative_start + chr_offset,
      linear_end = chr_cumulative_end + chr_offset
    ) %>%
    select(-chr_cumulative_start, -chr_cumulative_end, -chr_offset)
  
  cat("Created", nrow(genome_df), "markers\n")
  
  return(list(
    genome_data = genome_data,
    genome_df = genome_df
  ))
}

#' Simulate a single inversion event
#' @param genome_vec Vector of gene/marker identifiers
#' @return List with new genome and inversion coordinates
simulate_single_inversion <- function(genome_vec) {
  n <- length(genome_vec)
  if (n < 2) {
    stop("Genome too small for inversion")
  }
  
  # Choose random inversion boundaries
  start_idx <- sample(1:(n-1), 1)
  end_idx <- sample((start_idx + 1):n, 1)
  
  # Perform inversion
  new_genome <- genome_vec
  new_genome[start_idx:end_idx] <- rev(genome_vec[start_idx:end_idx])
  
  return(list(
    new_genome = new_genome,
    inversion_start = start_idx,
    inversion_end = end_idx,
    inversion_size = end_idx - start_idx + 1
  ))
}

#' Simulate multiple inversion events
#' @param genome_vec Vector of gene/marker identifiers
#' @param num_inversions Number of inversions to simulate
#' @return List with final genome, inversion log, and history
simulate_multiple_inversions <- function(genome_vec, num_inversions = 1) {
  cat("Simulating", num_inversions, "inversion(s)\n")
  
  current_genome <- genome_vec
  inversions_log <- tibble(
    inversion_id = integer(),
    start_idx = integer(),
    end_idx = integer(),
    size = integer(),
    original_start_marker = character(),
    original_end_marker = character()
  )
  
  genome_history <- list()
  genome_history[[1]] <- current_genome  # Store original
  
  for (i in 1:num_inversions) {
    result <- simulate_single_inversion(current_genome)
    
    # Log this inversion
    inversions_log <- bind_rows(inversions_log, tibble(
      inversion_id = i,
      start_idx = result$inversion_start,
      end_idx = result$inversion_end,
      size = result$inversion_size,
      original_start_marker = genome_vec[result$inversion_start],
      original_end_marker = genome_vec[result$inversion_end]
    ))
    
    current_genome <- result$new_genome
    genome_history[[i + 1]] <- current_genome
    
    cat("  Inversion", i, ":", result$inversion_start, "to", result$inversion_end, 
        "(", result$inversion_size, "markers )\n")
  }
  
  return(list(
    original_genome = genome_vec,
    final_genome = current_genome,
    inversions_log = inversions_log,
    genome_history = genome_history
  ))
}

#' Create dotplot comparing original and inverted genomes
#' @param original_genome Original gene order
#' @param inverted_genome Inverted gene order
#' @param output_path Path to save plot
#' @param title Plot title
create_inversion_dotplot <- function(original_genome, inverted_genome, 
                                     output_path, title = "Genome Inversion Dotplot") {
  cat("Creating dotplot...\n")
  
  # Create plotting data
  plot_data <- tibble(
    original_index = 1:length(original_genome),
    inverted_index = match(original_genome, inverted_genome),
    marker_id = original_genome
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = original_index, y = inverted_index)) +
    geom_point(color = "darkred", size = 1.2, alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, color = "gray70", linetype = "dashed") +
    theme_minimal() +
    labs(
      title = title,
      subtitle = paste("Total markers:", length(original_genome)),
      x = "Original Genome Position (Index)",
      y = "Inverted Genome Position (Index)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  # Save plot
  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat("Dotplot saved to:", output_path, "\n")
  
  return(p)
}

#' Save inverted genome to FASTA format with efficient sequence reconstruction
#' @param genome_data Original genome data with sequences
#' @param original_genome Original marker order
#' @param inverted_genome New order of markers after inversions
#' @param genome_df Genome dataframe with marker information
#' @param output_path Path to save inverted FASTA
#' @param inversions_log Log of inversions performed
save_inverted_fasta <- function(genome_data, original_genome, inverted_genome, genome_df, output_path, inversions_log) {
  cat("Saving inverted genome to FASTA with efficient sequence reconstruction...\n")
  
  cat("Using efficient inversion-based reconstruction...\n")
  cat("Inversion details:\n")
  print(inversions_log)
  
  # Start with original sequences
  reconstructed_sequences <- genome_data$sequence
  names(reconstructed_sequences) <- genome_data$chromosome
  
  # Apply each inversion efficiently
  for (inv_idx in 1:nrow(inversions_log)) {
    inversion <- inversions_log[inv_idx, ]
    
    cat("Applying inversion", inv_idx, ":\n")
    cat("  Start marker index:", inversion$start_idx, "\n")
    cat("  End marker index:", inversion$end_idx, "\n")
    cat("  Size:", inversion$size, "markers\n")
    
    # Get the chromosomes and positions for start and end markers
    start_marker_info <- genome_df %>% 
      filter(marker_id == inversion$original_start_marker)
    
    end_marker_info <- genome_df %>% 
      filter(marker_id == inversion$original_end_marker)
    
    # Take first match if multiple
    if (nrow(start_marker_info) > 0) {
      start_marker_info <- start_marker_info[1, ]
    }
    if (nrow(end_marker_info) > 0) {
      end_marker_info <- end_marker_info[1, ]
    }
    
    if (nrow(start_marker_info) == 0 || nrow(end_marker_info) == 0) {
      cat("  Warning: Could not find marker positions, skipping inversion\n")
      next
    }
    
    cat("  Start marker:", start_marker_info$marker_id, "\n")
    cat("  End marker:", end_marker_info$marker_id, "\n")
    
    # Check if inversion spans multiple chromosomes
    if (start_marker_info$chromosome != end_marker_info$chromosome) {
      cat("  Warning: Inversion spans multiple chromosomes - using complex reconstruction\n")
      
      # For multi-chromosome inversions, use the marker-based approach
      # Get affected chromosomes
      affected_chroms <- unique(c(start_marker_info$chromosome, end_marker_info$chromosome))
      
      for (chrom in affected_chroms) {
        # Get markers for this chromosome in original order
        chrom_markers_orig <- genome_df %>%
          filter(chromosome == chrom) %>%
          arrange(start_bp) %>%
          mutate(orig_order = row_number())
        
        # Get markers in inverted order
        chrom_markers_inv <- tibble(
          new_order = 1:length(inverted_genome),
          marker_id = inverted_genome
        ) %>%
          inner_join(chrom_markers_orig, by = "marker_id") %>%
          arrange(new_order)
        
        if (nrow(chrom_markers_inv) > 0) {
          # Reconstruct this chromosome's sequence
          original_seq <- reconstructed_sequences[[chrom]]
          new_seq <- ""
          
          for (i in 1:nrow(chrom_markers_inv)) {
            marker <- chrom_markers_inv[i, ]
            start_bp <- max(1, marker$start_bp)
            end_bp <- min(marker$end_bp, nchar(original_seq))
            
            if (start_bp <= end_bp && start_bp <= nchar(original_seq)) {
              segment <- substr(original_seq, start_bp, end_bp)
              new_seq <- paste0(new_seq, segment)
            }
          }
          
          if (nchar(new_seq) > 0) {
            reconstructed_sequences[[chrom]] <- new_seq
            cat("    Reconstructed", chrom, ":", nchar(new_seq), "bp\n")
          }
        }
      }
      
    } else {
      # Single chromosome inversion - use efficient method
      chrom <- start_marker_info$chromosome
      original_seq <- reconstructed_sequences[[chrom]]
      
      # Convert marker indices to base pair positions
      start_bp <- start_marker_info$start_bp
      end_bp <- end_marker_info$end_bp
      
      # Ensure positions are valid
      start_bp <- max(1, min(start_bp, nchar(original_seq)))
      end_bp <- max(start_bp, min(end_bp, nchar(original_seq)))
      
      cat("  Inverting", chrom, "from bp", start_bp, "to", end_bp, "\n")
      
      if (start_bp < end_bp) {
        # Extract segments
        before_segment <- if(start_bp > 1) substr(original_seq, 1, start_bp - 1) else ""
        inversion_segment <- substr(original_seq, start_bp, end_bp)
        after_segment <- if(end_bp < nchar(original_seq)) substr(original_seq, end_bp + 1, nchar(original_seq)) else ""
        
        # Reverse the inversion segment
        reversed_segment <- paste(rev(strsplit(inversion_segment, "")[[1]]), collapse = "")
        
        # Reconstruct the sequence
        new_seq <- paste0(before_segment, reversed_segment, after_segment)
        reconstructed_sequences[[chrom]] <- new_seq
        
        cat("    Inverted", nchar(inversion_segment), "bp in", chrom, "\n")
        cat("    New sequence length:", nchar(new_seq), "bp\n")
      }
    }
  }
  
  # Save reconstructed sequences
  sequences <- DNAStringSet(reconstructed_sequences)
  writeXStringSet(sequences, output_path)
  
  # Save marker order file for verification
  marker_order_file <- gsub("\\.fasta$", "_marker_order.csv", output_path)
  
  # Map original positions to inverted positions
  original_to_inverted_pos <- match(original_genome, inverted_genome)
  
  inverted_marker_df <- tibble(
    original_position = 1:length(original_genome),
    marker_id = original_genome,
    inverted_position = original_to_inverted_pos
  ) %>%
    left_join(
      genome_df %>% select(marker_id, chromosome, start_bp, end_bp),
      by = "marker_id"
    ) %>%
    arrange(original_position)
  
  # Verify the inversion is captured
  inversions_detected <- sum(diff(inverted_marker_df$inverted_position) < 0)
  cat("Detected", inversions_detected, "decreasing positions in marker order\n")
  
  write_csv(inverted_marker_df, marker_order_file)
  cat("Marker order saved to:", marker_order_file, "\n")
  
  # Save inversion log
  inversion_log_file <- gsub("\\.fasta$", "_inversion_log.csv", output_path)
  write_csv(inversions_log, inversion_log_file)
  cat("Inversion log saved to:", inversion_log_file, "\n")
  
  cat("Inverted FASTA with efficiently reconstructed sequences saved to:", output_path, "\n")
  
  # Calculate and report size differences
  original_total_size <- sum(nchar(genome_data$sequence))
  reconstructed_total_size <- sum(nchar(reconstructed_sequences))
  
  cat("Original genome size:", original_total_size, "bp\n")
  cat("Reconstructed genome size:", reconstructed_total_size, "bp\n")
  
  if (abs(original_total_size - reconstructed_total_size) < 1000) {
    cat("✓ Sequence reconstruction successful (size difference < 1kb)\n")
  } else {
    cat("⚠ Warning: Significant size difference in reconstruction\n")
  }
  
  return(list(
    marker_order_file = marker_order_file,
    inversion_log_file = inversion_log_file,
    inversions_in_order = inversions_detected,
    original_size = original_total_size,
    reconstructed_size = reconstructed_total_size
  ))
}

################################################################################
# MAIN EXECUTION FUNCTION
################################################################################

#' Main function to run the genome inversion simulation
#' @param config Configuration list with parameters
run_genome_inversion_simulation <- function(config = CONFIG) {
  cat("=== GENOME INVERSION SIMULATOR ===\n")
  cat("Starting simulation with", config$num_inversions, "inversion(s)\n\n")
  
  # Step 1: Load and process genome
  genome_result <- load_and_process_genome(
    config$input_fasta_path, 
    config$marker_size_bp
  )
  
  genome_data <- genome_result$genome_data
  genome_df <- genome_result$genome_df
  
  # Step 2: Simulate inversions (if enabled)
  if (config$simulate_inversions) {
    marker_vector <- genome_df$marker_id
    
    simulation_result <- simulate_multiple_inversions(
      marker_vector, 
      config$num_inversions
    )
    
    # Step 3: Create dotplot
    dotplot <- create_inversion_dotplot(
      simulation_result$original_genome,
      simulation_result$final_genome,
      config$dotplot_output,
      paste("Manual Inversion Dotplot -", config$num_inversions, "Inversion(s)")
    )
    
    # Step 4: Save inverted FASTA
    save_result <- save_inverted_fasta(
      genome_data,
      simulation_result$original_genome,
      simulation_result$final_genome,
      genome_df,
      config$inverted_fasta_output,
      simulation_result$inversions_log
    )
    
    cat("Marker order verification: found", save_result$inversions_in_order, "inversions in saved order\n")
    
    # Step 5: Return results for potential further analysis
    cat("\n=== SIMULATION COMPLETED ===\n")
    cat("Files created:\n")
    cat("  - Dotplot:", config$dotplot_output, "\n")
    cat("  - Inverted FASTA:", config$inverted_fasta_output, "\n")
    
    return(list(
      config = config,
      genome_data = genome_data,
      genome_df = genome_df,
      simulation_result = simulation_result,
      dotplot = dotplot
    ))
  } else {
    cat("Simulation disabled. Only genome loading performed.\n")
    return(list(
      config = config,
      genome_data = genome_data,
      genome_df = genome_df
    ))
  }
}

################################################################################
# RUN THE SIMULATION
################################################################################

# Execute the simulation
if (interactive() || !exists("SKIP_EXECUTION")) {
  results <- run_genome_inversion_simulation(CONFIG)
  
  # Print summary
  if (CONFIG$simulate_inversions) {
    cat("\nInversion Summary:\n")
    print(results$simulation_result$inversions_log)
  }
}