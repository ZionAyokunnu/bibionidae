# R Script: Full Genome Inversion Pipeline with Dotplot, Simulation, Detection, and Output
################################################################################
# Load required libraries
library(dplyr)                  # For data manipulation
library(tibble)                 # For tibble data frames
library(tidyr)                  # For data reshaping
library(ggplot2)                # For plotting
library(Biostrings)             # For FASTA file handling
library(stringr)                # For string manipulation

################################################################################
# Parameters - Toggle these as needed
use_busco <- FALSE                      # TRUE to filter BUSCO markers
simulate_manual_inversion <- TRUE       # TRUE to simulate inversions
num_inversions <- 1                     # Number of random inversions
input_fasta_path <- "GCA_910594885.2_idBibMarc1.2_genomic.fna" # Reference FASTA
inverted_fasta_output <- "Bibio_marci_inverted.fasta"           # Output inverted FASTA
inversion_summary_output <- "inversion_summary.csv"            # Inversion summary output CSV

################################################################################
# Load and Prepare Genome Data
get_scaffold_sizes <- function(fasta_path) {
  seqs <- readDNAStringSet(fasta_path)
  tibble(
    chromosome = names(seqs),
    chromsome_size_b = as.numeric(width(seqs)),
    sequence = as.character(seqs)
  )
}

genome_data <- get_scaffold_sizes(input_fasta_path)

# Expand into gene-level dummy markers (1 marker per 1000 bp block)
genome_df <- genome_data %>%
  mutate(marker_count = floor(chromsome_size_b / 1000)) %>%
  uncount(marker_count, .id = "marker_index") %>%
  group_by(chromosome) %>%
  mutate(
    start = marker_index * 1000,
    end = start + 999,
    marker = paste0(chromosome, "_", marker_index),
    linear_start = cumsum(if_else(row_number() == 1, 0, lag(end - start + 1))),
    linear_end = linear_start + (end - start + 1)
  ) %>%
  ungroup()

################################################################################
# Manual Inversion Simulation (if enabled)
invert_segment <- function(genome_vec) {
  n <- length(genome_vec)
  i <- sample(1:(n - 1), 1)
  j <- sample((i + 1):n, 1)
  genome_vec[i:j] <- rev(genome_vec[i:j])
  return(list(new_genome = genome_vec, inversion = c(i, j)))
}

simulate_inversions <- function(genome_vec, num_inversions = 10) {
  genome_states <- list()
  inversions_log <- matrix(NA, nrow = num_inversions, ncol = 2)
  colnames(inversions_log) <- c("start_idx", "end_idx")
  
  for (k in 1:num_inversions) {
    result <- invert_segment(genome_vec)
    genome_vec <- result$new_genome
    inversions_log[k, ] <- result$inversion
    genome_states[[k]] <- genome_vec
  }
  return(list(final_genome = genome_vec, inversions = inversions_log, history = genome_states))
}

if (simulate_manual_inversion) {
  gene_vector <- genome_df$marker
  sim_result <- simulate_inversions(gene_vector, num_inversions)
  
  # Plot inversion
  df_plot <- data.frame(
    original_index = 1:length(gene_vector),
    final_index = match(gene_vector, sim_result$final_genome),
    gene = gene_vector
  )
  
  p <- ggplot(df_plot, aes(x = original_index, y = final_index)) +
    geom_point(color = "darkred", size = 1.5) +
    theme_minimal() +
    labs(title = "Manual Inversion Dotplot", x = "Original Index", y = "Inverted Index")
  ggsave("manual_inversion_dotplot.png", p)
  
  # Frequency analysis
  inverted_genes <- c()
  for (k in 1:nrow(sim_result$inversions)) {
    i <- sim_result$inversions[k, 1]
    j <- sim_result$inversions[k, 2]
    inverted_genes <- c(inverted_genes, gene_vector[i:j])
  }
  gene_freq <- as.data.frame(table(inverted_genes)) %>% arrange(desc(Freq))
  write.csv(gene_freq, "inversion_frequency.csv", row.names = FALSE)
  
  # Summary output
  inversion_df <- data.frame(
    inversion_id = 1:nrow(sim_result$inversions),
    start_idx = sim_result$inversions[, 1],
    end_idx = sim_result$inversions[, 2],
    num_genes = sim_result$inversions[, 2] - sim_result$inversions[, 1] + 1
  )
  
  inversion_df$start_gene <- gene_vector[inversion_df$start_idx]
  inversion_df$end_gene <- gene_vector[inversion_df$end_idx]
  inversion_df$size_mb <- inversion_df$num_genes * 1000 / 1e6
  
  write.csv(inversion_df, inversion_summary_output, row.names = FALSE)
}

################################################################################
# Save Final Inverted Genome to FASTA Format
if (simulate_manual_inversion) {
  inverted_df <- genome_df
  inverted_df$marker <- sim_result$final_genome
  save_order <- match(sim_result$final_genome, genome_df$marker)
  inverted_df <- inverted_df[save_order, ]
  
  # Join sequences from original data for saving
  chrom_fasta <- genome_data %>% select(chromosome, sequence)
  writeXStringSet(DNAStringSet(chrom_fasta$sequence), inverted_fasta_output)
}

################################################################################
# Compare Inverted vs Original - Count Inversions Automatically
compare_gene_orders <- function(original_vec, altered_vec) {
  differences <- which(original_vec != altered_vec)
  inv_count <- 0
  if (length(differences) == 0) return(0)
  i <- 1
  while (i < length(differences)) {
    j <- i + 1
    while (j <= length(differences) && differences[j] == differences[j-1] + 1) j <- j + 1
    inv_count <- inv_count + 1
    i <- j
  }
  return(inv_count)
}

if (simulate_manual_inversion) {
  inv_count <- compare_gene_orders(gene_vector, sim_result$final_genome)
  cat("Detected inversion blocks:", inv_count, "\n")
}
