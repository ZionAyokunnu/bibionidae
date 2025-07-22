#' Detect inversions between two genome arrangements (improved version)
#' @param comparison_df Comparison dataframe with original and inverted positions
#' @return List with inversion detection results
detect_inversions_improved <- function(comparison_df) {
  cat("Detecting inversions (improved algorithm)...\n")
  
  # Sort by original position to ensure proper order
  comparison_df <- comparison_df %>% arrange(original_position)
  
  # Initialize inversions tibble
  inversions <- tibble(
    inversion_id = integer(),
    start_pos = integer(),
    end_pos = integer(),
    size = integer(),
    start_marker = character(),
    end_marker = character(),
    type = character()
  )
  
  if (nrow(comparison_df) < 2) {
    return(list(inversions = inversions, comparison = comparison_df))
  }
  
  # Detect inversions by finding regions where inverted_position decreases
  i <- 1
  inversion_count <- 0
  
  while (i < nrow(comparison_df)) {
    # Check if this is the start of an inversion (decreasing sequence)
    if (comparison_df$inverted_position[i] > comparison_df$inverted_position[i + 1]) {
      inversion_count <- inversion_count + 1
      start_pos <- i
      
      # Find the end of this inversion
      j <- i + 1
      while (j < nrow(comparison_df) && 
             comparison_df$inverted_position[j] <= comparison_df$inverted_position[j - 1]) {
        j <- j + 1
      }
      end_pos <- j - 1
      
      # Record this inversion if it's large enough
      inversion_size <- end_pos - start_pos + 1
      if (inversion_size >= 2) {  # Minimum size filter
        inversions <- bind_rows(inversions, tibble(
          inversion_id = inversion_count,
          start_pos = start_pos,
          end_pos = end_pos,
          size = inversion_size,
          start_marker = comparison_df$marker_id[start_pos],
          end_marker = comparison_df$marker_id[end_pos],
          type = "inversion"
        ))
      }
      
      i <- j
    } else {
      i <- i + 1
    }
  }
  
  cat("  Detected", nrow(inversions), "inversion blocks\n")
  
  return(list(
    inversions = inversions,
    comparison = comparison_df
  ))
}

#' Count inversion events per marker (improved version)
#' @param inversions Inversion detection results
#' @param comparison_df Comparison dataframe
#' @return tibble with inversion frequencies
calculate_inversion_frequency_improved <- function(inversions, comparison_df) {
  cat("Calculating inversion frequencies (improved)...\n")
  
  if (nrow(inversions) == 0) {
    return(tibble(
      marker_id = character(),
      frequency = integer(),
      total_events = integer()
    ))
  }
  
  # Get all markers involved in inversions
  involved_markers <- c()
  for (i in 1:nrow(inversions)) {
    start_pos <- inversions$start_pos[i]
    end_pos <- inversions$end_pos[i]
    involved_markers <- c(involved_markers, 
                          comparison_df$marker_id[start_pos:end_pos])
  }
  
  # Count frequencies
  frequency_df <- tibble(marker_id = involved_markers) %>%
    count(marker_id, name = "frequency") %>%
    mutate(total_events = nrow(inversions)) %>%
    arrange(desc(frequency))
  
  cat("  ", nrow(frequency_df), "markers involved in inversions\n")
  
  return(frequency_df)
}

#' Create comprehensive inversion analysis (improved version)
#' @param inversions Detected inversions
#' @param markers_df Markers dataframe
#' @param config Analysis configuration
#' @return tibble with detailed analysis
create_detailed_analysis_improved <- function(inversions, markers_df, config) {
  cat("Creating detailed analysis (improved)...\n")
  
  if (nrow(inversions) == 0) {
    return(tibble(
      inversion_id = integer(),
      chromosome = character(),
      start_bp = integer(),
      end_bp = integer(),
      size_markers = integer(),
      size_mb = numeric(),
      num_genes = integer(),
      inversion_rate_per_mb = numeric(),
      type = character(),
      inversion_events = integer()
    ))
  }
  
  # Create mapping from positions to genomic coordinates
  position_to_coords <- markers_df %>%
    mutate(linear_position = row_number()) %>%
    select(linear_position, chromosome, start_bp, end_bp, marker_id)
  
  # Enhance inversions with genomic coordinates
  detailed_analysis <- inversions %>%
    left_join(
      position_to_coords %>% rename(start_pos = linear_position),
      by = "start_pos"
    ) %>%
    rename(start_chromosome = chromosome, start_bp_coord = start_bp, start_end_bp = end_bp) %>%
    left_join(
      position_to_coords %>% rename(end_pos = linear_position),
      by = "end_pos"
    ) %>%
    rename(end_chromosome = chromosome, end_bp_coord = start_bp, end_end_bp = end_bp) %>%
    mutate(
      chromosome = coalesce(start_chromosome, end_chromosome),
      start_bp = start_bp_coord,
      end_bp = end_end_bp,
      size_markers = size,
      size_mb = size * config$marker_size_bp / 1e6,
      num_genes = size,  # Assuming 1 gene per marker for simplicity
      inversion_rate_per_mb = ifelse(size_mb > 0, 1 / size_mb, 0),
      inversion_events = if_else(config$detect_inversion_events, 1, 0)
    ) %>%
    select(
      inversion_id, chromosome, start_bp, end_bp, 
      size_markers, size_mb, num_genes, inversion_rate_per_mb, 
      type, inversion_events
    )
  
  cat("  Detailed analysis for", nrow(detailed_analysis), "inversions\n")
  
  return(detailed_analysis)
}################################################################################
# MODULE 2: GENOME INVERSION ANALYZER
# Purpose: Compare original and inverted FASTA files, detect inversions, 
#          create analysis plots, and generate detailed CSV reports
################################################################################

# Load required libraries
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(stringr)
library(readr)

################################################################################
# CONFIGURATION PARAMETERS
################################################################################
ANALYSIS_CONFIG <- list(
  # Input files
  original_fasta_path = "GCA_910594885.2_idBibMarc1.2_genomic.fna",
  inverted_fasta_path = "Bibio_marci_inverted.fasta",
  
  # Output files
  inversion_summary_csv = "inversion_summary.csv",
  inversion_frequency_csv = "inversion_frequency.csv",
  detailed_analysis_csv = "detailed_inversion_analysis.csv",
  comparison_dotplot_png = "comparison_dotplot.png",
  inversion_heatmap_png = "inversion_heatmap.png",
  
  # Analysis parameters
  marker_size_bp = 5000,
  detect_inversion_events = TRUE,  # Boolean toggle for event detection
  min_inversion_size = 2,          # Minimum markers for inversion detection
  
  # Plotting parameters
  plot_width = 12,
  plot_height = 8,
  dpi = 300
)

################################################################################
# CORE ANALYSIS FUNCTIONS
################################################################################

#' Process FASTA file into comparable marker format
#' @param fasta_path Path to FASTA file
#' @param marker_size Size of each marker in base pairs
#' @return tibble with processed genome markers
process_fasta_to_markers <- function(fasta_path, marker_size = 1000) {
  cat("Processing FASTA file:", fasta_path, "\n")
  
  # Read sequences
  seqs <- readDNAStringSet(fasta_path)
  
  # Create marker tibble
  genome_data <- tibble(
    chromosome = names(seqs),
    chromosome_size_bp = as.numeric(width(seqs)),
    sequence = as.character(seqs)
  )
  
  # Generate markers
  markers_df <- genome_data %>%
    mutate(marker_count = floor(chromosome_size_bp / marker_size)) %>%
    filter(marker_count > 0) %>%
    uncount(marker_count, .id = "marker_index") %>%
    group_by(chromosome) %>%
    mutate(
      start_bp = (marker_index - 1) * marker_size + 1,
      end_bp = marker_index * marker_size,
      marker_id = paste0(chromosome, "_", sprintf("%06d", marker_index)),
      size_mb = marker_size / 1e6
    ) %>%
    ungroup() %>%
    arrange(chromosome, marker_index) %>%
    mutate(
      linear_position = row_number()
    )
  
  cat("  Generated", nrow(markers_df), "markers\n")
  
  return(list(
    genome_data = genome_data,
    markers_df = markers_df
  ))
}

#' Detect inversions between two genome arrangements
#' @param original_markers Original genome markers
#' @param inverted_markers Inverted genome markers
#' @return List with inversion detection results
detect_inversions <- function(original_markers, inverted_markers) {
  cat("Detecting inversions...\n")
  
  # Create comparison dataframe
  comparison_df <- tibble(
    original_position = 1:length(original_markers),
    inverted_position = match(original_markers, inverted_markers),
    marker_id = original_markers
  ) %>%
    filter(!is.na(inverted_position))
  
  # Detect inversion blocks
  inversions <- tibble(
    inversion_id = integer(),
    start_pos = integer(),
    end_pos = integer(),
    size = integer(),
    start_marker = character(),
    end_marker = character(),
    type = character()
  )
  
  # Simple inversion detection algorithm
  n <- nrow(comparison_df)
  if (n < 2) return(list(inversions = inversions, comparison = comparison_df))
  
  i <- 1
  inversion_count <- 0
  
  while (i < n) {
    # Look for decreasing sequences (potential inversions)
    if (comparison_df$inverted_position[i] > comparison_df$inverted_position[i + 1]) {
      inversion_count <- inversion_count + 1
      start_pos <- i
      
      # Find the end of this inversion
      j <- i + 1
      while (j <= n && 
             (j == n || comparison_df$inverted_position[j] < comparison_df$inverted_position[j - 1])) {
        j <- j + 1
      }
      end_pos <- j - 1
      
      # Record this inversion
      inversions <- bind_rows(inversions, tibble(
        inversion_id = inversion_count,
        start_pos = start_pos,
        end_pos = end_pos,
        size = end_pos - start_pos + 1,
        start_marker = comparison_df$marker_id[start_pos],
        end_marker = comparison_df$marker_id[end_pos],
        type = "inversion"
      ))
      
      i <- j
    } else {
      i <- i + 1
    }
  }
  
  cat("  Detected", nrow(inversions), "inversion blocks\n")
  
  return(list(
    inversions = inversions,
    comparison = comparison_df
  ))
}

#' Count inversion events per marker
#' @param inversions Inversion detection results
#' @param markers_df Markers dataframe
#' @return tibble with inversion frequencies
calculate_inversion_frequency <- function(inversions, markers_df) {
  cat("Calculating inversion frequencies...\n")
  
  if (nrow(inversions) == 0) {
    return(tibble(
      marker_id = character(),
      frequency = integer(),
      total_events = integer()
    ))
  }
  
  # Get all markers involved in inversions
  involved_markers <- c()
  for (i in 1:nrow(inversions)) {
    start_pos <- inversions$start_pos[i]
    end_pos <- inversions$end_pos[i]
    involved_markers <- c(involved_markers, 
                          markers_df$marker_id[start_pos:end_pos])
  }
  
  # Count frequencies
  frequency_df <- tibble(marker_id = involved_markers) %>%
    count(marker_id, name = "frequency") %>%
    mutate(total_events = nrow(inversions)) %>%
    arrange(desc(frequency))
  
  cat("  ", nrow(frequency_df), "markers involved in inversions\n")
  
  return(frequency_df)
}

#' Create comprehensive inversion analysis
#' @param inversions Detected inversions
#' @param markers_df Markers dataframe
#' @param config Analysis configuration
#' @return tibble with detailed analysis
create_detailed_analysis <- function(inversions, markers_df, config) {
  cat("Creating detailed analysis...\n")
  
  if (nrow(inversions) == 0) {
    return(tibble(
      inversion_id = integer(),
      chromosome = character(),
      start_bp = integer(),
      end_bp = integer(),
      size_markers = integer(),
      size_mb = numeric(),
      num_genes = integer(),
      inversion_rate_per_mb = numeric(),
      type = character()
    ))
  }
  
  # Enhance inversions with genomic coordinates
  detailed_analysis <- inversions %>%
    left_join(
      markers_df %>% select(linear_position, chromosome, start_bp, end_bp) %>%
        rename(start_pos = linear_position),
      by = "start_pos"
    ) %>%
    rename(start_chromosome = chromosome, start_bp_coord = start_bp, start_end_bp = end_bp) %>%
    left_join(
      markers_df %>% select(linear_position, chromosome, start_bp, end_bp) %>%
        rename(end_pos = linear_position),
      by = "end_pos"
    ) %>%
    rename(end_chromosome = chromosome, end_bp_coord = start_bp, end_end_bp = end_bp) %>%
    mutate(
      chromosome = start_chromosome,  # Assuming inversion within same chromosome
      start_bp = start_bp_coord,
      end_bp = end_end_bp,
      size_markers = size,
      size_mb = size * config$marker_size_bp / 1e6,
      num_genes = size,  # Assuming 1 gene per marker for simplicity
      inversion_rate_per_mb = 1 / size_mb,  # Inversions per Mb
      inversion_events = if_else(config$detect_inversion_events, 1, 0)
    ) %>%
    select(
      inversion_id, chromosome, start_bp, end_bp, 
      size_markers, size_mb, num_genes, inversion_rate_per_mb, 
      type, inversion_events
    )
  
  cat("  Detailed analysis for", nrow(detailed_analysis), "inversions\n")
  
  return(detailed_analysis)
}

#' Create comparison dotplot
#' @param comparison_df Comparison dataframe
#' @param output_path Output path for plot
#' @param title Plot title
create_comparison_dotplot <- function(comparison_df, output_path, title = "Genome Comparison") {
  cat("Creating comparison dotplot...\n")
  
  # Create the plot
  p <- ggplot(comparison_df, aes(x = original_position, y = inverted_position)) +
    geom_point(color = "steelblue", size = 0.8, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", alpha = 0.7) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = paste("Markers compared:", nrow(comparison_df)),
      x = "Original Genome Position",
      y = "Inverted Genome Position"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  # Save plot
  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat("Comparison dotplot saved to:", output_path, "\n")
  
  return(p)
}

#' Create inversion heatmap
#' @param frequency_df Frequency dataframe
#' @param markers_df Markers dataframe
#' @param output_path Output path for heatmap
create_inversion_heatmap <- function(frequency_df, markers_df, output_path) {
  cat("Creating inversion heatmap...\n")
  
  if (nrow(frequency_df) == 0) {
    cat("No inversions detected, skipping heatmap\n")
    return(NULL)
  }
  
  # Join frequency with position data
  heatmap_data <- markers_df %>%
    left_join(frequency_df, by = "marker_id") %>%
    replace_na(list(frequency = 0)) %>%
    mutate(
      chromosome = factor(chromosome),
      freq_category = case_when(
        frequency == 0 ~ "0",
        frequency == 1 ~ "1",
        frequency == 2 ~ "2",
        frequency >= 3 ~ "3+"
      )
    )
  
  # Create heatmap
  p <- ggplot(heatmap_data, aes(x = linear_position, y = chromosome, fill = freq_category)) +
    geom_tile(height = 0.8) +
    scale_fill_manual(
      values = c("0" = "white", "1" = "lightblue", "2" = "blue", "3+" = "darkblue"),
      name = "Inversion\nFrequency"
    ) +
    theme_minimal() +
    labs(
      title = "Inversion Frequency Heatmap",
      x = "Linear Genome Position",
      y = "Chromosome"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      panel.grid = element_blank()
    )
  
  # Save heatmap
  ggsave(output_path, p, width = 14, height = 6, dpi = 300)
  cat("Inversion heatmap saved to:", output_path, "\n")
  
  return(p)
}

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

#' Main function to run genome inversion analysis
#' @param config Analysis configuration
run_genome_inversion_analysis <- function(config = ANALYSIS_CONFIG) {
  cat("=== GENOME INVERSION ANALYZER ===\n")
  cat("Analyzing inversions between original and inverted genomes\n\n")
  
  # Step 1: Process both FASTA files
  cat("Step 1: Processing FASTA files\n")
  original_result <- process_fasta_to_markers(config$original_fasta_path, config$marker_size_bp)
  
  original_markers <- original_result$markers_df$marker_id
  
  # Check if marker order file exists (from simulation)
  use_marker_file <- isTRUE(config$use_marker_order_file) && 
    !is.null(config$marker_order_file) && 
    file.exists(config$marker_order_file)
  
  if (use_marker_file) {
    cat("Using marker order file:", config$marker_order_file, "\n")
    
    # Read the marker order from simulation
    marker_order_df <- read_csv(config$marker_order_file, show_col_types = FALSE)
    
    # Create comparison directly from the marker order file
    comparison_df <- marker_order_df %>%
      select(original_position, marker_id, inverted_position) %>%
      arrange(original_position)
    
    cat("  Loaded marker order with", nrow(comparison_df), "markers\n")
    
    # Quick check for inversions
    decreasing_count <- sum(diff(comparison_df$inverted_position) < 0)
    cat("  Found", decreasing_count, "decreasing positions in marker order\n")
    
  } else {
    cat("Independent mode: Processing both FASTA files\n")
    cat("Processing inverted FASTA file:", config$inverted_fasta_path, "\n")
    inverted_result <- process_fasta_to_markers(config$inverted_fasta_path, config$marker_size_bp)
    inverted_markers <- inverted_result$markers_df$marker_id
    
    # Create comparison dataframe between original and inverted
    cat("Creating comparison between original and inverted genomes...\n")
    
    # Method 1: Try to detect structural differences by comparing marker sequences
    # If FASTA files are truly different, we need a different approach
    
    # For now, create a simple comparison - this assumes the FASTA files contain
    # the same markers but potentially in different orders
    comparison_df <- tibble(
      original_position = 1:length(original_markers),
      inverted_position = match(original_markers, inverted_markers),
      marker_id = original_markers
    ) %>%
      filter(!is.na(inverted_position))
    
    # Check if this is a meaningful comparison
    identical_order <- all(comparison_df$original_position == comparison_df$inverted_position)
    
    if (identical_order) {
      cat("WARNING: Original and inverted genomes appear identical.\n")
      cat("For true independent analysis, FASTA files should contain rearranged sequences.\n")
      cat("Consider using marker order file from simulation, or provide truly rearranged FASTA.\n")
    }
  }
  
  cat("  Generated comparison for", nrow(comparison_df), "markers\n")
  
  # Step 2: Detect inversions using improved algorithm
  cat("\nStep 2: Detecting inversions\n")
  detection_result <- detect_inversions_improved(comparison_df)
  
  # Step 3: Calculate frequencies
  cat("\nStep 3: Calculating frequencies\n")
  frequency_result <- calculate_inversion_frequency_improved(
    detection_result$inversions, 
    comparison_df
  )
  
  # Step 4: Create detailed analysis
  cat("\nStep 4: Creating detailed analysis\n")
  detailed_analysis <- create_detailed_analysis_improved(
    detection_result$inversions,
    original_result$markers_df,
    config
  )
  
  # Step 5: Create visualizations
  cat("\nStep 5: Creating visualizations\n")
  comparison_plot <- create_comparison_dotplot(
    comparison_df,
    config$comparison_dotplot_png,
    "Original vs Inverted Genome Comparison"
  )
  
  heatmap_plot <- create_inversion_heatmap(
    frequency_result,
    original_result$markers_df,
    config$inversion_heatmap_png
  )
  
  # Step 6: Save CSV outputs
  cat("\nStep 6: Saving analysis results\n")
  
  # Inversion summary
  if (nrow(detailed_analysis) > 0) {
    write_csv(detailed_analysis, config$inversion_analyser_summary_csv)
    cat("Inversion summary saved to:", config$inversion_analyser_summary_csv, "\n")
    
    # Frequency analysis
    write_csv(frequency_result, config$inversion_analyser_frequency_csv)
    cat("Frequency analysis saved to:", config$inversion_analyser_frequency_csv, "\n")
    
    # Detailed analysis with all metrics
    comprehensive_analysis <- detailed_analysis %>%
      mutate(
        total_inversions = nrow(detailed_analysis),
        total_markers = length(original_markers),
        genome_size_mb = sum(original_result$markers_df$size_mb, na.rm = TRUE),
        inversion_density = total_inversions / genome_size_mb
      )
    
    write_csv(comprehensive_analysis, config$detailed_analysis_csv)
    cat("Detailed analysis saved to:", config$detailed_analysis_csv, "\n")
  } else {
    cat("No inversions detected - creating empty output files\n")
    
    # Create empty tibbles with proper column structure
    empty_summary <- tibble(
      inversion_id = integer(0),
      chromosome = character(0),
      start_bp = integer(0),
      end_bp = integer(0),
      size_markers = integer(0),
      size_mb = numeric(0),
      num_genes = integer(0),
      inversion_rate_per_mb = numeric(0),
      type = character(0),
      inversion_events = integer(0)
    )
    
    empty_frequency <- tibble(
      marker_id = character(0),
      frequency = integer(0),
      total_events = integer(0)
    )
    
    # Write empty files safely
    tryCatch({
      write_csv(empty_summary, config$inversion_analyser_summary_csv)
      cat("Empty summary saved to:", config$inversion_analyser_summary_csv, "\n")
    }, error = function(e) {
      cat("Error saving summary:", e$message, "\n")
    })
    
    tryCatch({
      write_csv(empty_frequency, config$inversion_analyser_frequency_csv)
      cat("Empty frequency saved to:", config$inversion_analyser_frequency_csv, "\n")
    }, error = function(e) {
      cat("Error saving frequency:", e$message, "\n")
    })
    
    tryCatch({
      write_csv(empty_summary, config$detailed_analysis_csv)
      cat("Empty detailed analysis saved to:", config$detailed_analysis_csv, "\n")
    }, error = function(e) {
      cat("Error saving detailed analysis:", e$message, "\n")
    })
  }
  
  # Step 7: Generate summary report
  cat("\n=== ANALYSIS COMPLETED ===\n")
  cat("Summary Statistics:\n")
  cat("  - Total inversions detected:", nrow(detection_result$inversions), "\n")
  cat("  - Total markers analyzed:", length(original_markers), "\n")
  cat("  - Markers involved in inversions:", nrow(frequency_result), "\n")
  cat("  - Average inversion size:", 
      if(nrow(detailed_analysis) > 0) round(mean(detailed_analysis$size_markers), 2) else 0, 
      "markers\n")
  
  if (nrow(detailed_analysis) > 0) {
    cat("  - Largest inversion:", max(detailed_analysis$size_markers), "markers\n")
    cat("  - Smallest inversion:", min(detailed_analysis$size_markers), "markers\n")
    cat("  - Total inverted sequence:", round(sum(detailed_analysis$size_mb), 2), "Mb\n")
  }
  
  cat("\nFiles created:\n")
  cat("  - Comparison dotplot:", config$comparison_dotplot_png, "\n")
  if (!is.null(heatmap_plot)) {
    cat("  - Inversion heatmap:", config$inversion_heatmap_png, "\n")
  }
  cat("  - Inversion summary:", config$inversion_analyser_summary_csv, "\n")
  cat("  - Frequency analysis:", config$inversion_analyser_frequency_csv, "\n")
  cat("  - Detailed analysis:", config$detailed_analysis_csv, "\n")
  
  return(list(
    config = config,
    original_data = original_result,
    detection_result = detection_result,
    frequency_result = frequency_result,
    detailed_analysis = detailed_analysis,
    comparison_plot = comparison_plot,
    heatmap_plot = heatmap_plot,
    comparison_df = comparison_df
  ))
}

################################################################################
# RUN THE ANALYSIS
################################################################################

# Execute the analysis
if (interactive() || !exists("SKIP_EXECUTION")) {
  analysis_results <- run_genome_inversion_analysis(ANALYSIS_CONFIG)
  
  # Print detailed summary
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("DETAILED INVERSION ANALYSIS SUMMARY\n")
  cat(paste(rep("=", 50), collapse = ""), "\n")
  
  if (nrow(analysis_results$detailed_analysis) > 0) {
    cat("Inversion Details:\n")
    print(analysis_results$detailed_analysis)
    
    cat("\nFrequency Analysis:\n")
    print(head(analysis_results$frequency_result, 10))
  } else {
    cat("No inversions detected in the comparison.\n")
  }
}