################################################################################
# COMPLETE GENOME INVERSION WORKFLOW
# This script demonstrates how to use both modules together
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
# WORKFLOW CONFIGURATION
################################################################################

WORKFLOW_CONFIG <- list(
  # Input files
  input_fasta = "GCA_910594885.2_idBibMarc1.2_genomic.fna",
  
  # Intermediate files (output from Module 1, input to Module 2)
  inverted_fasta = "Bibio_marci_inverted.fasta",
  simulation_dotplot = "simulation_dotplot.png",
  
  # Final output files from Module 2
  analysis_dotplot = "analysis_comparison_dotplot.png",
  analysis_heatmap = "analysis_inversion_heatmap.png",
  inversion_summary = "final_inversion_summary.csv",
  inversion_frequency = "final_inversion_frequency.csv",
  detailed_analysis = "final_detailed_analysis.csv",
  
  # Simulation parameters
  num_inversions = 3,
  marker_size_bp = 1000,
  
  # Analysis parameters
  detect_inversion_events = TRUE,
  min_inversion_size = 2
)

################################################################################
# WORKFLOW EXECUTION FUNCTIONS
################################################################################

#' Run the complete workflow
#' @param config Workflow configuration
run_complete_workflow <- function(config = WORKFLOW_CONFIG) {
  cat("################################################################################\n")
  cat("COMPLETE GENOME INVERSION WORKFLOW\n")
  cat("################################################################################\n\n")
  
  # ============================================================================
  # PHASE 1: SIMULATION (Module 1)
  # ============================================================================
  cat("PHASE 1: GENOME INVERSION SIMULATION\n")
  cat("="*50 + "\n")
  
  # Configure Module 1
  simulation_config <- list(
    input_fasta_path = config$input_fasta,
    inverted_fasta_output = config$inverted_fasta,
    dotplot_output = config$simulation_dotplot,
    simulate_inversions = TRUE,
    num_inversions = config$num_inversions,
    marker_size_bp = config$marker_size_bp,
    use_busco = FALSE,
    save_intermediate_states = TRUE
  )
  
  # Source Module 1 functions (assuming they're available)
  # In practice, you would source("module1_genome_inversion_simulator.R")
  
  cat("Running inversion simulation...\n")
  cat("Input file:", config$input_fasta, "\n")
  cat("Number of inversions:", config$num_inversions, "\n")
  cat("Marker size:", config$marker_size_bp, "bp\n\n")
  
  # Run simulation
  simulation_results <- run_genome_inversion_simulation(simulation_config)
  
  cat("Simulation completed successfully!\n")
  cat("Generated files:\n")
  cat("  - Inverted FASTA:", config$inverted_fasta, "\n")
  cat("  - Simulation dotplot:", config$simulation_dotplot, "\n\n")
  
  # ============================================================================
  # PHASE 2: ANALYSIS (Module 2)
  # ============================================================================
  cat("PHASE 2: GENOME INVERSION ANALYSIS\n")
  cat("="*50 + "\n")
  
  # Configure Module 2
  analysis_config <- list(
    original_fasta_path = config$input_fasta,
    inverted_fasta_path = config$inverted_fasta,
    inversion_summary_csv = config$inversion_summary,
    inversion_frequency_csv = config$inversion_frequency,
    detailed_analysis_csv = config$detailed_analysis,
    comparison_dotplot_png = config$analysis_dotplot,
    inversion_heatmap_png = config$analysis_heatmap,
    marker_size_bp = config$marker_size_bp,
    detect_inversion_events = config$detect_inversion_events,
    min_inversion_size = config$min_inversion_size,
    plot_width = 12,
    plot_height = 8,
    dpi = 300
  )
  
  # Source Module 2 functions (assuming they're available)
  # In practice, you would source("module2_genome_inversion_analyzer.R")
  
  cat("Running inversion analysis...\n")
  cat("Comparing:", config$input_fasta, "vs", config$inverted_fasta, "\n")
  cat("Inversion event detection:", config$detect_inversion_events, "\n\n")
  
  # Run analysis
  analysis_results <- run_genome_inversion_analysis(analysis_config)
  
  cat("Analysis completed successfully!\n")
  cat("Generated files:\n")
  cat("  - Analysis dotplot:", config$analysis_dotplot, "\n")
  cat("  - Inversion heatmap:", config$analysis_heatmap, "\n")
  cat("  - Inversion summary:", config$inversion_summary, "\n")
  cat("  - Frequency analysis:", config$inversion_frequency, "\n")
  cat("  - Detailed analysis:", config$detailed_analysis, "\n\n")
  
  # ============================================================================
  # PHASE 3: COMPREHENSIVE SUMMARY
  # ============================================================================
  cat("PHASE 3: COMPREHENSIVE WORKFLOW SUMMARY\n")
  cat("="*50 + "\n")
  
  # Generate comprehensive summary
  generate_workflow_summary(simulation_results, analysis_results, config)
  
  return(list(
    config = config,
    simulation_results = simulation_results,
    analysis_results = analysis_results
  ))
}

#' Generate comprehensive workflow summary
#' @param sim_results Results from simulation
#' @param ana_results Results from analysis
#' @param config Workflow configuration
generate_workflow_summary <- function(sim_results, ana_results, config) {
  cat("WORKFLOW EXECUTION SUMMARY\n")
  cat("-" * 30 + "\n")
  
  # Simulation summary
  cat("SIMULATION PHASE:\n")
  if (!is.null(sim_results$simulation_result)) {
    cat("  ✓ Simulated inversions:", nrow(sim_results$simulation_result$inversions_log), "\n")
    cat("  ✓ Total markers processed:", length(sim_results$simulation_result$original_genome), "\n")
    cat("  ✓ Genome size:", sum(sim_results$genome_data$chromosome_size_bp), "bp\n")
  }
  
  # Analysis summary
  cat("\nANALYSIS PHASE:\n")
  if (!is.null(ana_results$detailed_analysis)) {
    cat("  ✓ Inversions detected:", nrow(ana_results$detailed_analysis), "\n")
    cat("  ✓ Markers involved:", nrow(ana_results$frequency_result), "\n")
    if (nrow(ana_results$detailed_analysis) > 0) {
      cat("  ✓ Average inversion size:", round(mean(ana_results$detailed_analysis$size_markers), 2), "markers\n")
      cat("  ✓ Total inverted sequence:", round(sum(ana_results$detailed_analysis$size_mb), 2), "Mb\n")
    }
  }
  
  # File inventory
  cat("\nGENERATED FILES:\n")
  files_to_check <- c(
    config$inverted_fasta,
    config$simulation_dotplot,
    config$analysis_dotplot,
    config$analysis_heatmap,
    config$inversion_summary,
    config$inversion_frequency,
    config$detailed_analysis
  )
  
  for (file in files_to_check) {
    if (file.exists(file)) {
      cat("  ✓", file, "\n")
    } else {
      cat("  ✗", file, "(missing)\n")
    }
  }
  
  cat("\nWORKFLOW COMPLETED SUCCESSFULLY!\n")
}

################################################################################
# BATCH PROCESSING FUNCTIONS
################################################################################

#' Run workflow with multiple inversion scenarios
#' @param base_config Base configuration
#' @param inversion_counts Vector of inversion counts to test
run_batch_analysis <- function(base_config = WORKFLOW_CONFIG, inversion_counts = c(1, 3, 5, 10)) {
  cat("BATCH ANALYSIS: Testing multiple inversion scenarios\n")
  cat("="*60 + "\n")
  
  batch_results <- list()
  
  for (i in seq_along(inversion_counts)) {
    n_inv <- inversion_counts[i]
    cat("\nBATCH", i, "of", length(inversion_counts), ": Testing", n_inv, "inversions\n")
    cat("-" * 40 + "\n")
    
    # Modify configuration for this batch
    batch_config <- base_config
    batch_config$num_inversions <- n_inv
    batch_config$inverted_fasta <- paste0("batch_", n_inv, "_inverted.fasta")
    batch_config$simulation_dotplot <- paste0("batch_", n_inv, "_simulation.png")
    batch_config$analysis_dotplot <- paste0("batch_", n_inv, "_analysis.png")
    batch_config$analysis_heatmap <- paste0("batch_", n_inv, "_heatmap.png")
    batch_config$inversion_summary <- paste0("batch_", n_inv, "_summary.csv")
    batch_config$inversion_frequency <- paste0("batch_", n_inv, "_frequency.csv")
    batch_config$detailed_analysis <- paste0("batch_", n_inv, "_detailed.csv")
    
    # Run workflow
    batch_results[[i]] <- run_complete_workflow(batch_config)
    
    cat("Batch", i, "completed\n")
  }
  
  cat("\nBATCH ANALYSIS COMPLETED\n")
  cat("Processed", length(inversion_counts), "scenarios\n")
  
  return(batch_results)
}

#' Compare results across different inversion counts
#' @param batch_results Results from batch analysis
compare_batch_results <- function(batch_results) {
  cat("BATCH COMPARISON SUMMARY\n")
  cat("="*40 + "\n")
  
  comparison_df <- tibble(
    scenario = integer(),
    simulated_inversions = integer(),
    detected_inversions = integer(),
    total_markers = integer(),
    inverted_markers = integer(),
    avg_inversion_size = numeric(),
    total_inverted_mb = numeric()
  )
  
  for (i in seq_along(batch_results)) {
    result <- batch_results[[i]]
    
    comparison_df <- bind_rows(comparison_df, tibble(
      scenario = i,
      simulated_inversions = result$config$num_inversions,
      detected_inversions = nrow(result$analysis_results$detailed_analysis),
      total_markers = length(result$simulation_results$simulation_result$original_genome),
      inverted_markers = nrow(result$analysis_results$frequency_result),
      avg_inversion_size = if(nrow(result$analysis_results$detailed_analysis) > 0) 
        round(mean(result$analysis_results$detailed_analysis$size_markers), 2) else 0,
      total_inverted_mb = if(nrow(result$analysis_results$detailed_analysis) > 0) 
        round(sum(result$analysis_results$detailed_analysis$size_mb), 2) else 0
    ))
  }
  
  print(comparison_df)
  
  # Save comparison
  write_csv(comparison_df, "batch_comparison_summary.csv")
  cat("\nBatch comparison saved to: batch_comparison_summary.csv\n")
  
  return(comparison_df)
}

################################################################################
# EXAMPLE USAGE
################################################################################

# Example 1: Run single workflow
cat("EXAMPLE 1: Single workflow execution\n")
cat("="*40 + "\n")

# Uncomment to run:
# workflow_results <- run_complete_workflow(WORKFLOW_CONFIG)

# Example 2: Run batch analysis
cat("\nEXAMPLE 2: Batch analysis\n")
cat("="*40 + "\n")

# Uncomment to run:
# batch_results <- run_batch_analysis(WORKFLOW_CONFIG, c(1, 2, 5))
# comparison_summary <- compare_batch_results(batch_results)

cat("\nTo execute the workflow, uncomment the desired example above.\n")
cat("Make sure both Module 1 and Module 2 scripts are available in your environment.\n")