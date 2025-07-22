# Comprehensive Pairwise Species Comparison Analysis
# Generates dotplots for all possible species combinations
################################################################################

# Load required libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gsheet)
library(Biostrings)
library(knitr)

################################################################################
# Setup
################################################################################
# Set working directory to current location
setwd("~/Documents/Bibionidae/")
root <- getwd()
if (!endsWith(root, "/")) root <- paste0(root, "/")

# Color palette for ALGs
pal <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
         "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", 
         "unassigned" = "grey")

# Create output directories
figures_dir <- file.path(root, "figures")
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir)
}

dotplots_dir <- file.path(root, "figures", "dotplots")
if (!dir.exists(dotplots_dir)) {
  dir.create(dotplots_dir)
}

summaries_dir <- file.path(root, "figures", "summaries")
if (!dir.exists(summaries_dir)) {
  dir.create(summaries_dir)
}

################################################################################
# Function to read BUSCO data with validation
################################################################################
read_valid_buscos <- function(filepath, species_name) {
  if (!file.exists(filepath)) {
    warning(paste("BUSCO file not found:", filepath))
    return(NULL)
  }
  
  lines <- readLines(filepath)
  valid_lines <- lines[sapply(strsplit(lines, "\t"), length) == 8]
  
  if (length(valid_lines) == 0) {
    warning(paste("No valid BUSCO data in:", filepath))
    return(NULL)
  }
  
  df <- read.table(text = valid_lines, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("marker", "status", "chromosome", "start", "end", "strand", "score", "length")
  df <- df[df$status == "Complete", ]
  df$species <- species_name
  df <- df[, c("marker", "chromosome", "start", "end", "species")]
  return(df)
}

################################################################################
# Function to get scaffold sizes from FASTA
################################################################################
get_scaffold_sizes <- function(fasta_path, species_name) {
  if (!file.exists(fasta_path)) {
    warning(paste("FASTA file not found:", fasta_path))
    return(NULL)
  }
  
  seqs <- readDNAStringSet(fasta_path)
  tibble::tibble(
    species = species_name,
    chromosome = names(seqs),
    chromsome_size_b = as.numeric(width(seqs))
  )
}

################################################################################
# Main Analysis Function
################################################################################
create_species_dotplot <- function(species1, species2, all_data, busco_data, syngraph_clean, output_prefix) {
  
  cat("\n=== Creating dotplot for:", species1, "vs", species2, "===\n")
  
  # Filter data for target species
  target_species <- c(species1, species2)
  target_species_data <- all_data %>% 
    filter(species %in% target_species)
  
  if (nrow(target_species_data) == 0) {
    cat("No genome data found for these species.\n")
    return(NULL)
  }
  
  # Merge with BUSCO data
  target_species_data <- left_join(target_species_data, busco_data, 
                                   by = c("species", "chromosome"))
  
  # Create n3-style mapping from syngraph
  n3 <- syngraph_clean %>%
    select(marker = 1, ALG) %>%
    distinct()
  
  # Merge with ALG assignments
  target_species_data <- left_join(target_species_data, n3, by = "marker")
  
  # Remove chromosomes without markers
  target_species_data <- target_species_data %>% 
    filter(!is.na(marker))
  
  if (nrow(target_species_data) == 0) {
    cat("No shared markers found between species.\n")
    return(NULL)
  }
  
  # Mark unassigned BUSCOs
  target_species_data$ALG[is.na(target_species_data$ALG)] <- "unassigned"
  
  # Convert to linear coordinates
  df_lin <- target_species_data %>% arrange(species, chromosome, start)
  df_lin$linear_start <- NA
  df_lin$linear_end <- NA
  
  chr_offset <- 0
  current_species <- ""
  
  for (i in 1:nrow(df_lin)) {
    if (i == 1 || df_lin$species[i] != df_lin$species[i - 1]) {
      chr_offset <- 0
      current_species <- df_lin$species[i]
    } else if (df_lin$chromosome[i] != df_lin$chromosome[i - 1]) {
      chr_offset <- chr_offset + as.integer(df_lin$chromsome_size_b[i - 1])
    }
    df_lin$linear_start[i] <- chr_offset + df_lin$start[i]
    df_lin$linear_end[i] <- chr_offset + df_lin$end[i]
  }
  
  # Remove duplicates
  duplicates <- df_lin %>%
    summarise(n = n(), .by = c(marker, species)) %>%
    filter(n > 1L)
  
  # Create ALG dataframe
  alg_df <- df_lin %>%
    select(marker, ALG) %>%
    distinct()
  
  # Pivot wider
  df_wide <- df_lin %>%
    filter(!marker %in% duplicates$marker) %>%
    select(marker, species, chromosome, linear_start, linear_end) %>%
    pivot_wider(
      id_cols = marker,
      names_from = species,
      values_from = c(chromosome, linear_start, linear_end),
      names_sep = "_"
    ) %>%
    left_join(alg_df, by = "marker") %>%
    relocate(ALG, .after = marker)
  
  # Filter for complete data
  start_x <- paste0("linear_start_", species1)
  start_y <- paste0("linear_start_", species2)
  
  df_wide_clean <- df_wide %>%
    filter(!is.na(!!sym(start_x)) & !is.na(!!sym(start_y)))
  
  if (nrow(df_wide_clean) == 0) {
    cat("No complete marker pairs found.\n")
    return(NULL)
  }
  
  # Chromosome info for grid lines
  chr_info_x <- df_lin %>%
    filter(species == species1) %>%
    arrange(chromosome) %>%
    distinct(chromosome, chromsome_size_b) %>%
    mutate(cum_end = cumsum(chromsome_size_b))
  
  chr_info_y <- df_lin %>%
    filter(species == species2) %>%
    arrange(chromosome) %>%
    distinct(chromosome, chromsome_size_b) %>%
    mutate(cum_end = cumsum(chromsome_size_b))
  
  # Chromosome labels
  chr_labels_x <- df_lin %>%
    filter(species == species1) %>%
    group_by(chromosome) %>%
    summarise(mid = mean(c(min(linear_start), max(linear_end))), .groups = 'drop')
  
  chr_labels_y <- df_lin %>%
    filter(species == species2) %>%
    group_by(chromosome) %>%
    summarise(mid = mean(c(min(linear_start), max(linear_end))), .groups = 'drop')
  
  # Calculate correlation
  cor_value <- cor(df_wide_clean[[start_x]], df_wide_clean[[start_y]], use = "complete.obs")
  
  # Create plot
  p <- ggplot(df_wide_clean) + 
    geom_point(aes(
      x = !!sym(start_x), 
      y = !!sym(start_y), 
      fill = ALG, 
      colour = ALG
    ), size = 2, alpha = 0.8) +
    geom_vline(data = chr_info_x, aes(xintercept = cum_end), 
               color = "grey80", linewidth = 0.2) +
    geom_hline(data = chr_info_y, aes(yintercept = cum_end), 
               color = "grey80", linewidth = 0.2) +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      expand = c(0.01, 0),
      labels = function(x) x / 1e6,
      name = paste0("\n", gsub("_", " ", species1), " (Mb)"),
      sec.axis = sec_axis(~., breaks = chr_labels_x$mid, labels = chr_labels_x$chromosome)
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      labels = function(y) y / 1e6,
      name = paste0(gsub("_", " ", species2), " (Mb)\n"),
      sec.axis = sec_axis(~., breaks = chr_labels_y$mid, labels = chr_labels_y$chromosome)
    ) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text.x.top = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5),
      axis.text.x.bottom = element_text(size = 12),
      axis.title.x.bottom = element_text(size = 14),
      axis.text.y.right = element_text(size = 10),
      axis.text.y.left = element_text(size = 12),
      axis.title.y.left = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    labs(
      title = paste("Synteny:", species1, "vs", species2),
      subtitle = paste("Correlation R =", round(cor_value, 3), 
                       "| Markers:", nrow(df_wide_clean))
    )
  
  # Save plot
  filename <- paste0(output_prefix, gsub(" ", "_", species1), "_vs_", gsub(" ", "_", species2))
  
  png_path <- file.path(root, "figures", "dotplots", paste0(filename, ".png"))
  svg_path <- file.path(root, "figures", "dotplots", paste0(filename, ".svg"))
  
  ggsave(png_path, plot = p, dpi = 300, height = 8, width = 10)
  ggsave(svg_path, plot = p, dpi = 300, height = 8, width = 10)
  
  # Return summary statistics
  return(list(
    species1 = species1,
    species2 = species2,
    n_markers = nrow(df_wide_clean),
    correlation = cor_value,
    n_chr_sp1 = length(unique(chr_info_x$chromosome)),
    n_chr_sp2 = length(unique(chr_info_y$chromosome)),
    genome_size_sp1 = sum(chr_info_x$chromsome_size_b) / 1e6,
    genome_size_sp2 = sum(chr_info_y$chromsome_size_b) / 1e6
  ))
}

################################################################################
# MAIN EXECUTION
################################################################################

cat("Starting comprehensive pairwise species comparison analysis...\n")

# Read Google Sheets data
cat("\nReading genome data from Google Sheets...\n")
all_genome_data <- read.csv(
  text = gsheet2text("https://docs.google.com/spreadsheets/d/1K01wVWkMW-m6yT9zDX8gDekp-OECubE-9HcmD8RnmkM/edit?gid=1940964825#gid=1940964825", 
                     format='csv'),
  stringsAsFactors = FALSE, 
  header = TRUE, 
  check.names = FALSE
)

# Get list of available species from your files
available_species <- c("Bibio_marci", "Dilophus_febrilis", "Dioctria_linearis", 
                       "Dioctria_rufipes", "Plecia_longiforceps")

# Filter for available species
genome_data_filtered <- all_genome_data %>%
  filter(species %in% available_species)

# For species with FASTA files, update with actual sizes
cat("\nReading FASTA files for accurate genome sizes...\n")

# Define FASTA files
fasta_files <- list(
  "Bibio_marci" = "GCA_910594885.2_idBibMarc1.2_genomic.fna",
  "Dilophus_febrilis" = "GCA_958336335.1_idDilFebr1.1_genomic.fna"
)

# Read FASTA data
fasta_data <- list()
for (species in names(fasta_files)) {
  if (file.exists(fasta_files[[species]])) {
    fasta_data[[species]] <- get_scaffold_sizes(fasta_files[[species]], species)
  }
}

# Combine FASTA data
if (length(fasta_data) > 0) {
  fasta_combined <- bind_rows(fasta_data)
  fasta_combined$chromosome <- sub(" .*", "", fasta_combined$chromosome)
  
  # Update genome data with FASTA information
  genome_data_final <- genome_data_filtered %>%
    filter(!species %in% names(fasta_files)) %>%
    bind_rows(fasta_combined)
} else {
  genome_data_final <- genome_data_filtered
}

# Load all BUSCO data
cat("\nLoading BUSCO data...\n")
busco_data_list <- list()
for (species in available_species) {
  busco_file <- paste0(species, ".tsv")
  if (file.exists(busco_file)) {
    busco_data_list[[species]] <- read_valid_buscos(busco_file, species)
  }
}

# Combine BUSCO data
all_busco_data <- bind_rows(busco_data_list)

# Load syngraph data
cat("\nLoading synteny graph data...\n")
syngraph_file <- file.path(root, "diptera.syngraph_tabulate.table.tsv")
if (!file.exists(syngraph_file)) {
  # Alternative file name
  syngraph_file <- file.path(root, "diptera.syngraph_tabulate.tsv")
}

if (!file.exists(syngraph_file)) {
  stop("Syngraph file not found. Please check that 'diptera.syngraph_tabulate.table.tsv' exists in your working directory.")
}

syngraph <- read.delim(syngraph_file, header = TRUE, sep = "\t", check.names = FALSE)

# Process syngraph for available species
syngraph_cols <- colnames(syngraph)
species_cols <- list()

for (species in available_species) {
  cols <- c(paste0(species, "_seq"), paste0(species, "_start"), paste0(species, "_end"))
  if (all(cols %in% syngraph_cols)) {
    species_cols[[species]] <- cols
  }
}

# Create ALG assignments
syngraph_clean <- syngraph
syngraph_clean$ALG <- paste0("M", seq_len(nrow(syngraph_clean)))

# Generate all pairwise combinations
species_with_data <- names(busco_data_list)
comparisons <- combn(species_with_data, 2, simplify = FALSE)

cat("\nGenerating", length(comparisons), "pairwise comparisons...\n")

# Run all comparisons
results <- list()
for (i in seq_along(comparisons)) {
  sp1 <- comparisons[[i]][1]
  sp2 <- comparisons[[i]][2]
  
  result <- create_species_dotplot(
    sp1, sp2, 
    genome_data_final, 
    all_busco_data, 
    syngraph_clean,
    output_prefix = paste0("dotplot_", i, "_")
  )
  
  if (!is.null(result)) {
    results[[length(results) + 1]] <- result
  }
}

# Create summary report
cat("\n\nCreating summary report...\n")

summary_df <- bind_rows(results)

# Generate markdown report
report <- c(
  "# Pairwise Species Synteny Analysis Summary",
  paste("\nGenerated on:", Sys.Date()),
  paste("\nTotal comparisons:", nrow(summary_df)),
  "\n## Overview",
  "\nThis analysis compared genome synteny between all pairs of available Bibionidae species using BUSCO markers.",
  "\n## Summary Statistics\n"
)

# Add summary table
summary_table <- summary_df %>%
  arrange(desc(correlation)) %>%
  mutate(
    Comparison = paste(species1, "vs", species2),
    Markers = n_markers,
    Correlation = round(correlation, 3),
    `Genome Size 1 (Mb)` = round(genome_size_sp1, 1),
    `Genome Size 2 (Mb)` = round(genome_size_sp2, 1)
  ) %>%
  select(Comparison, Markers, Correlation, 
         `Genome Size 1 (Mb)`, `Genome Size 2 (Mb)`)

report <- c(report, 
            knitr::kable(summary_table, format = "markdown")
)

# Add interpretation
report <- c(report,
            "\n## Interpretation Guide",
            "\n### Correlation Values:",
            "- **R > 0.9**: Very high synteny, likely recent divergence",
            "- **R = 0.7-0.9**: High synteny with some rearrangements",
            "- **R = 0.5-0.7**: Moderate synteny, significant chromosomal evolution",
            "- **R < 0.5**: Low synteny, extensive rearrangements",
            "\n### Dotplot Patterns:",
            "- **Diagonal lines**: Conserved syntenic blocks",
            "- **Scattered dots**: Translocations or assembly issues",
            "- **Perpendicular patterns**: Inversions",
            "- **Gaps**: Missing data or lineage-specific regions"
)

# Species-specific notes
if (any(grepl("Bibio_marci.*Dilophus_febrilis", summary_df$Comparison))) {
  report <- c(report,
              "\n### Bibio marci vs Dilophus febrilis:",
              "These genera represent different subfamilies within Bibionidae, showing moderate synteny conservation."
  )
}

# Write report
report_path <- file.path(root, "figures", "synteny_analysis_summary.md")
writeLines(report, report_path)

# Also save as CSV
csv_path <- file.path(root, "figures", "synteny_summary_stats.csv")
write.csv(summary_df, csv_path, row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Dotplots saved in:", file.path(root, "figures", "dotplots"), "\n")
cat("Summary report:", report_path, "\n")
cat("Summary statistics:", csv_path, "\n")

# Print final summary
cat("\nComparisons completed:\n")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %s vs %s: %d markers, R = %.3f\n",
              summary_df$species1[i], 
              summary_df$species2[i],
              summary_df$n_markers[i],
              summary_df$correlation[i]))
}