# Figure 1A - Bibio_marci vs Dilophus_febrilis analyses
################################################################################
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Biostrings)

################################################################################
# Set working directory
root <- paste0(getwd(), "/")

################################################################################
# Color palette
pal <- c("M1" = "#1573afff", "M2" = "#e59d38ff", "M3" = "#f0e354ff", 
         "M4" = "#169e73ff", "M5" = "#60b5e1ff", "M6" = "black", "unassigned" = "grey")

target_species <- c("Bibio_marci", "Dilophus_febrilis")

################################################################################
# Function to get scaffold sizes from a FASTA file
get_scaffold_sizes <- function(fasta_path) {
  seqs <- readDNAStringSet(fasta_path)
  tibble::tibble(
    chromosome = names(seqs),
    chromsome_size_b = as.numeric(width(seqs))
  )
}

# Get sizes for each species
sizes_bibio <- get_scaffold_sizes("GCA_910594885.2_idBibMarc1.2_genomic.fna")
sizes_dilophus <- get_scaffold_sizes("GCA_958336335.1_idDilFebr1.1_genomic.fna")

# Combine into single dataframe
all_genome_data <- bind_rows(
  tibble(species = "Bibio_marci", accession = "GCA_910594885.2_idBibMarc1.2") %>%
    bind_cols(sizes_bibio),
  tibble(species = "Dilophus_febrilis", accession = "GCA_958336335.1_idDilFebr1.1") %>%
    bind_cols(sizes_dilophus)
)

# Filter for target species
target_species_data <- all_genome_data %>% 
  filter(species %in% target_species) %>%
  select(species, accession, chromosome, chromsome_size_b)

# Clean chromosome names
target_species_data$chromosome <- sub(" .*", "", target_species_data$chromosome)
target_species_data$accession <- sub("\\.[0-9]", "", target_species_data$accession)

################################################################################
# Load BUSCO data
busco_files_list <- c("Bibio_marci.tsv", "Dilophus_febrilis.tsv")

read_valid_buscos <- function(filepath, species_name) {
  lines <- readLines(filepath)
  valid_lines <- lines[sapply(strsplit(lines, "\t"), length) == 8]
  df <- read.table(text = valid_lines, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("marker", "status", "chromosome", "start", "end", "strand", "score", "length")
  df <- df[df$status == "Complete", ]  # filter to only "Complete"
  df$species <- species_name
  df <- df[, c("marker", "chromosome", "start", "end", "species")]
  return(df)
}

buscos <- bind_rows(
  read_valid_buscos(busco_files_list[1], "Bibio_marci"),
  read_valid_buscos(busco_files_list[2], "Dilophus_febrilis")
)

################################################################################
# Load syngraph table
syngraph <- read.delim("diptera.syngraph_tabulate.table.tsv", header = TRUE, sep = "\t", check.names = FALSE)

# Extract BUSCO marker IDs
markers <- syngraph[[1]]  # 1st column = BUSCO IDs

# Keep only Bibio + Dilophus columns + marker ID
cols_to_keep <- c("Bibio_marci_seq", "Bibio_marci_start", "Bibio_marci_end", 
                  "Dilophus_febrilis_seq", "Dilophus_febrilis_start", "Dilophus_febrilis_end")

syngraph_subset <- syngraph[, c(1, which(colnames(syngraph) %in% cols_to_keep))]

# Filter rows where both species are present
syngraph_clean <- syngraph_subset %>%
  filter(!is.na(Bibio_marci_seq) & !is.na(Dilophus_febrilis_seq))

# Create artificial ALG clusters
syngraph_clean$marker <- markers[match(syngraph_clean[[1]], markers)]
syngraph_clean$ALG <- paste0("M", seq_len(nrow(syngraph_clean)))

# Create n3-style mapping (marker vs ALG)
n3 <- syngraph_clean[, c("marker", "ALG")]

################################################################################
# Merge everything with the main table
target_species_data <- left_join(target_species_data, buscos, by = c("species", "chromosome"))
target_species_data <- left_join(target_species_data, n3, by = "marker")

# Remove chromosomes that have no markers
target_species_data <- target_species_data %>% filter(!is.na(marker))

# Mark buscos that are not in ALGs as unassigned
target_species_data$ALG[is.na(target_species_data$ALG)] <- "unassigned"

################################################################################
# Convert coordinates into linear format
df_lin <- target_species_data %>% arrange(species, chromosome, start)
df_lin$linear_start <- NA
df_lin$linear_end <- NA

# Calculate linear positions
chr_offset <- 0
current_species <- ""

for (i in 1:nrow(df_lin)) {
  if (i == 1 || df_lin$species[i] != df_lin$species[i - 1]) {
    # New species → reset offset
    chr_offset <- 0
    current_species <- df_lin$species[i]
  } else if (df_lin$chromosome[i] != df_lin$chromosome[i - 1]) {
    # Same species but new chromosome → increase offset
    chr_offset <- chr_offset + as.integer(df_lin$chromsome_size_b[i - 1])
  }
  
  # Assign linearized positions
  df_lin$linear_start[i] <- chr_offset + df_lin$start[i]
  df_lin$linear_end[i]   <- chr_offset + df_lin$end[i]
}

# Remove duplicated markers
duplicates <- df_lin %>%
  summarise(n = n(), .by = c(marker, species)) %>%
  filter(n > 1L)

# Create ALG dataframe for merging
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

################################################################################
# Plotting
sp_x <- target_species[1]
sp_y <- target_species[2]

# Filter chromosome info for species X and Y
chr_info_x <- df_lin %>%
  filter(species == sp_x) %>%
  arrange(chromosome) %>%
  distinct(chromosome, chromsome_size_b) %>%
  mutate(cum_end = cumsum(chromsome_size_b))

chr_info_y <- df_lin %>%
  filter(species == sp_y) %>%
  arrange(chromosome) %>%
  distinct(chromosome, chromsome_size_b) %>%
  mutate(cum_end = cumsum(chromsome_size_b))

# Chromosome label midpoints
chr_labels_x <- df_lin %>%
  filter(species == sp_x) %>%
  group_by(chromosome) %>%
  summarise(mid = mean(c(min(linear_start), max(linear_end))), .groups = 'drop')

chr_labels_y <- df_lin %>%
  filter(species == sp_y) %>%
  group_by(chromosome) %>%
  summarise(mid = mean(c(min(linear_start), max(linear_end))), .groups = 'drop')

# Build column names dynamically
start_x <- paste0("linear_start_", sp_x)
end_x   <- paste0("linear_end_", sp_x)
start_y <- paste0("linear_start_", sp_y)
end_y   <- paste0("linear_end_", sp_y)

# Buffer for positioning
buffer <- 1e5

# Remove any rows with NA values in coordinates
df_wide_clean <- df_wide %>%
  filter(!is.na(!!sym(start_x)) & !is.na(!!sym(start_y)))

# Build plot
p <- ggplot(df_wide_clean) + 
  geom_point(aes(
    x = !!sym(start_x) - buffer, 
    y = !!sym(start_y) - buffer, 
    fill = ALG, 
    colour = ALG
  ), size = 2) +
  geom_vline(data = chr_info_x, aes(xintercept = cum_end), color = "grey80", linewidth = 0.2) +
  geom_hline(data = chr_info_y, aes(yintercept = cum_end), color = "grey80", linewidth = 0.2) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    expand = c(0,0),
    labels = function(x) x / 1e6,
    name = paste0("\n", gsub("_", " ", sp_x), " (Mb)"),
    sec.axis = sec_axis(~., breaks = chr_labels_x$mid, labels = chr_labels_x$chromosome)
  ) +
  scale_y_continuous(
    expand = c(0,0),
    labels = function(y) y / 1e6,
    name = paste0(gsub("_", " ", sp_y), " (Mb)\n"),
    sec.axis = sec_axis(~., breaks = chr_labels_y$mid, labels = chr_labels_y$chromosome)
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x.top = element_text(angle = 90, size = 14, hjust = 0, vjust = 0.5),
    axis.text.x.bottom = element_text(size = 14),
    axis.title.x.bottom = element_text(size = 18),
    axis.text.y.right = element_text(size = 14),
    axis.text.y.left = element_text(size = 14),
    axis.title.y.left = element_text(size = 18),
    legend.position = "none"
  )

# Create figures directory if it doesn't exist
if (!dir.exists(paste0(root, "figures"))) {
  dir.create(paste0(root, "figures"))
}

# Save plots
ggsave(paste0(root, "figures/two_species_dotplot.svg"), plot = p, dpi = 600, height = 10.25, width = 11.58)
ggsave(paste0(root, "figures/two_species_dotplot.png"), plot = p, dpi = 600, height = 10.25, width = 11.58)

# Print summary
cat("Plot saved successfully!\n")
cat("Number of markers plotted:", nrow(df_wide_clean), "\n")
cat("ALG categories:", unique(df_wide_clean$ALG), "\n")