# Save as fixed_extract_synteny.R
library(dplyr)
library(ggplot2)
library(stringr)

# Fixed data loading function
load_syngraph_data <- function() {
  cat("Loading syngraph data...\n")
  
  # Load rearrangement data
  rearr_df <- read.delim('multi_species_analysis/diptera_results_clean.rearrangements.tsv', sep='\t')
  
  # Load BUSCO data properly (skip comments and use correct headers)
  read_busco_file <- function(file_path) {
    # Read all lines
    lines <- readLines(file_path)
    
    # Find the header line (starts with "# Busco id")
    header_line <- grep("^# Busco id", lines)
    
    if(length(header_line) > 0) {
      # Skip to data lines (after header)
      data_lines <- lines[(header_line + 1):length(lines)]
      # Remove any remaining comment lines
      data_lines <- data_lines[!grepl("^#", data_lines)]
      # Remove empty lines
      data_lines <- data_lines[data_lines != ""]
      
      # Parse the data
      data_list <- strsplit(data_lines, "\t")
      # Ensure all rows have same number of columns
      valid_rows <- sapply(data_list, length) >= 8
      data_list <- data_list[valid_rows]
      
      if(length(data_list) > 0) {
        # Convert to data frame
        busco_df <- data.frame(
          busco_id = sapply(data_list, function(x) x[1]),
          status = sapply(data_list, function(x) x[2]),
          sequence = sapply(data_list, function(x) x[3]),
          start = as.numeric(sapply(data_list, function(x) x[4])),
          end = as.numeric(sapply(data_list, function(x) x[5])),
          strand = sapply(data_list, function(x) x[6]),
          score = as.numeric(sapply(data_list, function(x) x[7])),
          length = as.numeric(sapply(data_list, function(x) x[8])),
          stringsAsFactors = FALSE
        )
        
        # Filter for Complete BUSCO genes only
        busco_df <- busco_df[busco_df$status == "Complete", ]
        
        return(busco_df)
      }
    }
    
    # Fallback: try simple tab-separated reading
    cat("Using fallback method for", basename(file_path), "\n")
    return(data.frame())
  }
  
  # Load BUSCO data for both species
  bibio_busco <- read_busco_file('busco-data/Bibio_marci.tsv')
  dilophus_busco <- read_busco_file('busco-data/Dilophus_febrilis.tsv')
  
  # Load network data
  nodes <- read.csv('multi_species_analysis/visualizations/nodes.csv')
  edges <- read.csv('multi_species_analysis/visualizations/edges.csv')
  
  cat(sprintf("Loaded: %d Bibio genes, %d Dilophus genes\n", 
              nrow(bibio_busco), nrow(dilophus_busco)))
  
  list(
    rearrangements = rearr_df,
    bibio_busco = bibio_busco,
    dilophus_busco = dilophus_busco,
    nodes = nodes,
    edges = edges
  )
}

# Simplified synteny block creation
create_synteny_blocks <- function(data) {
  cat("Creating synteny blocks from BUSCO orthologs...\n")
  
  if(nrow(data$bibio_busco) == 0 || nrow(data$dilophus_busco) == 0) {
    cat("ERROR: No BUSCO data loaded!\n")
    return(data.frame())
  }
  
  # Find common BUSCO genes
  common_buscos <- inner_join(
    data$bibio_busco, 
    data$dilophus_busco, 
    by = "busco_id",
    suffix = c("_bibio", "_dilophus")
  )
  
  cat(sprintf("Found %d common BUSCO genes\n", nrow(common_buscos)))
  
  if(nrow(common_buscos) == 0) {
    cat("ERROR: No common BUSCO genes found!\n")
    return(data.frame())
  }
  
  # Print sample of common genes
  cat("Sample common genes:\n")
  print(head(common_buscos[,c("busco_id", "sequence_bibio", "sequence_dilophus")], 3))
  
  # Create synteny blocks
  synteny_blocks <- common_buscos %>%
    # Remove genes with missing coordinates
    filter(!is.na(start_bibio), !is.na(start_dilophus)) %>%
    # Group by chromosome pairs that have multiple genes
    group_by(sequence_bibio, sequence_dilophus) %>%
    filter(n() >= 2) %>%  # At least 2 genes per block
    summarise(
      busco_count = n(),
      s1_start = min(start_bibio),
      s1_end = max(end_bibio),
      s2_start = min(start_dilophus),
      s2_end = max(end_dilophus),
      busco_genes = paste(busco_id[1:min(5, n())], collapse = ","),
      .groups = 'drop'
    ) %>%
    mutate(
      block_id = row_number(),
      s1_chr = paste0("Chr", dense_rank(sequence_bibio)),
      s2_chr = paste0("Chr", dense_rank(sequence_dilophus)),
      block_size = s1_end - s1_start,
      color = rainbow(n())[block_id]
    )
  
  cat(sprintf("Created %d synteny blocks\n", nrow(synteny_blocks)))
  
  return(synteny_blocks)
}

# Test data loading
test_data_loading <- function() {
  cat("=== Testing Data Loading ===\n")
  
  # Test BUSCO file reading
  cat("Testing Bibio BUSCO file...\n")
  lines <- readLines('busco-data/Bibio_marci.tsv')
  cat("First few lines:\n")
  cat(paste(head(lines, 5), collapse="\n"), "\n\n")
  
  # Find header
  header_line <- grep("^# Busco id", lines)
  if(length(header_line) > 0) {
    cat("Header found at line:", header_line, "\n")
    cat("Header:", lines[header_line], "\n")
    
    # Show first data line
    if(length(lines) > header_line) {
      cat("First data line:", lines[header_line + 1], "\n")
    }
  } else {
    cat("No header found!\n")
  }
}

# Run test first
test_data_loading()

# Then try to load data
cat("\n=== Attempting to Load Data ===\n")
data <- load_syngraph_data()

if(nrow(data$bibio_busco) > 0 && nrow(data$dilophus_busco) > 0) {
  cat("✓ Data loaded successfully!\n")
  synteny_blocks <- create_synteny_blocks(data)
  
  if(nrow(synteny_blocks) > 0) {
    cat("✓ Synteny blocks created successfully!\n")
    cat("Ready to create visualization!\n")
  } else {
    cat("✗ No synteny blocks created\n")
  }
} else {
  cat("✗ Failed to load BUSCO data\n")
}




#extract synteny block successfully, but does not create visualization.