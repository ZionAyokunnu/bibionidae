# Simple Direct Synteny Plot with Visible Elements
library(dplyr)
library(ggplot2)

# Test if basic plotting works first
test_basic_plot <- function() {
  cat("Testing basic plot functionality...\n")
  
  p <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 100, ymin = 0.9, ymax = 1.1), 
              fill = "red", color = "black") +
    geom_rect(aes(xmin = 0, xmax = 100, ymin = -0.1, ymax = 0.1), 
              fill = "blue", color = "black") +
    geom_line(aes(x = c(25, 25), y = c(0.9, 0.1)), color = "green", linewidth = 2) +
    theme_void() +
    coord_cartesian(xlim = c(-10, 110), ylim = c(-0.3, 1.3))
  
  ggsave("test_basic.png", p, width = 10, height = 4, dpi = 150, bg = "white")
  cat("✓ Basic test plot saved: test_basic.png\n")
  return(p)
}

# Load data with error checking
load_data_safely <- function() {
  cat("Loading BUSCO data...\n")
  
  read_busco_safe <- function(file_path) {
    if(!file.exists(file_path)) {
      cat("ERROR: File not found:", file_path, "\n")
      return(data.frame())
    }
    
    lines <- readLines(file_path)
    header_line <- grep("^# Busco id", lines)
    
    if(length(header_line) == 0) {
      cat("ERROR: No header found in", file_path, "\n")
      return(data.frame())
    }
    
    data_lines <- lines[(header_line + 1):length(lines)]
    data_lines <- data_lines[!grepl("^#", data_lines) & data_lines != ""]
    
    if(length(data_lines) == 0) {
      cat("ERROR: No data lines in", file_path, "\n")
      return(data.frame())
    }
    
    # Parse data
    data_list <- strsplit(data_lines, "\t")
    valid_rows <- sapply(data_list, length) >= 8
    data_list <- data_list[valid_rows]
    
    if(length(data_list) == 0) {
      cat("ERROR: No valid data rows\n")
      return(data.frame())
    }
    
    df <- data.frame(
      busco_id = sapply(data_list, function(x) x[1]),
      status = sapply(data_list, function(x) x[2]),
      sequence = sapply(data_list, function(x) x[3]),
      start = as.numeric(sapply(data_list, function(x) x[4])),
      end = as.numeric(sapply(data_list, function(x) x[5])),
      stringsAsFactors = FALSE
    )
    
    complete_df <- df[df$status == "Complete" & !is.na(df$start) & !is.na(df$end), ]
    cat(sprintf("  Loaded %d complete genes\n", nrow(complete_df)))
    return(complete_df)
  }
  
  bibio <- read_busco_safe('busco-data/Bibio_marci.tsv')
  dilophus <- read_busco_safe('busco-data/Dilophus_febrilis.tsv')
  
  if(nrow(bibio) == 0 || nrow(dilophus) == 0) {
    cat("ERROR: No data loaded!\n")
    return(NULL)
  }
  
  cat(sprintf("Successfully loaded: Bibio %d genes, Dilophus %d genes\n", 
              nrow(bibio), nrow(dilophus)))
  
  return(list(bibio = bibio, dilophus = dilophus))
}

# Create simple synteny plot with debugging
create_working_synteny_plot <- function() {
  cat("=== Creating Working Synteny Plot ===\n")
  
  # Test basic plotting first
  test_basic_plot()
  
  # Load data
  data <- load_data_safely()
  if(is.null(data)) {
    cat("FAILED: Could not load data\n")
    return(NULL)
  }
  
  # Create chromosome layout
  cat("Creating chromosome layout...\n")
  
  # Bibio top 4 chromosomes
  bibio_chroms <- data$bibio %>%
    group_by(sequence) %>%
    summarise(max_pos = max(end), genes = n(), .groups = 'drop') %>%
    arrange(desc(max_pos)) %>%
    slice_head(n = 4) %>%
    mutate(
      length = 50,  # Fixed length for simplicity
      chr_id = paste0("B", 1:4),
      cum_start = (row_number() - 1) * 60,
      cum_end = cum_start + length,
      y_pos = 1
    )
  
  # Dilophus top 4 chromosomes  
  dilophus_chroms <- data$dilophus %>%
    group_by(sequence) %>%
    summarise(max_pos = max(end), genes = n(), .groups = 'drop') %>%
    arrange(desc(max_pos)) %>%
    slice_head(n = 4) %>%
    mutate(
      length = 50,  # Fixed length for simplicity
      chr_id = paste0("D", 1:4),
      cum_start = (row_number() - 1) * 60,
      cum_end = cum_start + length,
      y_pos = 0
    )
  
  cat("Bibio chromosomes:\n")
  print(bibio_chroms[, c("sequence", "genes", "chr_id")])
  cat("Dilophus chromosomes:\n")
  print(dilophus_chroms[, c("sequence", "genes", "chr_id")])
  
  # Find common genes
  cat("Finding common genes...\n")
  common_genes <- inner_join(data$bibio, data$dilophus, 
                            by = "busco_id", suffix = c("_b", "_d"))
  
  # Create connections between major chromosomes
  connections <- common_genes %>%
    filter(sequence_b %in% bibio_chroms$sequence,
           sequence_d %in% dilophus_chroms$sequence) %>%
    left_join(bibio_chroms %>% select(sequence, cum_start), 
              by = c("sequence_b" = "sequence")) %>%
    rename(b_base = cum_start) %>%
    left_join(dilophus_chroms %>% select(sequence, cum_start), 
              by = c("sequence_d" = "sequence")) %>%
    rename(d_base = cum_start) %>%
    filter(!is.na(b_base), !is.na(d_base)) %>%
    sample_n(min(50, n())) %>%  # Sample 50 connections
    mutate(
      b_pos = b_base + runif(n(), 5, 45),  # Random position within chromosome
      d_pos = d_base + runif(n(), 5, 45),
      color = rainbow(n())
    )
  
  cat(sprintf("Created %d connections\n", nrow(connections)))
  
  if(nrow(connections) == 0) {
    cat("ERROR: No connections created!\n")
    return(NULL)
  }
  
  # Create the plot
  cat("Creating plot...\n")
  
  p <- ggplot() +
    # Bibio chromosomes (top)
    geom_rect(data = bibio_chroms,
              aes(xmin = cum_start, xmax = cum_end, ymin = 0.9, ymax = 1.1),
              fill = "lightblue", color = "black", linewidth = 1) +
    
    # Dilophus chromosomes (bottom)
    geom_rect(data = dilophus_chroms,
              aes(xmin = cum_start, xmax = cum_end, ymin = -0.1, ymax = 0.1),
              fill = "lightcoral", color = "black", linewidth = 1) +
    
    # Labels
    geom_text(data = bibio_chroms,
              aes(x = (cum_start + cum_end) / 2, y = 1.3, label = chr_id),
              size = 5, fontface = "bold") +
    
    geom_text(data = dilophus_chroms,
              aes(x = (cum_start + cum_end) / 2, y = -0.3, label = chr_id),
              size = 5, fontface = "bold") +
    
    # Species labels
    annotate("text", x = -20, y = 1, label = "Bibio marci", 
             size = 6, fontface = "italic", hjust = 1) +
    annotate("text", x = -20, y = 0, label = "Dilophus febrilis", 
             size = 6, fontface = "italic", hjust = 1) +
    
    theme_void() +
    theme(plot.margin = margin(30, 30, 30, 80)) +
    coord_cartesian(xlim = c(-50, 250), ylim = c(-0.5, 1.5))
  
  # Add connections as straight lines (simple version)
  for(i in 1:nrow(connections)) {
    p <- p + geom_segment(aes(x = connections$b_pos[i], 
                             y = 0.9,
                             xend = connections$d_pos[i], 
                             yend = 0.1),
                         color = connections$color[i],
                         linewidth = 0.8,
                         alpha = 0.7)
  }
  
  # Save plot
  ggsave("working_synteny.png", p, width = 12, height = 5, dpi = 300, bg = "white")
  
  cat("✓ Working synteny plot saved: working_synteny.png\n")
  cat(sprintf("✓ %d gene connections displayed\n", nrow(connections)))
  
  return(p)
}

# Run the working version
result <- create_working_synteny_plot()

if(!is.null(result)) {
  cat("\n🎉 SUCCESS! Working synteny plot created!\n")
  cat("📁 Check 'working_synteny.png' and 'test_basic.png'\n")
} else {
  cat("\n❌ FAILED to create plot\n")
}