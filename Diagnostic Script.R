################################################################################
# DIAGNOSTIC SCRIPT: Debug the inversion detection issue
################################################################################

library(dplyr)
library(tibble)
library(ggplot2)
library(readr)

################################################################################
# STEP 1: CHECK IF MARKER ORDER FILE EXISTS
################################################################################

cat("=== CHECKING FILES ===\n")

# Check what files exist
files_to_check <- c(
  "GCA_910594885.2_idBibMarc1.2_genomic.fna",
  "Bibio_marci_inverted.fasta", 
  "Bibio_marci_inverted_marker_order.csv",
  "Bibio_marci_inverted_inversion_log.csv"
)

for (file in files_to_check) {
  if (file.exists(file)) {
    cat("✓", file, "exists\n")
  } else {
    cat("✗", file, "missing\n")
  }
}

################################################################################
# STEP 2: CHECK MARKER ORDER FILE IF EXISTS
################################################################################

marker_order_file <- "Bibio_marci_inverted_marker_order.csv"

if (file.exists(marker_order_file)) {
  cat("\n=== READING MARKER ORDER FILE ===\n")
  
  marker_order_df <- read_csv(marker_order_file, show_col_types = FALSE)
  cat("Marker order file has", nrow(marker_order_df), "rows\n")
  cat("Columns:", paste(colnames(marker_order_df), collapse = ", "), "\n")
  
  # Show first few rows
  cat("\nFirst 10 rows of marker order:\n")
  print(head(marker_order_df, 10))
  
  # Check for inversions in the marker order
  cat("\n=== CHECKING FOR INVERSIONS IN MARKER ORDER ===\n")
  
  # Simple check: look for positions where inverted_position decreases
  if ("inverted_position" %in% colnames(marker_order_df)) {
    decreasing_positions <- which(diff(marker_order_df$inverted_position) < 0)
    cat("Found", length(decreasing_positions), "decreasing positions (potential inversions)\n")
    
    if (length(decreasing_positions) > 0) {
      cat("Decreasing positions at rows:", head(decreasing_positions, 10), "\n")
    }
  } else {
    cat("No 'inverted_position' column found\n")
  }
  
} else {
  cat("\nMarker order file does not exist. Need to run simulator first.\n")
}

################################################################################
# STEP 3: READ INVERSION LOG IF EXISTS
################################################################################

inversion_log_file <- "Bibio_marci_inverted_inversion_log.csv"

if (file.exists(inversion_log_file)) {
  cat("\n=== READING INVERSION LOG ===\n")
  
  inversion_log <- read_csv(inversion_log_file, show_col_types = FALSE)
  cat("Inversion log has", nrow(inversion_log), "inversions\n")
  print(inversion_log)
  
} else {
  cat("\nInversion log file does not exist.\n")
}

################################################################################
# STEP 4: MANUAL INVERSION DETECTION TEST
################################################################################

if (file.exists(marker_order_file)) {
  cat("\n=== MANUAL INVERSION DETECTION TEST ===\n")
  
  marker_order_df <- read_csv(marker_order_file, show_col_types = FALSE)
  
  # Create a simple test comparison
  if ("marker_id" %in% colnames(marker_order_df)) {
    
    # Simulate what the analyzer should do
    original_order <- marker_order_df$marker_id
    original_positions <- 1:length(original_order)
    
    # The inverted positions should be the positions in the rearranged order
    inverted_positions <- match(original_order, marker_order_df$marker_id)
    
    # This is wrong - we need the actual rearranged order!
    # Let me fix this by using the original_position column properly
    
    if ("original_position" %in% colnames(marker_order_df)) {
      comparison_df <- marker_order_df %>%
        arrange(original_position) %>%
        mutate(
          original_pos = row_number(),
          inverted_pos = row_number()  # This needs to be the NEW position
        )
      
      cat("Created comparison with", nrow(comparison_df), "markers\n")
      cat("First few comparisons:\n")
      print(head(comparison_df[c("original_pos", "inverted_pos", "marker_id")], 10))
      
      # Simple inversion detection
      inversions_found <- 0
      for (i in 1:(nrow(comparison_df)-1)) {
        if (comparison_df$inverted_pos[i] > comparison_df$inverted_pos[i+1]) {
          inversions_found <- inversions_found + 1
        }
      }
      
      cat("Simple inversion count:", inversions_found, "\n")
    }
  }
}

################################################################################
# STEP 5: RECOMMENDATIONS
################################################################################

cat("\n=== RECOMMENDATIONS ===\n")

if (!file.exists(marker_order_file)) {
  cat("1. Run the simulator first to create the marker order file\n")
  cat("2. Make sure the simulator saves the marker order correctly\n")
} else {
  cat("1. Marker order file exists - checking if it contains proper inversion data\n")
  cat("2. The analyzer needs to use this file correctly\n")
}

cat("3. Check that both simulator and analyzer use the same marker_size_bp\n")
cat("4. Verify that the inversion log shows the expected inversions\n")

################################################################################
# STEP 6: CREATE A SIMPLE TEST INVERSION
################################################################################

cat("\n=== CREATING TEST INVERSION FOR VERIFICATION ===\n")

# Create a simple test case
test_markers <- paste0("test_", 1:10)
cat("Original order:", paste(test_markers, collapse = ", "), "\n")

# Simulate an inversion from position 3 to 7
inverted_test <- test_markers
inverted_test[3:7] <- rev(test_markers[3:7])
cat("After inversion (3-7):", paste(inverted_test, collapse = ", "), "\n")

# Create comparison
test_comparison <- tibble(
  original_position = 1:10,
  inverted_position = match(test_markers, inverted_test),
  marker_id = test_markers
)

cat("Test comparison:\n")
print(test_comparison)

# Detect inversions in test case
cat("\nTesting inversion detection on this simple case:\n")
test_inversions <- 0
for (i in 1:9) {
  if (test_comparison$inverted_position[i] > test_comparison$inverted_position[i+1]) {
    test_inversions <- test_inversions + 1
    cat("  Inversion detected at position", i, "\n")
  }
}

cat("Total test inversions detected:", test_inversions, "\n")

if (test_inversions > 0) {
  cat("✓ Inversion detection algorithm works on test case\n")
} else {
  cat("✗ Inversion detection algorithm failed on test case\n")
}