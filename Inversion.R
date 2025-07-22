# Extract Bibio_marci genome markers in linear order
bibio_df <- df_lin %>%                   # Start from your full dataframe of linearized genes
  filter(species == "Bibio_marci") %>%   # Keep only rows where the species is Bibio_marci
  arrange(linear_start)                  # Sort genes in genome order by their linear_start position

# Get the gene order as a vector (actual gene IDs from Bibio_marci)
bibio_genome <- bibio_df$marker          # Extract the marker column as a vector of gene names/IDs

# Function to perform a single inversion on a gene vector
invert_segment <- function(genome_vec) {
  n <- length(genome_vec)                           # Get genome length
  i <- sample(1:(n - 1), 1)                          # Randomly pick a start index for inversion
  j <- sample((i + 1):n, 1)                          # Randomly pick an end index greater than i
  genome_vec[i:j] <- rev(genome_vec[i:j])           # Reverse the gene segment between i and j
  return(list(new_genome = genome_vec,              # Return the modified genome
              inversion = c(i, j)))                 # Log the coordinates of the inversion
}

# Function to simulate multiple inversion events
simulate_inversions <- function(genome_vec, num_inversions = 10) {
  genome_states <- list()                                        # Will store genome after each inversion
  inversions_log <- matrix(NA, nrow = num_inversions, ncol = 2)  # Pre-allocate inversion log
  colnames(inversions_log) <- c("start_idx", "end_idx")          # Name the log columns
  
  for (k in 1:num_inversions) {
    result <- invert_segment(genome_vec)         # Run one inversion
    genome_vec <- result$new_genome              # Update genome to reflect that inversion
    inversions_log[k, ] <- result$inversion      # Record the start/end indices of that inversion
    genome_states[[k]] <- genome_vec             # Save the state of genome after this inversion
  }
  
  return(list(final_genome = genome_vec,         # Return final genome after all inversions
              inversions = inversions_log,       # Full log of inversion positions
              history = genome_states))          # All intermediate genome states
}

# Simulate 30 inversions on Bibio_marci genome
sim_result <- simulate_inversions(bibio_genome, num_inversions = 1)  # Run simulation with 30 events

# Print the inversion start and end positions (1-based gene indices)
print(sim_result$inversions)

# Construct a dataframe for plotting the original vs final gene positions
df_plot <- data.frame(
  original_index = 1:length(bibio_genome),                            # Original positions 1:N
  final_index = match(bibio_genome, sim_result$final_genome),         # Where each gene ended up in final genome
  gene = bibio_genome                                                 # Gene names (optional use in labels)
)

# Generate the dotplot comparing original to final gene order
p <- ggplot(df_plot, aes(x = original_index, y = final_index)) +
  geom_point(size = 1.5, color = "darkred") +                         # Red dots represent gene shifts
  theme_minimal() +                                                   # Minimal theme for clean plot
  labs(
    title = "Manual Inversion Simulation on Bibio_marci",             # Plot title
    x = "Original Genome Index",                                      # X-axis label
    y = "Post-Inversion Index"                                        # Y-axis label
  )
print(p)                                                              # Explicitly print the plot in script

# Inversion frequency analysis — which genes were affected most
inverted_genes <- c()                                         # Init empty list of affected genes
for (k in 1:nrow(sim_result$inversions)) {
  i <- sim_result$inversions[k, 1]                            # Get start index of inversion
  j <- sim_result$inversions[k, 2]                            # Get end index of inversion
  inverted_genes <- c(inverted_genes, bibio_genome[i:j])     # Append all gene names in the inverted region
}

# Count how many times each gene was involved in any inversion
gene_freq <- as.data.frame(table(inverted_genes)) %>%         # Tabulate frequencies of gene appearance
  arrange(desc(Freq))                                         # Sort from most to least affected genes

# View top 10 most frequently inverted genes
head(gene_freq, 10)
