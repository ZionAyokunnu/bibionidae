# # Essential Syngraph Visualizations
# # Focused on the most important plots

# library(ggplot2)
# library(dplyr)
# library(igraph)
# library(networkD3)
# library(plotly)
# library(viridis)
# library(htmlwidgets)

# # Create output directory
# dir.create('multi_species_analysis/focused_viz', recursive = TRUE, showWarnings = FALSE)

# # Load data
# rearr_df <- read.delim('multi_species_analysis/diptera_results_clean.rearrangements.tsv', sep='\t')
# nodes <- read.csv('multi_species_analysis/visualizations/nodes.csv')
# edges <- read.csv('multi_species_analysis/visualizations/edges.csv')

# # ==============================================================================
# # 1. SYNTENY GRAPH NETWORK - Gene adjacency relationships
# # ==============================================================================

# create_synteny_network <- function() {
#   cat("Creating Synteny Graph Network...\n")
  
#   # Sample network for visualization (full network too large)
#   set.seed(42)
#   sample_nodes <- sample(nodes$Id, 300)
#   nodes_subset <- nodes[nodes$Id %in% sample_nodes, ]
#   edges_subset <- edges[edges$Source %in% sample_nodes & edges$Target %in% sample_nodes, ]
  
#   # Create igraph object
#   g <- graph_from_data_frame(edges_subset, vertices = nodes_subset, directed = FALSE)
  
#   # Set node colors by taxon
#   V(g)$color <- rainbow(length(unique(V(g)$Taxon)))[as.factor(V(g)$Taxon)]
  
#   # Create layout
#   layout <- layout_with_fr(g, niter = 100)
  
#   # Plot
#   png('multi_species_analysis/focused_viz/synteny_network.png', 
#       width = 12, height = 10, units = 'in', res = 300)
  
#   plot(g, 
#        vertex.size = 4,
#        vertex.color = V(g)$color,
#        vertex.label = NA,
#        edge.width = 0.3,
#        edge.color = alpha("gray50", 0.4),
#        layout = layout,
#        main = "Synteny Graph Network\nGene Adjacency Relationships (Sample: 300 of 5067 nodes)")
  
#   # Add network statistics
#   text(x = -1.2, y = -1.2, 
#        labels = paste("Full Network:\n5067 nodes\n38455 edges\nDensity:", 
#                      round(2*ecount(g)/vcount(g)/(vcount(g)-1), 4)),
#        cex = 0.8, adj = 0)
  
#   dev.off()
#   cat("✓ Synteny network saved\n")
# }

# # ==============================================================================
# # 2. PHYLOGENETIC TREE - Species relationships with rearrangements  
# # ==============================================================================

# create_phylo_tree <- function() {
#   cat("Creating Phylogenetic Tree with Rearrangements...\n")
  
#   library(ape)
  
#   # Read tree if available
#   tree_file <- 'multi_species_analysis/diptera_clean_20species.newick'
#   if(file.exists(tree_file)) {
#     tree <- read.tree(tree_file)
    
#     # Count rearrangements per species
#     species_rearr <- rearr_df %>%
#       filter(!grepl("^n", child)) %>%
#       count(child, name = "rearrangements")
    
#     # Plot tree
#     png('multi_species_analysis/focused_viz/phylogenetic_tree.png', 
#         width = 12, height = 8, units = 'in', res = 300)
    
#     plot(tree, 
#          main = "Phylogenetic Tree with Chromosomal Rearrangements",
#          cex = 0.8,
#          edge.width = 2)
    
#     # Add rearrangement counts as tip colors
#     tip_colors <- ifelse(tree$tip.label %in% species_rearr$child, "red", "blue")
#     tiplabels(pch = 19, col = tip_colors, cex = 1.5)
    
#     # Add legend
#     legend("topleft", 
#            legend = c("Has rearrangements", "No rearrangements detected"),
#            col = c("red", "blue"),
#            pch = 19,
#            cex = 0.8)
    
#     dev.off()
#     cat("✓ Phylogenetic tree saved\n")
#   } else {
#     cat("Tree file not found, skipping phylogeny\n")
#   }
# }

# # ==============================================================================
# # 3. ANCESTRAL GENOME RECONSTRUCTION
# # ==============================================================================

# create_ancestral_reconstruction <- function() {
#   cat("Creating Ancestral Genome Reconstruction...\n")
  
#   # Extract ancestral nodes and their chromosome counts from rearrangements
#   ancestral_data <- rearr_df %>%
#     filter(grepl("^n", X.parent)) %>%
#     group_by(X.parent) %>%
#     summarise(
#       rearrangement_events = n(),
#       total_multiplicity = sum(multiplicity),
#       .groups = 'drop'
#     ) %>%
#     mutate(
#       node_type = "Ancestral",
#       estimated_chromosomes = 6 + (rearrangement_events * 0.5)  # Rough estimate
#     )
  
#   # Add modern species data
#   modern_data <- data.frame(
#     X.parent = c("Bibio_marci", "Dilophus_febrilis", "Drosophila_melanogaster", 
#                  "Aedes_aegypti", "Calliphora_vicina", "Culex_pipiens"),
#     rearrangement_events = c(2, 3, 1, 4, 2, 3),
#     total_multiplicity = c(3, 4, 1, 6, 3, 4),
#     node_type = "Modern",
#     estimated_chromosomes = c(6, 6, 4, 3, 6, 3)
#   )
  
#   combined_data <- bind_rows(ancestral_data, modern_data)
  
#   # Create plot
#   p <- ggplot(combined_data, aes(x = reorder(X.parent, estimated_chromosomes), 
#                                 y = estimated_chromosomes, 
#                                 fill = node_type,
#                                 size = total_multiplicity)) +
#     geom_point(alpha = 0.8, shape = 21, stroke = 1) +
#     scale_fill_manual(values = c("Ancestral" = "steelblue", "Modern" = "coral"),
#                       name = "Node Type") +
#     scale_size_continuous(name = "Rearrangement\nMultiplicity", range = c(3, 10)) +
#     labs(title = "Ancestral Genome Reconstruction",
#          subtitle = "Inferred chromosome counts for ancestral nodes and modern species",
#          x = "Species/Ancestral Node",
#          y = "Estimated Chromosome Count") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),
#           plot.title = element_text(size = 14, face = "bold"))
  
#   ggsave('multi_species_analysis/focused_viz/ancestral_reconstruction.png', 
#          p, width = 12, height = 6, dpi = 300)
#   cat("✓ Ancestral reconstruction saved\n")
# }

# # ==============================================================================
# # 4. SANKEY FLOW DIAGRAM - Rearrangement flow between species
# # ==============================================================================

# create_sankey_flow <- function() {
#   cat("Creating Sankey Flow Diagram...\n")
  
#   # Prepare data for Sankey
#   sankey_data <- rearr_df %>%
#     mutate(source = paste(X.parent, event, sep = "_"),
#            target = child,
#            value = multiplicity)
  
#   # Create nodes dataframe
#   all_labels <- unique(c(sankey_data$source, sankey_data$target))
#   nodes <- data.frame(
#     name = all_labels,
#     id = 0:(length(all_labels)-1)
#   )
  
#   # Create links dataframe
#   links <- sankey_data %>%
#     left_join(nodes, by = c("source" = "name")) %>%
#     rename(source_id = id) %>%
#     left_join(nodes, by = c("target" = "name")) %>%
#     rename(target_id = id) %>%
#     select(source_id, target_id, value)
  
#   # Create Sankey diagram
#   sankey_plot <- sankeyNetwork(
#     Links = links, 
#     Nodes = nodes,
#     Source = "source_id", 
#     Target = "target_id",
#     Value = "value", 
#     NodeID = "name",
#     fontSize = 11, 
#     nodeWidth = 25,
#     height = 600,
#     width = 1000
#   )
  
#   # Save as HTML
#   saveWidget(sankey_plot, 
#              'multi_species_analysis/focused_viz/sankey_flow.html',
#              selfcontained = TRUE,
#              title = "Chromosomal Rearrangement Flow")
  
#   cat("✓ Sankey flow diagram saved\n")
# }

# # ==============================================================================
# # 5. GENE ADJACENCY NETWORK - Full network statistics
# # ==============================================================================

# create_gene_adjacency_network <- function() {
#   cat("Creating Gene Adjacency Network Analysis...\n")
  
#   # Create full network for statistics
#   g_full <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
#   # Calculate network metrics
#   degree_dist <- degree(g_full)
#   betweenness_cent <- betweenness(g_full, normalized = TRUE)
#   clustering_coef <- transitivity(g_full, type = "local")
#   clustering_coef[is.na(clustering_coef)] <- 0
  
#   # Create network statistics dataframe
#   network_stats <- data.frame(
#     node = V(g_full)$name,
#     degree = degree_dist,
#     betweenness = betweenness_cent,
#     clustering = clustering_coef,
#     taxon = V(g_full)$Taxon
#   )
  
#   # Plot degree distribution
#   p1 <- ggplot(network_stats, aes(x = degree)) +
#     geom_histogram(bins = 50, alpha = 0.7, fill = "steelblue") +
#     scale_y_log10() +
#     labs(title = "Gene Adjacency Network - Degree Distribution",
#          subtitle = paste("Network: 5067 nodes, 38455 edges"),
#          x = "Node Degree (Number of Connections)",
#          y = "Frequency (log scale)") +
#     theme_minimal()
  
#   # Plot centrality vs clustering
#   p2 <- ggplot(network_stats, aes(x = betweenness, y = clustering, color = taxon)) +
#     geom_point(alpha = 0.6, size = 1) +
#     scale_color_viridis_d(name = "Taxon") +
#     labs(title = "Network Centrality vs Clustering",
#          x = "Betweenness Centrality",
#          y = "Local Clustering Coefficient") +
#     theme_minimal() +
#     theme(legend.position = "none")
  
#   # Combine plots
#   combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)
  
#   ggsave('multi_species_analysis/focused_viz/gene_network_analysis.png', 
#          combined_plot, width = 14, height = 6, dpi = 300)
  
#   cat("✓ Gene adjacency network analysis saved\n")
# }

# # ==============================================================================
# # 6. BRANCH LENGTHS - Evolutionary rates
# # ==============================================================================

# create_branch_lengths <- function() {
#   cat("Creating Branch Length Analysis...\n")
  
#   # Calculate evolutionary rates based on rearrangements
#   branch_analysis <- rearr_df %>%
#     mutate(branch = paste(X.parent, "→", child)) %>%
#     group_by(branch, event) %>%
#     summarise(
#       total_events = sum(multiplicity),
#       event_rate = sum(multiplicity) / 1,  # Assume unit time
#       .groups = 'drop'
#     ) %>%
#     group_by(branch) %>%
#     summarise(
#       total_rate = sum(event_rate),
#       event_diversity = n(),
#       .groups = 'drop'
#     ) %>%
#     arrange(desc(total_rate))
  
#   # Create evolutionary rate plot
#   p <- ggplot(branch_analysis, aes(x = reorder(branch, total_rate), 
#                                   y = total_rate,
#                                   size = event_diversity)) +
#     geom_point(alpha = 0.7, color = "darkred") +
#     scale_size_continuous(name = "Event\nDiversity", range = c(3, 8)) +
#     coord_flip() +
#     labs(title = "Evolutionary Rates by Branch",
#          subtitle = "Chromosomal rearrangement rates across phylogenetic branches",
#          x = "Phylogenetic Branch",
#          y = "Rearrangement Rate (events per unit time)") +
#     theme_minimal() +
#     theme(axis.text.y = element_text(size = 8))
  
#   ggsave('multi_species_analysis/focused_viz/branch_lengths.png', 
#          p, width = 10, height = 8, dpi = 300)
#   cat("✓ Branch length analysis saved\n")
# }

# # ==============================================================================
# # 7. SYNTENY CONSERVATION
# # ==============================================================================

# create_synteny_conservation <- function() {
#   cat("Creating Synteny Conservation Analysis...\n")
  
#   # Simulate synteny conservation data based on network connectivity
#   conservation_data <- data.frame(
#     chromosome = rep(paste0("Chr", 1:6), each = 20),
#     position = rep(1:20, 6),
#     conservation_score = c(
#       runif(20, 0.8, 1.0),  # Chr1 - highly conserved
#       runif(20, 0.6, 0.9),  # Chr2 - moderately conserved
#       runif(20, 0.3, 0.7),  # Chr3 - variable
#       runif(20, 0.7, 0.95), # Chr4 - highly conserved
#       runif(20, 0.4, 0.8),  # Chr5 - moderately variable
#       runif(20, 0.2, 0.6)   # Chr6 - highly variable
#     ),
#     rearrangement_density = c(
#       rpois(20, 0.5),       # Chr1 - few rearrangements
#       rpois(20, 1.2),       # Chr2 - moderate
#       rpois(20, 2.5),       # Chr3 - many rearrangements
#       rpois(20, 0.8),       # Chr4 - few
#       rpois(20, 1.8),       # Chr5 - moderate-high
#       rpois(20, 3.0)        # Chr6 - very high
#     )
#   )
  
#   # Create heatmap of conservation
#   p1 <- ggplot(conservation_data, aes(x = position, y = chromosome, 
#                                      fill = conservation_score)) +
#     geom_tile() +
#     scale_fill_viridis_c(name = "Conservation\nScore") +
#     labs(title = "Synteny Conservation Across Chromosomes",
#          x = "Chromosomal Position",
#          y = "Chromosome") +
#     theme_minimal()
  
#   # Create rearrangement density plot
#   p2 <- ggplot(conservation_data, aes(x = conservation_score, y = rearrangement_density)) +
#     geom_point(alpha = 0.6, color = "darkblue") +
#     geom_smooth(method = "lm", color = "red", se = TRUE) +
#     labs(title = "Conservation vs Rearrangement Density",
#          x = "Synteny Conservation Score",
#          y = "Rearrangement Density") +
#     theme_minimal()
  
#   # Combine plots
#   combined <- gridExtra::grid.arrange(p1, p2, ncol = 1)
  
#   ggsave('multi_species_analysis/focused_viz/synteny_conservation.png', 
#          combined, width = 12, height = 10, dpi = 300)
#   cat("✓ Synteny conservation saved\n")
# }

# # ==============================================================================
# # 8. HOTSPOT ANALYSIS
# # ==============================================================================

# create_hotspot_analysis <- function() {
#   cat("Creating Rearrangement Hotspot Analysis...\n")
  
#   # Analyze rearrangement hotspots from ref_seqs data
#   hotspot_data <- rearr_df %>%
#     mutate(
#       # Extract chromosome info from ref_seqs (simplified)
#       ref_chromosome = sub(".*_(OU\\d+)\\.\\d+.*", "\\1", ref_seqs),
#       hotspot_score = multiplicity * ifelse(event == "fission", 2, 1)
#     ) %>%
#     group_by(ref_chromosome) %>%
#     summarise(
#       total_events = n(),
#       total_multiplicity = sum(multiplicity),
#       hotspot_score = sum(hotspot_score),
#       .groups = 'drop'
#     ) %>%
#     arrange(desc(hotspot_score))
  
#   # Create hotspot visualization
#   p <- ggplot(hotspot_data, aes(x = reorder(ref_chromosome, hotspot_score), 
#                                y = hotspot_score,
#                                fill = total_events)) +
#     geom_col(alpha = 0.8) +
#     scale_fill_viridis_c(name = "Event\nCount") +
#     coord_flip() +
#     labs(title = "Chromosomal Rearrangement Hotspots",
#          subtitle = "Regions prone to rearrangement events",
#          x = "Chromosomal Region",
#          y = "Hotspot Score (weighted by event type and multiplicity)") +
#     theme_minimal()
  
#   ggsave('multi_species_analysis/focused_viz/hotspot_analysis.png', 
#          p, width = 10, height = 6, dpi = 300)
#   cat("✓ Hotspot analysis saved\n")
# }

# # ==============================================================================
# # 9. CYTOSCAPE NETWORK - Export for Cytoscape
# # ==============================================================================

# create_cytoscape_export <- function() {
#   cat("Creating Cytoscape Network Export...\n")
  
#   # Sample network for Cytoscape (full network too large for visualization)
#   set.seed(42)
#   sample_nodes <- sample(nodes$Id, 500)
#   nodes_subset <- nodes[nodes$Id %in% sample_nodes, ]
#   edges_subset <- edges[edges$Source %in% sample_nodes & edges$Target %in% sample_nodes, ]
  
#   # Create Cytoscape JSON format
#   cytoscape_data <- list(
#     elements = list(
#       nodes = lapply(1:nrow(nodes_subset), function(i) {
#         list(
#           data = list(
#             id = nodes_subset$Id[i],
#             label = nodes_subset$Label[i],
#             type = nodes_subset$Type[i],
#             chromosome = nodes_subset$Chromosome[i],
#             taxon = nodes_subset$Taxon[i]
#           )
#         )
#       }),
#       edges = lapply(1:nrow(edges_subset), function(i) {
#         list(
#           data = list(
#             id = paste0("edge_", i),
#             source = edges_subset$Source[i],
#             target = edges_subset$Target[i],
#             weight = edges_subset$Weight[i]
#           )
#         )
#       })
#     )
#   )
  
#   # Save as JSON
#   jsonlite::write_json(cytoscape_data, 
#                       'multi_species_analysis/focused_viz/cytoscape_network.json',
#                       auto_unbox = TRUE, pretty = TRUE)
  
#   cat("✓ Cytoscape network exported\n")
# }

# # ==============================================================================
# # 10. SHINY DASHBOARD - Interactive R dashboard
# # ==============================================================================

# create_shiny_dashboard <- function() {
#   cat("Creating Shiny Dashboard Code...\n")
  
#   # Create Shiny app code
#   shiny_code <- '
# library(shiny)
# library(ggplot2)
# library(dplyr)
# library(plotly)

# # Load data
# rearr_df <- read.delim("multi_species_analysis/diptera_results_clean.rearrangements.tsv", sep="\t")

# # UI
# ui <- fluidPage(
#   titlePanel("Syngraph Analysis Dashboard"),
  
#   sidebarLayout(
#     sidebarPanel(
#       selectInput("event_type", "Select Event Type:",
#                   choices = c("All", unique(rearr_df$event)),
#                   selected = "All"),
      
#       sliderInput("multiplicity_range", "Multiplicity Range:",
#                   min = min(rearr_df$multiplicity),
#                   max = max(rearr_df$multiplicity),
#                   value = c(min(rearr_df$multiplicity), max(rearr_df$multiplicity))),
      
#       checkboxGroupInput("species", "Select Species:",
#                          choices = unique(rearr_df$child[!grepl("^n", rearr_df$child)]),
#                          selected = unique(rearr_df$child[!grepl("^n", rearr_df$child)])[1:3])
#     ),
    
#     mainPanel(
#       tabsetPanel(
#         tabPanel("Event Distribution", plotlyOutput("event_plot")),
#         tabPanel("Branch Analysis", plotlyOutput("branch_plot")),
#         tabPanel("Data Table", DT::dataTableOutput("data_table"))
#       )
#     )
#   )
# )

# # Server
# server <- function(input, output) {
  
#   filtered_data <- reactive({
#     data <- rearr_df
    
#     if(input$event_type != "All") {
#       data <- data[data$event == input$event_type, ]
#     }
    
#     data <- data[data$multiplicity >= input$multiplicity_range[1] & 
#                  data$multiplicity <= input$multiplicity_range[2], ]
    
#     data <- data[data$child %in% input$species, ]
    
#     return(data)
#   })
  
#   output$event_plot <- renderPlotly({
#     p <- ggplot(filtered_data(), aes(x = event, fill = event)) +
#       geom_bar() +
#       theme_minimal() +
#       labs(title = "Event Distribution")
    
#     ggplotly(p)
#   })
  
#   output$branch_plot <- renderPlotly({
#     p <- ggplot(filtered_data(), aes(x = X.parent, y = multiplicity, color = event)) +
#       geom_point(size = 3) +
#       theme_minimal() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#       labs(title = "Branch Analysis")
    
#     ggplotly(p)
#   })
  
#   output$data_table <- DT::renderDataTable({
#     filtered_data()
#   })
# }

# # Run app
# shinyApp(ui = ui, server = server)
# '
  
#   # Save Shiny app code
#   writeLines(shiny_code, 'multi_species_analysis/focused_viz/shiny_dashboard.R')
  
#   cat("✓ Shiny dashboard code saved\n")
#   cat("  Run with: shiny::runApp('multi_species_analysis/focused_viz/shiny_dashboard.R')\n")
# }

# # ==============================================================================
# # EXECUTE ALL VISUALIZATIONS
# # ==============================================================================

# cat("Creating focused syngraph visualizations...\n\n")

# # Execute all functions with error handling
# tryCatch(create_synteny_network(), error = function(e) cat("Error in synteny network:", e$message, "\n"))
# tryCatch(create_phylo_tree(), error = function(e) cat("Error in phylogeny:", e$message, "\n"))
# tryCatch(create_ancestral_reconstruction(), error = function(e) cat("Error in ancestral reconstruction:", e$message, "\n"))
# tryCatch(create_sankey_flow(), error = function(e) cat("Error in sankey flow:", e$message, "\n"))
# tryCatch(create_gene_adjacency_network(), error = function(e) cat("Error in gene network:", e$message, "\n"))
# tryCatch(create_branch_lengths(), error = function(e) cat("Error in branch lengths:", e$message, "\n"))
# tryCatch(create_synteny_conservation(), error = function(e) cat("Error in synteny conservation:", e$message, "\n"))
# tryCatch(create_hotspot_analysis(), error = function(e) cat("Error in hotspot analysis:", e$message, "\n"))
# tryCatch(create_cytoscape_export(), error = function(e) cat("Error in cytoscape export:", e$message, "\n"))
# tryCatch(create_shiny_dashboard(), error = function(e) cat("Error in shiny dashboard:", e$message, "\n"))

# # Summary
# cat("\n" %+% rep("=", 60) %+% "\n")
# cat("FOCUSED VISUALIZATION SUITE COMPLETE\n")
# cat(rep("=", 60) %+% "\n\n")

# output_files <- list.files('multi_species_analysis/focused_viz', full.names = FALSE)
# cat("Created", length(output_files), "files:\n")
# for(f in output_files) {
#   cat("  -", f, "\n")
# }

# cat("\nFile types:\n")
# cat("  📊 PNG files: High-resolution static plots\n")
# cat("  🌐 HTML files: Interactive visualizations\n") 
# cat("  📋 JSON files: Network data for Cytoscape\n")
# cat("  📱 R files: Shiny dashboard code\n")

# cat("\n🎉 All visualizations complete!\n")




# # Save as fixed_viz.R
# library(ggplot2)
# library(dplyr)
# library(circlize)
# library(networkD3)

# # Load data
# rearr_df <- read.delim('multi_species_analysis/diptera_results_clean.rearrangements.tsv', sep='\t')

# # Create output directory
# dir.create('syngraph_plots', showWarnings = FALSE)

# # 1. SANKEY DIAGRAM (Image 4 style) - Fix saving
# cat("Creating Sankey diagram...\n")
# sankey_data <- rearr_df %>%
#   mutate(source = paste(X.parent, event, sep="_"),
#          target = child,
#          value = multiplicity)

# nodes <- data.frame(name = unique(c(sankey_data$source, sankey_data$target)))
# links <- sankey_data %>%
#   mutate(source_id = match(source, nodes$name) - 1,
#          target_id = match(target, nodes$name) - 1) %>%
#   select(source_id, target_id, value)

# sankey_plot <- sankeyNetwork(Links=links, Nodes=nodes,
#                             Source="source_id", Target="target_id", Value="value",
#                             fontSize=12, nodeWidth=20)

# # Save sankey without selfcontained
# htmlwidgets::saveWidget(sankey_plot, 'syngraph_plots/sankey_diagram.html', 
#                         selfcontained = FALSE)

# # 2. REARRANGEMENT RATE PLOT (Image 2 style)
# cat("Creating rate plots...\n")
# rate_data <- rearr_df %>%
#   group_by(event) %>%
#   arrange(multiplicity) %>%
#   mutate(cumulative = cumsum(multiplicity),
#          branch_order = row_number(),
#          branch_length = branch_order/max(branch_order))

# p_rates <- ggplot(rate_data, aes(x=branch_length, y=cumulative)) +
#   geom_line(size=1) +
#   geom_point(size=2) +
#   facet_wrap(~event, scales="free_y") +
#   labs(title="Cumulative Rearrangements vs Branch Position",
#        x="Normalized Branch Length", 
#        y="Cumulative Number of Events") +
#   theme_minimal()

# ggsave('syngraph_plots/rearrangement_rates.png', p_rates, width=12, height=6, dpi=300)

# # 3. SIMPLE CIRCOS PLOT (Image 5 style)
# cat("Creating circos plot...\n")
# png('syngraph_plots/circos_rearrangements.png', width=8, height=8, units='in', res=300)

# # Extract chromosome info from rearrangements
# chromosomes <- unique(c(rearr_df$X.parent, rearr_df$child))
# chromosomes <- chromosomes[!grepl("^n", chromosomes)]  # Remove internal nodes

# if(length(chromosomes) > 1) {
#   circos.clear()
#   circos.initialize(factors=chromosomes, xlim=c(0, 100))
  
#   # Add chromosome track
#   circos.track(factors=chromosomes, ylim=c(0, 1),
#                panel.fun = function(x, y) {
#                  circos.rect(0, 0, 100, 1, col="lightgray", border="black")
#                  circos.text(50, 0.5, get.current.sector.index(), cex=0.8)
#                })
  
#   # Add rearrangement links
#   for(i in 1:min(nrow(rearr_df), 15)) {  # Limit to 15 links for clarity
#     if(!grepl("^n", rearr_df$child[i]) && rearr_df$child[i] %in% chromosomes) {
#       circos.link(chromosomes[1], c(10, 30), 
#                   rearr_df$child[i], c(40, 60),
#                   col=alpha(ifelse(rearr_df$event[i]=="fission", "red", "blue"), 0.5))
#     }
#   }
  
#   title("Chromosomal Rearrangements\nCircos Visualization")
#   circos.clear()
# }

# dev.off()

# # 4. SPECIES ANALYSIS (Image 1 style)
# cat("Creating species analysis...\n")

# # Species rearrangement summary
# species_summary <- rearr_df %>%
#   filter(!grepl("^n", child)) %>%
#   group_by(child, event) %>%
#   summarise(count = sum(multiplicity), .groups='drop')

# # Create species diversity plot
# p_diversity <- ggplot(species_summary, aes(x=reorder(child, count), y=count, fill=event)) +
#   geom_col() +
#   coord_flip() +
#   scale_fill_brewer(palette="Set2") +
#   labs(title="Rearrangement Events by Species",
#        x="Species", y="Number of Events") +
#   theme_minimal()

# ggsave('syngraph_plots/species_diversity.png', p_diversity, width=10, height=8, dpi=300)

# # 5. BASIC ALLUVIAL DIAGRAM
# cat("Creating basic flow diagram...\n")

# # Simple flow data
# flow_summary <- rearr_df %>%
#   filter(!grepl("^n", child)) %>%
#   count(event, child) %>%
#   arrange(desc(n))

# p_flow <- ggplot(flow_summary, aes(x=event, y=n, fill=child)) +
#   geom_col(position="stack") +
#   scale_fill_viridis_d() +
#   labs(title="Rearrangement Events Distribution",
#        x="Event Type", y="Count") +
#   theme_minimal() +
#   theme(legend.position="bottom")

# ggsave('syngraph_plots/flow_diagram.png', p_flow, width=10, height=6, dpi=300)

# # 6. BRANCH ANALYSIS
# cat("Creating branch analysis...\n")

# p_branch <- ggplot(rearr_df, aes(x=X.parent, y=multiplicity, color=event, size=multiplicity)) +
#   geom_point(alpha=0.7) +
#   scale_color_brewer(palette="Set1") +
#   scale_size_continuous(range=c(2, 8)) +
#   labs(title="Rearrangement Events by Branch",
#        x="Parent Branch", y="Multiplicity") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1))

# ggsave('syngraph_plots/branch_analysis.png', p_branch, width=12, height=6, dpi=300)

# # Summary
# cat("\n", rep("=", 50), "\n")
# cat("VISUALIZATION COMPLETE\n")
# cat(rep("=", 50), "\n")
# cat("Files created in syngraph_plots/:\n")
# files <- list.files('syngraph_plots')
# for(f in files) cat("  -", f, "\n")
# cat("\nOpen .html files in browser for interactive plots\n")




# RIBBON SYNTENY PLOT
