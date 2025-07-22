# Create better visualizations with R
library(ggplot2)
library(dplyr)
library(plotly)
library(networkD3)
library(treemap)
library(RColorBrewer)

# Load data
df <- read.delim("Documents/Bibionidae/multi_species_analysis/diptera_results_clean.rearrangements.tsv", sep="\t")


# 1. Professional bar chart
p1 <- ggplot(df, aes(x=event, fill=event)) +
  geom_bar(alpha=0.8) +
  scale_fill_brewer(palette="Set2") +
  labs(title="Chromosomal Rearrangement Events", 
       x="Event Type", y="Count") +
  theme_minimal() +
  theme(legend.position="none")

ggsave('Documents/Bibionidae/multi_species_analysis/visualizations/events_ggplot.png', p1, width=8, height=6, dpi=300)

# 2. Treemap of events by branch
treemap_data <- df %>%
  count(X.parent, child, event) %>%
  mutate(branch = paste(X.parent, "→", child))

treemap(treemap_data, 
        index=c("branch", "event"),
        vSize="n",
        type="index",
        palette="Set3",
        title="Rearrangements by Branch")

# 3. Sankey diagram with networkD3
sankey_data <- df %>%
  mutate(source = paste(X.parent, event, sep="_"),
         target = child) %>%
  select(source, target, multiplicity)

# Create nodes
nodes <- data.frame(
  name = c(unique(sankey_data$source), unique(sankey_data$target))
)

# Create links
links <- sankey_data %>%
  mutate(source = match(source, nodes$name) - 1,
         target = match(target, nodes$name) - 1,
         value = multiplicity)

# Create Sankey
sankeyNetwork(Links = links, Nodes = nodes, 
              Source = "source", Target = "target", 
              Value = "value", NodeID = "name",
              fontSize = 12, nodeWidth = 30)