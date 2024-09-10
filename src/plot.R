library(ggtree)
library(ggridges)
library(aplot)
library(ggplot2)
library(dplyr)
library(purrr)

# Generate a minimal example phylogenetic tree
tree <- rtree(3)  # A small tree with 3 tips

# Plot the phylogenetic tree
p_tree <- ggtree(tree) + theme_tree()

# Generate minimal example data for the ridgeline plot
set.seed(123)
tip_labels <- tree$tip.label
data <- data.frame(
  label = rep(tip_labels, each = 50),
  value = unlist(lapply(tip_labels, function(x) rnorm(50, mean = which(tip_labels == x))))
)

# Compute density for each group
density_data <- data %>%
  group_by(label) %>%
  do({
    dens <- density(.$value)
    data.frame(x = dens$x, y = dens$y)
  }) %>%
  ungroup() %>%
  mutate(y = y * 1000)  # Rescale y for better visibility

# Create the ridgeline plot
p_ridges <- ggplot(density_data, aes(x = x, y = label, height = y, fill = label)) +
  geom_density_ridges(stat = "identity", scale = 1, alpha = 0.8) +
  theme_ridges() +
  theme(legend.position = "none")

# Combine the plots using aplot's insert_right
aligned_plot <- insert_right(p_tree, p_ridges, width = 1.5)

# Display the combined plot
print(aligned_plot)