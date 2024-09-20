library(csv)
library(ggplot2)

# Read the CSV file
bic <- read.csv("/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/VeTrans/by_protein/G/G_prev_in/extracted_bic.csv")

# Print the bic object to check its contents
print(bic)

# Remove rows where BIC is NA
bic <- bic[!is.na(bic), ]

# Ensure that n_haps is treated as numeric and sort by n_haps
bic$n_haps <- as.numeric(as.character(bic$n_haps))
bic <- bic[order(bic$n_haps), ]

# Create the plot with points only (no lines) and dynamic axis and title
bic_plot <- ggplot(bic, aes(x = n_haps, y = BIC)) +
  geom_point() +
  labs(title="G_prev_in", x="# Haplotypes", y="BIC") +
  theme_minimal()

# Define the path for saving the plot
output_file <- "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/VeTrans/by_protein/G/G_prev_in/bic_plot.jpeg"

# Save the plot as a .jpeg file
ggsave(output_file, plot = bic_plot, width = 8, height = 6)

# No need to plot the graph in the R script
