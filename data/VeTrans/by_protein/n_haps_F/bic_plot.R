library(csv)
library(ggplot2)

# Read the CSV file
bic = read.csv("/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/VeTrans/by_protein/n_haps_F/extracted_bic.csv")

# Create the plot
bic_plot = ggplot(bic, aes(x = Run, y = BIC)) +
  geom_point() +
  geom_line() +
  labs(title="BIC for n_haps", x="n_haps", y="BIC") +
  theme_minimal()

# Define the path for saving the plot
output_file = "data/VeTrans/by_protein/n_haps_F/bic_plot.jpeg"

# Save the plot as a .jpeg file
ggsave(output_file, plot = bic_plot, width = 8, height = 6)

# No need to plot the graph in the R script