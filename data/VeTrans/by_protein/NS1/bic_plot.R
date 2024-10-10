library(ggplot2)

# Read the CSV file
bic <- read.csv("/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/VeTrans/by_protein/NS1/extracted_bic.csv")

# Remove rows where BIC is NA
bic <- bic[!is.na(bic), ]

# Ensure n_haps is numeric and sort by n_haps
bic$n_haps <- as.numeric(as.character(bic$n_haps))
bic <- bic[order(bic$n_haps), ]

# Create the plot
bic_plot <- ggplot(bic, aes(x = n_haps, y = BIC)) +
  geom_point() +
  labs(title="NS1", x="# Haplotypes", y="BIC") +
  theme_minimal()

# Save the plot
ggsave("/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/VeTrans/by_protein/NS1/bic_plot.jpeg", plot = bic_plot, width = 8, height = 6)

# Find the row with the lowest BIC
min_bic <- bic[which.min(bic$BIC), ]

# Print the file with the lowest BIC
cat("File with the lowest BIC:", min_bic$file, "\n")
