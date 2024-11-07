# Julia setup
using Pkg
using Logging

# Set the Julia environment path
Julia_env_path = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/env/"

# Create the folder if it doesn't exist using built-in Julia functions
if !isdir(Julia_env_path)
    mkpath(Julia_env_path)
end

# Activate the Julia environment
Pkg.activate(Julia_env_path)
Pkg.instantiate()

using IJulia
using CSV
using DataFrames
using RCall
using Images
using FileIO
using FilePathsBase
using WebIO
using CategoricalArrays
using Distributions
using Optim
using Roots

# Load R dependencies
R"""
# Set CRAN mirror quietly
chooseCRANmirror(ind=1)

# Install necessary packages quietly
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools", quiet = TRUE)

library(devtools, quietly = TRUE)

# Install packages quietly with dependencies
install.packages(c("ggplot2", "viridis", "tidyverse", "hrbrthemes", "plotly", "htmlwidgets", "ggridges", "aplot"), dependencies = TRUE, quiet = TRUE)

# Function to check and install missing packages quietly
install_if_missing <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        BiocManager::install(package, quiet = TRUE)
    }
}

# Ensure BiocManager is installed quietly
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
}

# Install additional packages quietly
install_if_missing("ggtree")
install_if_missing("treeio")
install_if_missing("ggtreeExtra")

# Load the packages quietly
suppressPackageStartupMessages({
  library("ggplot2")
  library("viridis")
  library("tidyverse")
  library("hrbrthemes")
  library("plotly")
  library("htmlwidgets")
  library("ggridges")
  library("ggtree")
  library("treeio")
  library("ggtreeExtra")
  library("aplot")
})
"""

# Display setup complete message
display("Julia and R setup complete")