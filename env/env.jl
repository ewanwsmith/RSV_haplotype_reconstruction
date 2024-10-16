# Julia setup
using Pkg, Logging

# Set the Julia environment path
Julia_env_path = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/env"

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

# load R dependencies
R"""
# Function to set CRAN mirror and ensure R can download packages
chooseCRANmirror(ind=1) # Selects a default CRAN mirror

# Install the devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# Load the devtools package
library(devtools)

# Install ggplot2 and viridis with dependencies
install.packages(c("ggplot2", "viridis", "tidyverse", "hrbrthemes", "plotly", "htmlwidgets", "ggridges", "aplot"), dependencies = TRUE)

# Function to check and install missing packages
install_if_missing <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        BiocManager::install(package)
    } else {
        message(paste(package, "is already installed."))
    }
}

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install ggtree and treeio if not already installed
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
# Output setup complete message after the process
display("Julia and R setup complete")
