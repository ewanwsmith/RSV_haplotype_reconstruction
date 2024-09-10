using DataFrames
using CSV

# Function to process Newick file using embedded R code
function newick_to_df(filepath::String)
    # Create the R command as a properly formatted string
    r_script = """
    library(ape)
    library(dplyr)

    # Read the Newick tree file
    tree <- ape::read.tree(file = '$filepath')

    # Extract the edges from the tree structure
    edges <- as.data.frame(tree\$edge)
    colnames(edges) <- c('From', 'To')

    # Get the node labels (sequence names)
    labels <- tree\$tip.label
    num_tips <- length(labels)
    all_labels <- c(labels, as.character((num_tips + 1):(num_tips + tree\$Nnode)))

    # Add the labels to the edges DataFrame
    edges\$FromLabel <- all_labels[edges\$From]
    edges\$ToLabel <- all_labels[edges\$To]

    # Print the DataFrame as CSV to standard output
    write.csv(edges, row.names = FALSE)
    """

    # Execute the R script via command line and capture the output
    cmd = `Rscript -e $r_script`
    csv_output = read(cmd, String)
    
    # Parse the CSV output into a DataFrame in Julia
    df = CSV.read(IOBuffer(csv_output), DataFrame)
    
    return df
end

# Path to the Newick tree file
treefile_path = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/iqtree/F_haplotypes/14/Consensus0_F_14_haplotypes_JC.treefile"

# Get the DataFrame
julia_df = newick_to_df(treefile_path)

# Display the DataFrame
println(julia_df)