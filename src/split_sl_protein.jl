#!/usr/bin/env julia

using Pkg

# Function to check if a package is installed
function is_installed(pkg::String)
    installed_pkgs = Pkg.installed()
    return haskey(installed_pkgs, pkg)
end

# Install required packages if not already installed
packages = ["RCall", "DataFrames", "CSV", "ArgParse"]
for pkg in packages
    if !is_installed(pkg)
        println("Installing $pkg...")
        Pkg.add(pkg)
    else
        println("$pkg is already installed.")
    end
end

# Load the required packages
using RCall
using DataFrames
using CSV
using ArgParse

# Ensure the required R libraries are installed
println("Checking and installing required R packages...")
R"""
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("Biostrings")
}
"""

# Function to extract data from a .fasta file
function pull_fasta(fasta_path::String)
    println("Reading and processing the fasta file...")
    # Read the fasta file using Biostrings in R
    R"""
    library(Biostrings)
    library(stringr)
    dna_seqs <- readDNAStringSet($fasta_path)
    headers <- names(dna_seqs)
    
    # Parse headers to extract Accession, Start, End, and Protein
    accession <- str_extract(headers, "^[^:]+")
    range <- str_extract(headers, "\\d+\\.\\.\\d+")
    start <- as.integer(sub("\\.\\..*", "", range))
    end <- as.integer(sub("^.*\\.\\.", "", range))
    protein <- str_extract(headers, "\\|\\s*(.*)$")
    protein <- str_trim(str_sub(protein, str_locate(protein, "\\|")[[1]] + 1)) # Trim pipe and whitespace
    """

    # Convert R data to Julia DataFrame
    accessions = try rcopy(R"accession") catch _ missing end
    starts = try rcopy(R"start") catch _ missing end
    ends = try rcopy(R"end") catch _ missing end
    proteins = try rcopy(R"protein") catch _ missing end
    sequences = try rcopy(R"sapply(dna_seqs, as.character)") catch _ missing end

    # Handle missing values by replacing with empty arrays if necessary
    accessions = ismissing(accessions) ? String[] : accessions
    starts = ismissing(starts) ? Int[] : starts
    ends = ismissing(ends) ? Int[] : ends
    proteins = ismissing(proteins) ? String[] : proteins
    sequences = ismissing(sequences) ? String[] : sequences

    # Create a DataFrame from the extracted data
    fasta_df = DataFrame(
        protein = proteins,
        start_pos = starts,
        end_pos = ends,
        sequence = sequences
    )

    println("Fasta file processed successfully.")
    println("Extracted data from fasta file:")
    println(fasta_df)
    return fasta_df
end

# Function to read .out file and split it based on protein information
function split_out_file(out_file::String, protein_info::DataFrame, out_dir::String)
    println("Reading the .out file...")
    # Read .out file line by line
    lines = readlines(out_file)
    
    # Initialize DataFrame to store locus and remainder
    out_data = DataFrame(locus = Int[], remainder = String[])
    
    for line in lines
        # Split the line into components
        parts = split(line)
        # Extract the first integer as locus
        locus = parse(Int, parts[1])
        # Join the rest as the remainder
        remainder = join(parts[2:end], " ")
        # Push to DataFrame
        push!(out_data, (locus, remainder))
    end

    println(".out file read successfully. Processing proteins...")
    println("Data from .out file:")
    println(out_data)
    
    for row in eachrow(protein_info)
        protein_name = replace(row[:protein], " " => "_")
        start_pos = row[:start_pos]
        end_pos = row[:end_pos]

        # Filter rows where locus is between start and stop positions
        filtered_data = filter(row -> row[:locus] >= start_pos && row[:locus] <= end_pos, out_data)

        if nrow(filtered_data) > 0
            # Create new .out file for the protein in the specified output directory
            new_out_file = joinpath(out_dir, replace(basename(out_file), ".out" => "_$protein_name.out"))
            open(new_out_file, "w") do io
                for row in eachrow(filtered_data)
                    println(io, "$(row.locus) $(row.remainder)")
                end
            end
            println("Written file: $new_out_file")
        else
            println("No data for protein $protein_name in range $start_pos - $end_pos")
        end
    end
    println("All proteins processed successfully.")
end

# Main function
function process_files(fasta_file::String, out_file::String, out_dir::String)
    println("Starting file processing...")
    # Extract protein information from the fasta file
    protein_info = pull_fasta(fasta_file)
    
    # Split the .out file based on protein information
    split_out_file(out_file, protein_info, out_dir)
    println("File processing completed.")
end

# Parse command line arguments
function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--ref"
        help = "Path to the reference fasta file"
        arg_type = String
    end

    @add_arg_table s begin
        "--in"
        help = "Path to the input .out file"
        arg_type = String
    end

    @add_arg_table s begin
        "--out"
        help = "Path to the output directory"
        arg_type = String
        default = ""
    end

    parsed_args = parse_args(s)

    fasta_file = parsed_args["ref"]
    out_file = parsed_args["in"]
    out_dir = parsed_args["out"] == "" ? dirname(out_file) : parsed_args["out"]

    process_files(fasta_file, out_file, out_dir)
end

main()