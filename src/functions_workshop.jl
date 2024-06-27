using RCall
using DataFrames
using CSV

# Function to extract data from a .fasta file
function pull_fasta(fasta_path::String)
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
    accessions = rcopy(R"accession")
    starts = rcopy(R"start")
    ends = rcopy(R"end")
    proteins = rcopy(R"protein")
    sequences = rcopy(R"sapply(dna_seqs, as.character)")

    # Create a DataFrame from the extracted data
    fasta_df = DataFrame(
        protein = proteins,
        start_pos = starts,
        end_pos = ends,
        sequence = sequences
    )

    return fasta_df
end

# Function to read .out file and split it based on protein information
function split_out_file(out_file::String, protein_info::DataFrame)
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

    for row in eachrow(protein_info)
        protein_name = replace(row[:protein], " " => "_")
        start_pos = row[:start_pos]
        end_pos = row[:end_pos]

        # Filter rows where locus is between start and stop positions
        filtered_data = filter(row -> row[:locus] >= start_pos && row[:locus] <= end_pos, out_data)

        if nrow(filtered_data) > 0
            # Create new .out file for the protein
            new_out_file = replace(out_file, ".out" => "_$protein_name.out")
            open(new_out_file, "w") do io
                for row in eachrow(filtered_data)
                    println(io, "$(row.locus) $(row.remainder)")
                end
            end
            println("Written file: $new_out_file")
        end
    end
end

# Main function
function process_files(fasta_file::String, out_file::String)
    # Extract protein information from the fasta file
    protein_info = pull_fasta(fasta_file)
    
    # Split the .out file based on protein information
    split_out_file(out_file, protein_info)
end

# Example usage
fasta_file = "data/ref/Reading_Frames_nucleotide.fa"
out_file = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/samfire/single_locus_trajectories10.out"
process_files(fasta_file, out_file)