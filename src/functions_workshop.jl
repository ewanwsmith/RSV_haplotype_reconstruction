using RCall
using DataFrames
using CSV

# Function to read the .fasta file and extract the required information using R's Bioconductor package
function read_fasta(fasta_file::String)
    R"""
    library(Biostrings)
    fasta <- readAAStringSet($fasta_file)
    fasta_df <- data.frame(
        protein_names = names(fasta),
        sequences = as.character(fasta)
    )
    fasta_df$start <- as.numeric(gsub(".*:(\\d+)-\\d+.*", "\\1", fasta_df$protein_names))
    fasta_df$stop <- as.numeric(gsub(".*:\\d+-(\\d+).*", "\\1", fasta_df$protein_names))
    fasta_df$protein_names <- gsub(" .*", "", fasta_df$protein_names)
    fasta_df
    """ |> DataFrame
end

# Function to read the .out file and filter non-synonymous mutations
function read_out(out_file::String)
    out_data = CSV.read(out_file, DataFrame; delim=' ')
    
    # Prepare data for R function
    out_df = DataFrame(
        locus = out_data.Column1,
        original_base = out_data.Column2,
        variant_base = out_data.Column3,
        additional_info = join.(eachrow(out_data[:, 4:end]), ' ')
    )
    
    # R script to filter non-synonymous mutations
    non_synonymous_mutations = R"""
    library(Biostrings)
    
    filter_non_synonymous <- function(df) {
        codon_table <- getGeneticCode("Standard")
        
        # Function to get the codon at a given position
        get_codon <- function(sequence, position) {
            start <- position - ((position - 1) %% 3)
            end <- start + 2
            return(substr(sequence, start, end))
        }
        
        # Add columns for codons and amino acids
        df$original_codon <- mapply(get_codon, as.character(df$original_base), df$locus)
        df$variant_codon <- mapply(get_codon, as.character(df$variant_base), df$locus)
        df$original_aa <- translate(DNAStringSet(df$original_codon))
        df$variant_aa <- translate(DNAStringSet(df$variant_codon))
        df$non_synonymous <- df$original_aa != df$variant_aa
        
        return(df[df$non_synonymous, ])
    }
    
    out_df <- as.data.frame($out_df)
    filtered_df <- filter_non_synonymous(out_df)
    """ |> DataFrame
    
    return non_synonymous_mutations(out_df)
end

# Function to compare fasta and out data
function compare_fasta_out(fasta_df::DataFrame, out_df::DataFrame)
    # Placeholder example of comparison logic
    joined_df = leftjoin(out_df, fasta_df, on = :locus => :start)
    return joined_df
end

# Main function to run the process and save the new .out file
function process_files(fasta_file::String, out_file::String)
    fasta_df = read_fasta(fasta_file)
    out_df = read_out(out_file)
    result_df = compare_fasta_out(fasta_df, out_df)
    
    # Save the filtered .out file with _ns suffix
    output_file = replace(out_file, ".out" => "_ns.out")
    CSV.write(output_file, result_df)
    
    return result_df
end

# Example usage
fasta_file = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/ref/Reading_Frames_nucleotide.fa"
out_file = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/samfire/single_locus_trajectories10.out"
result_df = process_files(fasta_file, out_file)
println(result_df)