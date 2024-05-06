# dependencies
using DelimitedFiles
using DataFrames

# try to read in .out files to a dataframe in as few steps as possible
function read_out_file(filename::AbstractString; delim::Char=' ')
    data = readdlm(filename, delim, Int)
    col_names = [:pos, :A, :C, :G, :T, :n]
    df = DataFrame(data, col_names)
    return df
end


# function to filter out positions with low coverage
function filter_coverage(df::DataFrame, threshold::Int)
    return filter(row -> row.n >= threshold, df)
end


# look for variants at >10% frequency
function filter_positions(df::DataFrame, proportion::Float64)
    # Validate that the required columns are present in the DataFrame
    required_columns = ["A", "C", "G", "T", "n"]
    for col in required_columns
        if !(col in names(df))
            error("DataFrame must contain the columns: A, C, G, T, and n")
        end
    end
    
    # Function to check if any of three smallest columns are greater than n * proportion
    function check_proportion(row)
        values = [row.A, row.C, row.G, row.T]
        sorted_indices = sortperm(values, rev=true)  # Sort indices based on values descending
        max_index = sorted_indices[1]                # Index of the maximum value
        remaining_values = deleteat!(values, max_index)  # Remove the maximum value
        
        # Check if any remaining value exceeds n * proportion
        return any(v -> v > row.n * proportion, remaining_values)
    end
    
    # Filter rows based on the check_proportion function
    return filter(check_proportion, df)
end

# function to pull out the variants in an easier format
function pull_variants(df::DataFrame)
    # Define the result DataFrame with the specified columns
    result_df = DataFrame(pos = Int[], original_base = String[], variant_base = String[], freq = Float64[])
    
    # Process each row of the input DataFrame
    for row in eachrow(df)
        pos = row[:pos]
        bases = [:A, :C, :G, :T]
        counts = [row.A, row.C, row.G, row.T]
        
        # Sort bases by their counts in descending order
        sorted_indices = sortperm(counts, rev=true)
        
        # Find the bases with the highest and second highest counts
        highest_base = bases[sorted_indices[1]]
        second_highest_base = bases[sorted_indices[2]]
        
        # Calculate the frequency of the second highest base
        freq = counts[sorted_indices[2]] / row.n
        
        # Append the results to the DataFrame, converting symbols to strings
        push!(result_df, (pos, string(highest_base), string(second_highest_base), freq))
    end
    
    return result_df
end

