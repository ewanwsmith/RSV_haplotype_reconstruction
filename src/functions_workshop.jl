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


#look for variants at >10% frequency
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
