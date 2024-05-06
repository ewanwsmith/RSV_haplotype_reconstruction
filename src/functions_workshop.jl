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
    
    # Filter rows where the maximum of A, C, G, T is less than n * proportion
    return filter(row -> max(row.A, row.C, row.G, row.T) < row.n * proportion, df)
end
