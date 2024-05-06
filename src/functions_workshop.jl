# try to read in .out files to a dataframe in as few steps as possible
using DelimitedFiles
using DataFrames

function read_out_file(filename::AbstractString; delim::Char=' ')
    data = readdlm(filename, delim, Int)
    col_names = [:pos, :A, :C, :G, :T, :n]
    df = DataFrame(data, col_names)
    return df
end
