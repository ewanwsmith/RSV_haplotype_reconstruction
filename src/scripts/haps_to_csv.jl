#!/usr/bin/env julia

using Pkg

# Install required packages if they are not already installed
function install_packages()
    packages = ["ArgParse", "CSV", "DataFrames", "FilePathsBase"]
    installed_pkgs = Pkg.dependencies()

    # Check if each package is installed and install if it is not
    for pkg in packages
        if !any(p -> p.name == pkg, values(installed_pkgs))
            Pkg.add(pkg)
        end
    end
end

install_packages()

using ArgParse
using CSV
using DataFrames
using FilePathsBase  # To handle file paths easily

# Function to parse the .out file
function parse_out_file(out_path::String)
    haplotypes = []
    frequencies = []

    open(out_path, "r") do file
        for line in eachline(file)
            # Split the line by a regex that captures the boundary between a number and a letter
            haplotype_freq_str = split(line, r"(?<=[0-9e\.-])(?=[A-Z])")
            for haplotype_freq in haplotype_freq_str
                # Extract the haplotype sequence
                haplotype = match(r"[A-Z]+", haplotype_freq).match
                # Replace the haplotype part with an empty string and then split the remaining frequencies
                freqs_str = replace(haplotype_freq, haplotype => "")
                freqs = parse.(Float64, filter(x -> !isempty(x), split(freqs_str, r"\s+")))
                push!(haplotypes, haplotype)
                push!(frequencies, freqs)
            end
        end
    end

    return haplotypes, frequencies
end

# Function to parse the Times.in file
function parse_times_file(times_path::String)
    timepoints = []

    open(times_path, "r") do file
        for line in eachline(file)
            push!(timepoints, parse(Int, line))
        end
    end

    return timepoints
end

# Function to convert data to CSV
function convert_to_csv(out_path::String, times_path::String, csv_path::String)
    haplotypes, frequencies = parse_out_file(out_path)
    timepoints = parse_times_file(times_path)
    
    if length(timepoints) != length(frequencies[1])
        throw(ArgumentError("Number of timepoints does not match the number of frequency columns"))
    end

    df = DataFrame()
    df[!, "Haplotype"] = haplotypes
    for (index, tp) in enumerate(timepoints)
        df[!, string(tp)] = [freq[index] for freq in frequencies]
    end

    CSV.write(csv_path, df)
end

# Main function to handle command-line arguments and call the processing functions
function main()
    # Define argument parser
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--in"
            help = "Path to the .out file"
            arg_type = String
            required = true
        "--times"
            help = "Path to the Times.in file"
            arg_type = String
            required = true
        "--out"
            help = "Path to the output CSV file"
            arg_type = String
            required = false
    end

    # Parse command-line arguments
    parsed_args = parse_args(s)

    # Extract arguments
    in_path = parsed_args["in"]
    times_path = parsed_args["times"]
    csv_path = get(parsed_args, "out", nothing)

    # Set default output path if not provided
    if isnothing(csv_path)
        base_dir = dirname(in_path)
        csv_filename = replace(basename(in_path), ".out" => ".csv")
        csv_path = joinpath(base_dir, csv_filename)
    end

    # Convert to CSV
    convert_to_csv(in_path, times_path, csv_path)
end

# Call the main function
main()