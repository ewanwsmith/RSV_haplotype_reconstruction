using Pkg

# Ensure necessary packages are installed
packages = ["BioSequences", "ArgParse"]
for pkg in packages
    if !haskey(Pkg.installed(), pkg)
        Pkg.add(pkg)
    end
end

using BioSequences
using ArgParse

function parse_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input"
            help = "Input FASTA file"
            required = true
        "--output"
            help = "Output FASTA file"
            required = true
        "--start"
            help = "Start position (1-based index)"
            required = true
            arg_type = Int
        "--end"
            help = "End position (1-based index)"
            required = true
            arg_type = Int
    end
    return parse_args(s)
end

function trim_sequences(input_file::String, output_file::String, start::Int, stop::Int)
    open(output_file, "w") do out
        for record in open(BioSequence.FASTA.Reader, input_file)
            seq = BioSequence.sequence(record)
            if start > 0 && stop <= length(seq) && start <= stop
                trimmed_seq = subseq(seq, start:stop)
                write(out, '>', BioSequence.identifier(record), '\n')
                write(out, string(trimmed_seq), '\n')
            else
                println("Warning: Skipping sequence ", BioSequence.identifier(record), " due to invalid range.")
            end
        end
    end
end

function main()
    args = parse_args()
    input_file = args["input"]
    output_file = args["output"]
    start = args["start"]
    stop = args["end"]

    trim_sequences(input_file, output_file, start, stop)
end

main()