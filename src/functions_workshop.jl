# try and pull two .fas from .dat file 
using BioSequences
using FASTX

function pull_frames(dat_file_path::String)
    # Open the .dat file to read
    file = open(dat_file_path, "r")
    fasta_nucleotides_path = replace(dat_file_path, ".dat" => "_nucleotide.fa")
    fasta_amino_acids_path = replace(dat_file_path, ".dat" => "_amino_acid.fa")
    fasta_nucleotides = FASTA.Writer(open(fasta_nucleotides_path, "w"))
    fasta_amino_acids = FASTA.Writer(open(fasta_amino_acids_path, "w"))

    try
        while !eof(file)
            line = readline(file)
            if isempty(line) || occursin(r"^\s*$", line)  # Skip empty lines or lines with only whitespace
                continue
            end
            fields = split(line)
            if length(fields) < 3
                continue  # Skip malformed header lines
            end
            annotation = join(fields[3:end], " ")
            nucleotide_sequence = readline(file)
            amino_acid_sequence = readline(file)

            # Convert string to DNA sequence
            dna_seq = LongDNASeq(nucleotide_sequence)

            # Writing to nucleotide fasta file
            write(fasta_nucleotides, FASTA.Record(annotation, dna_seq))
            # Writing to amino acid fasta file
            write(fasta_amino_acids, FASTA.Record(annotation, amino_acid_sequence))
        end
    finally
        close(fasta_nucleotides)
        close(fasta_amino_acids)
        close(file)
    end
end




pull_frames("/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/ref/Reading_Frames.dat")