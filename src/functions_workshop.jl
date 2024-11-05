using Turing, DataFrames, Distributions, Random

# File paths for input
times_path = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/samfire/Times.in"
haps_path = "/Users/e.smith.5/Documents/PhD/RSV_project/RSV_haplotype_reconstruction/data/VeTrans/by_protein/F/F_prev_in/Inference_29_3456.out"

# Read and create haps_wide DataFrame
times = readlines(times_path)
numeric_column_headers = parse.(Int, times)
haps_lines = readlines(haps_path)[2:end]
haps_data = [split(line) for line in haps_lines]
haplotypes = [row[1] for row in haps_data]
numeric_values = [parse.(Float64, row[2:end]) for row in haps_data]
haps_wide = DataFrame(Haplotype = haplotypes)
for (i, header) in enumerate(numeric_column_headers)
    haps_wide[!, Symbol(header)] = [row[i] for row in numeric_values]
end

# Define population evolution function with broadcasting
function population_evolution(q::Vector{Float64}, s::Vector{Float64}, steps::Int)
    q_values = [copy(q)]
    for t in 1:steps
        denominator = sum(s .* q)
        q_new = (s .* q) ./ denominator  # Broadcast division
        push!(q_values, q_new)
        q = q_new
    end
    return q_values
end

# Bayesian MCMC model definition
@model function bayesian_model(q_initial::Vector{Float64}, steps::Int)
    # Priors for selection coefficients `s_i`
    s ~ filldist(Uniform(0.0, 1.0), length(q_initial))

    # Predict frequencies across time points
    q_values = population_evolution(q_initial, s, steps)

    # Likelihood for observed data
    for t in 1:steps
        for i in 1:length(q_initial)
            # Ensure q_initial is Float64
            q_initial[i] ~ Normal(q_values[t][i], 0.1)
        end
    end
end

# Run MCMC for each row in haps_wide
steps = ncol(haps_wide) - 1
chains = []
for row in eachrow(haps_wide)
    q_initial = collect([row[i] for i in 2:ncol(haps_wide)])  # Extract row values as Array
    model = bayesian_model(q_initial, steps)  # No need to declare `model` as local
    sampled_chain = sample(model, NUTS(), 1000)  # Use NUTS or MH as needed
    push!(chains, sampled_chain)
end

# Display the chains for each haplotype
for chain in chains
    display(chain)
end