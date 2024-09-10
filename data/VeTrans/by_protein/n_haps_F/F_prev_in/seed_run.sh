#!/bin/bash

# Initialize the log file with headers
echo "file,BIC" > seed_run.out

# Default values
traj_in="Multi_locus_trajectories.out"
n_samples=2
max_haps=20 # Default value

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --traj_in)
        traj_in="$2"
        shift
        ;;
        --n_samples)
        n_samples="$2"
        shift
        ;;
        --max_haps)
        max_haps="$2"
        shift
        ;;
        *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
    shift
done

# Initial run of the calculation
while read seed
do
    ./VeTrans reconstruct_haplotypes_multi --traj_in "$traj_in" --n_samples "$n_samples" --n_haps 4 --seed $seed &
done < seeds.in
wait

one=1
# Loop from 5 to the value of max_haps
for ((i=5; i<=max_haps; i++))
do
    j=$((i - one))
    
    # Extract all BIC values and their corresponding files
    for file in Inference_$j*
    do
        # Ensure each file is processed independently
        bic_value=$(grep "Final BIC" $file | awk '{print $3}')
        
        # Log the file name and BIC to the seed_run.out
        if [ ! -z "$bic_value" ]; then
            echo "$file,$bic_value" >> seed_run.out
        fi
    done
    
    # Extract the best inference file (highest BIC)
    best_inference=$(grep "Final BIC" Inference_$j* | sort -nrk3 | head -n1 | cut -f1 -d ":")
    echo $best_inference > best_inference.in
    
    while read best
    do
        ./cut_inference.sh $best
    done < best_inference.in
    
    while read seed
    do
        ./VeTrans reconstruct_haplotypes_multi --traj_in "$traj_in" --n_samples "$n_samples" --n_haps $i --seed $seed &
    done < seeds.in
    wait
done