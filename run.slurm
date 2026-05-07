#!/bin/bash
#SBATCH --job-name=qtt_purification
#SBATCH --output=qtt_run_%j.out
#SBATCH --error=qtt_run_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=cfp
#SBATCH --account=cfp
#SBATCH --nodelist=cfp11
#SBATCH --gres=gpu:1

# Ensure output directory exists


# Execute the Julia script. 
# --project=. ensures it uses the Project.toml in the current directory.
julia --project=. tests/test_timing_2D.jl