#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=rmsf
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1                  # Number of tasks per node
#SBATCH --cpus-per-task=8           # Number of CPU cores per task
#SBATCH --time=7-00:00:00           # Time limit 7 days
#SBATCH --mem=32G                   # Memory per node

# Wesley Gomersall
# Thu, 21 Aug 2025 13:50:45 -0700
# Calculate RMSF from trajectories

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

# Associative array for trajectories and their reference files
declare -A directories_references

# Add like so: (do not use the trajectory's path, input path of its directory)
directories_references[/path/to/trajectorydir/]='/path/to/reference.pdb'

# Example: 
# directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134502_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
# directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134526_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

for key in "${!directories_references[@]}"; do 
    cd $key
    python3 ~/openmm_scripts/mdtraj_metrics/rmsf.py -i $key -r ${directories_references[$key]}
done
