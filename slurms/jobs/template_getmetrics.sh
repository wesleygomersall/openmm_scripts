#!/bin/bash

# Template to initiate slurm jobs which use mdtraj to get metrics from MD
# trajectory pdb. 

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/getmetrics.sbatch'
T_NAME='output.pdb'
OUTPUTNAME='mdtrajanalysis.csv'

# Associative array for trajectories and their reference files
declare -A directories_references

# Add like so: (do not use the trajectory's path, input path of its directory)
directories_references[/path/to/trajectorydir/]='/path/to/reference.pdb'

# Example: 
# directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134502_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
# directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134526_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
