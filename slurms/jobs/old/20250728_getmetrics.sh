#!/bin/bash

# Template to initiate slurm jobs which use mdtraj to get metrics from MD
# trajectory pdb. 

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/getmetrics.sbatch'
T_NAME='output.pdb'
OUTPUTNAME='mdtrajanalysis.csv'

# Associative array for trajectories and their reference files
declare -A directories_references

# Add like so: (do not use the trajectory's path, input path of its directory)
# directories_references[/path/to/trajectorydir/]='/path/to/reference.pdb'

# files
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T154838_c16test_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T154838_c16test_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'

directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T150204_ala_newforcepos1_1_NP]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T150204_ala_oldforce-1_2_NP]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T150204_ala_oldforce-1_1_NP]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T134653_ala_oldforce-1_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T154838_alatest_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T134653_ala_oldforce-1_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T154838_alatest_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T134653_ala_newforcepos1_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T150204_ala_newforcepos1_2_NP]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250728_SMD/20250728T134653_ala_newforcepos1_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
