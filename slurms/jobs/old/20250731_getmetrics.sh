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

directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_ala_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_ala_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'

directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16g2_1]='/home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16g2_2]='/home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16g2_3]='/home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb'

directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_m8g2_1]='/home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_m8g2_2]='/home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_m8g2_3]='/home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb'

for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
