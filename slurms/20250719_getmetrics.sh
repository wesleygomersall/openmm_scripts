#!/bin/bash

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/getmetrics.sbatch'
T_NAME='output.pdb'
OUTPUTNAME='mdtrajanalysis.csv'

declare -A directories_references
directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134502_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134527_c16_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134527_c16_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134526_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134538_ala_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250719_SMD/20250719T134557_ala_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
