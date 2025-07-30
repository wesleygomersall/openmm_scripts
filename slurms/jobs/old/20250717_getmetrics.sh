#!/bin/bash

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/getmetrics.sbatch'
T_NAME='output.pdb'
OUTPUTNAME='mdtrajanalysis.csv'

declare -A directories_references
directories_references[/projects/parisahlab/wesg/20250717_SMD/20250717T151842_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250717_SMD/20250717T151853_ala_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250717_SMD/20250717T151902_ala_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250717_SMD/20250717T151828_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250717_SMD/20250717T151834_c16_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250717_SMD/20250717T151834_c16_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'


for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
