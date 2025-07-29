#!/bin/bash

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/getmetrics.sbatch'
T_NAME='output.pdb'
OUTPUTNAME='mdtrajanalysis.csv'

declare -A directories_references
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T205103_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T205103_ala_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T205103_ala_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T215351_pmpnn_1]='/home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T215351_pmpnn_2]='/home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T215351_pmpnn_3]='/home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T205103_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T205103_c16_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T205103_c16_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T215351_c16m12_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T215351_c16m12_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250707T215351_c16m12_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m13_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m13_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m13_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m14_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m14_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m14_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb'

directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m15_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m15_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb'
directories_references[/projects/parisahlab/wesg/20250707_SMD/20250708T133217_c16m15_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb'

for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
