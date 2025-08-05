#!/bin/bash

# This is a repeat of the script ran on 2025-07-31, without the pull force. 
# Changes: PULLFORCE=0 and added 'noforce_' to the output file. 

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250804_SMD
PULLFORCE=0

mkdir -p $OUTPUTDIR

# Alanine chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR noforce_ala_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR noforce_ala_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR noforce_ala_3 $PULLFORCE

# C16 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR noforce_c16_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR noforce_c16_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR noforce_c16_3 $PULLFORCE

# C16G2 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb $OUTPUTDIR noforce_c16g2_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb $OUTPUTDIR noforce_c16g2_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb $OUTPUTDIR noforce_c16g2_3 $PULLFORCE

# M8G2 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR noforce_m8g2_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR noforce_m8g2_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR noforce_m8g2_3 $PULLFORCE
