#!/bin/bash

# These are runs in a nonperiodic system in which there is a force pulling the peptide in direction determined by centers of mass and a force defining
# a spherical boundary on the entire system.

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250819_SMD
PULLFORCE=5

mkdir -p $OUTPUTDIR

# C16-11.33
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11.33_fixed.pdb $OUTPUTDIR c16_11.33_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11.33_fixed.pdb $OUTPUTDIR c16_11.33_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11.33_fixed.pdb $OUTPUTDIR c16_11.33_3 $PULLFORCE

# C16-11.66
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11.66_fixed.pdb $OUTPUTDIR c16_11.66_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11.66_fixed.pdb $OUTPUTDIR c16_11.66_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11.66_fixed.pdb $OUTPUTDIR c16_11.66_3 $PULLFORCE

# C16-11 
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11_fixed.pdb $OUTPUTDIR c16_11_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11_fixed.pdb $OUTPUTDIR c16_11_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-11_fixed.pdb $OUTPUTDIR c16_11_3 $PULLFORCE

# C16-12
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR c16_12_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR c16_12_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR c16_12_3 $PULLFORCE

# C16-13
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb $OUTPUTDIR c16_13_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb $OUTPUTDIR c16_13_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb $OUTPUTDIR c16_13_3 $PULLFORCE

# C16-14
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb $OUTPUTDIR c16_14_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb $OUTPUTDIR c16_14_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb $OUTPUTDIR c16_14_3 $PULLFORCE

# C16-15
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb $OUTPUTDIR c16_15_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb $OUTPUTDIR c16_15_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb $OUTPUTDIR c16_15_3 $PULLFORCE

# PMPNN neg
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb $OUTPUTDIR pmpnn_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb $OUTPUTDIR pmpnn_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb $OUTPUTDIR pmpnn_3 $PULLFORCE

# Alanine chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_3 $PULLFORCE

# C16 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_3 $PULLFORCE

# C16G2 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb $OUTPUTDIR c16g2_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb $OUTPUTDIR c16g2_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb $OUTPUTDIR c16g2_3 $PULLFORCE

# M8G2 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR m8g2_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR m8g2_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR m8g2_3 $PULLFORCE

PULLFORCE=0
# Alanine chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR noforce_ala_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR noforce_ala_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR noforce_ala_3 $PULLFORCE

# C16 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR noforce_c16_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR noforce_c16_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR noforce_c16_3 $PULLFORCE
