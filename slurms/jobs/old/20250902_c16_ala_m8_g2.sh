#!/bin/bash

# These are runs in a nonperiodic system in which there is a force pulling the peptide in direction determined by centers of mass and a force defining
# a spherical boundary on the entire system.

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250902_SMD
PULLFORCE=5

mkdir -p $OUTPUTDIR

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

# M8 chain 
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR m8_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR m8_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR m8_3 $PULLFORCE

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
