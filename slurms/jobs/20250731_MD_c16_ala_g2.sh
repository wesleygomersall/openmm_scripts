#!/bin/bash

# Template to initiate slurm jobs which run MD simulations with openmm

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250731_SMD
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

# M8G2 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR m8g2_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR m8g2_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb $OUTPUTDIR m8g2_3 $PULLFORCE
