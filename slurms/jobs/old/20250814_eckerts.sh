#!/bin/bash

# These are runs in a nonperiodic system in which there is a force pulling the peptide in direction determined by centers of mass and a force defining
# a spherical boundary on the entire system.

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250814_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=5
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
