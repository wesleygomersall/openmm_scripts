#!/bin/bash

# Wesley Gomersall
# Mon, 07 Jul 2025 21:28:11 -0700
#
# Launch slurm jobs on Talapas:
# Run steered molecular dynamics openmm script
# Triplicate of one PMPNN negative control and triplicate of some ComD-C16
#   mutants from and inspired by Eckert et al. 2006. 

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts 
PULLFORCE=5

sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR c16m12_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR c16m12_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb $OUTPUTDIR c16m12_3 $PULLFORCE

sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb $OUTPUTDIR pmpnn_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb $OUTPUTDIR pmpnn_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb $OUTPUTDIR pmpnn_3 $PULLFORCE

# Put these in the next script:
# /home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb
# /home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb
# /home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb
