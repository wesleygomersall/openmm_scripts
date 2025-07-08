#!/bin/bash

# Wesley Gomersall
# Mon, 07 Jul 2025 19:21:29 -0700
# Mon, 07 Jul 2025 20:32:35 -0700
#
# Launch slurm jobs on Talapas:
# Run steered molecular dynamics openmm script
# Triplicate of ComD-C16 and triplicate of ComD-alanines

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts 
PULLFORCE=5

sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_3 $PULLFORCE

sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_3 $PULLFORCE
