#!/bin/bash

# Wesley Gomersall
# Thu, 17 Jul 2025 15:16:31 -0700
#
# Launch slurm jobs on Talapas:
# Run steered molecular dynamics openmm script
# 5x replicates of ComD-C16 and ComD-alanines, add explicit solvent 
# I changed a couple things with this run: 
#   nonbondedMethod=PME again
#   Added in the constraints with Hbonds again
#   Only calculating the displacement and aCarbon distances when I am writing

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd_nohold.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts 
PULLFORCE=5

sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_3 $PULLFORCE

sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_3 $PULLFORCE
