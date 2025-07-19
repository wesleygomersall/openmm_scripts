#!/bin/bash

# Wesley Gomersall
# Mon, 07 Jul 2025 20:32:35 -0700
#
# Launch slurm jobs on Talapas:
# I wanted to replicate the run from 20250707.
# I retrieved the version of the python script from GitHub.
# That version was saved as pullcomplex_md_V20250707.py or something. 
# I must create a copy of the slurm script and use that as $SCRIPT here. 

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd_nohold2.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts 
PULLFORCE=5

sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_3 $PULLFORCE

sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_1 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_2 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_3 $PULLFORCE
