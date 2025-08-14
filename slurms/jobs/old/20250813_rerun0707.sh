#!/bin/bash

SCRIPT=/home/wesg/openmm_scripts/slurms/0.5ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250813_0707_test_SMD
PULLFORCE=20

mkdir -p $OUTPUTDIR

# Alanine chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_directf_container1 $PULLFORCE

# C16 chain
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_directf_container1 $PULLFORCE
