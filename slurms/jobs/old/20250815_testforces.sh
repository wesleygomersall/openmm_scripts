#!/bin/bash

SCRIPT=/home/wesg/openmm_scripts/slurms/0.5ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250815_tests_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=20
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_test1_20 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_test1_20 $PULLFORCE

PULLFORCE=10
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_test1_10 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_test1_10 $PULLFORCE

PULLFORCE=5
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_test1_5 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_test1_5 $PULLFORCE
