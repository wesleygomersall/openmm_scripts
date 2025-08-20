#!/bin/bash

SCRIPT=/home/wesg/openmm_scripts/slurms/0.5ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250818_tests_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=20
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_testshort_20 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_testshort_20 $PULLFORCE

PULLFORCE=10
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_testshort_10 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_testshort_10 $PULLFORCE

PULLFORCE=5
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR ala_testshort_5 $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR c16_testshort_5 $PULLFORCE
