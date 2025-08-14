#!/bin/bash

# Test the new sink attractor within periodic boundary

SCRIPT=/home/wesg/openmm_scripts/slurms/0.5ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250807_sinktest2_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=40
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_ala $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_c16 $PULLFORCE

PULLFORCE=30
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_ala $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_c16 $PULLFORCE

PULLFORCE=20
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_ala $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_c16 $PULLFORCE

PULLFORCE=10
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_ala $PULLFORCE
sbatch $SCRIPT /home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb $OUTPUTDIR testsink${PULLFORCE}_c16 $PULLFORCE
