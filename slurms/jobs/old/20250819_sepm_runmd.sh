#!/bin/bash

# Run MD with and without holding for target SepM in complex with CSP21. There
# is no force added to the peptide in order to attempt to replicate RMSF
# calculations from Liu 2024 paper. 

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd_hold.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250819_sepm_SMD
PULLFORCE=0

HOLDFILE=/home/wesg/openmm_scripts/sepM_hold.txt

mkdir -p $OUTPUTDIR

INPUTPDB=/home/wesg/proteinpdbs/sepmcsp21/sepmcsp21_model_fixed.pdb

sbatch $SCRIPT $INPUTPDB $OUTPUTDIR sepm_hold_1 $PULLFORCE $HOLDFILE
sbatch $SCRIPT $INPUTPDB $OUTPUTDIR sepm_hold_2 $PULLFORCE $HOLDFILE
sbatch $SCRIPT $INPUTPDB $OUTPUTDIR sepm_hold_3 $PULLFORCE $HOLDFILE

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
sbatch $SCRIPT $INPUTPDB $OUTPUTDIR sepm_nohold_1 $PULLFORCE
sbatch $SCRIPT $INPUTPDB $OUTPUTDIR sepm_nohold_2 $PULLFORCE
sbatch $SCRIPT $INPUTPDB $OUTPUTDIR sepm_nohold_3 $PULLFORCE



