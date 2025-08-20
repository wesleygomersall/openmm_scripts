#!/bin/bash

# Template to initiate slurm jobs which run MD simulations with openmm
# This template runs MD with constrained protein by setting mass of some atoms to zero.  

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd_hold.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts 
PULLFORCE=5

mkdir -p $OUTPUTDIR

# Run triplicate
# sbatch $SCRIPT /path/to/input.pdb $OUTPUTDIR creative_run_name_1 $PULLFORCE
# sbatch $SCRIPT /path/to/input.pdb $OUTPUTDIR creative_run_name_2 $PULLFORCE
# sbatch $SCRIPT /path/to/input.pdb $OUTPUTDIR creative_run_name_3 $PULLFORCE

# Run triplicate
# sbatch $SCRIPT /path/to/input.pdb $OUTPUTDIR creative_run_name_1 $PULLFORCE
# sbatch $SCRIPT /path/to/input.pdb $OUTPUTDIR creative_run_name_2 $PULLFORCE
# sbatch $SCRIPT /path/to/input.pdb $OUTPUTDIR creative_run_name_3 $PULLFORCE
