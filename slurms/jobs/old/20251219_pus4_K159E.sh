#!/bin/bash

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20251219_pus4mutants_pull5_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=5

INPUT=/home/wesg/openmm_scripts/data/20251219_pus4_MDinput/pus4_tRNA_K159E_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251219_pus4_MDinput/ChainInfo/K159E_chaininfo

sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep3 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep4 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep5 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep6 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep7 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep8 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_K159E_rep9 $PULLFORCE $CHINFO

