#!/bin/bash

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20260115_comd_2shortmuts_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=5

INPUT=/home/wesg/openmm_scripts/data/20260112_comd_shortmutants_MDinput/comdn230-c16-12_seed3_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20260112_comd_shortmutants_MDinput/ChainInfo/comdn230_12.chinfo

sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep3 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep4 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep5 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep6 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep7 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep8 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8_rep9 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20260112_comd_shortmutants_MDinput/comdn230-c16-m8_14_seed3_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20260112_comd_shortmutants_MDinput/ChainInfo/comdn230_m8_14.chinfo

sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep3 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep4 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep5 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep6 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep7 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep8 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR comdn230_M8-14_rep9 $PULLFORCE $CHINFO
