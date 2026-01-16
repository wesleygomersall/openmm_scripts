#!/bin/bash

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20251210_pus4mutants_pull5_SMD

mkdir -p $OUTPUTDIR

PULLFORCE=5

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_1_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/1_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_1rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_1rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_1rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_2_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/2_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_2rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_2rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_2rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_3_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/3_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_3rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_3rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_3rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_4_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/4_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_4rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_4rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_4rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_5_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/5_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_5rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_5rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_5rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_6_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/6_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_6rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_6rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_6rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_7_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/7_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_7rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_7rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_7rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_8_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/8_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_8rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_8rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_8rep3 $PULLFORCE $CHINFO

INPUT=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/pus4_tRNA_9_fixed.pdb
CHINFO=/home/wesg/openmm_scripts/data/20251203_pus4_MDinput/ChainInfo/9_chaininfo
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_9rep1 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_9rep2 $PULLFORCE $CHINFO
sbatch $SCRIPT $INPUT $OUTPUTDIR pus4_9rep3 $PULLFORCE $CHINFO

exit 
