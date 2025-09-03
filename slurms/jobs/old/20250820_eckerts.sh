#!/bin/bash

# These are runs in a nonperiodic system in which there is a force pulling the
# peptide in direction determined by centers of mass and a force defining a
# spherical boundary on the entire system. Two residues are held in the ComD
# protein specified in file 

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd_hold.sbatch
OUTPUTDIR=/home/wesg/openmm_scripts/20250820_SMD
PULLFORCE=5
HOLDINGFILE=/home/wesg/openmm_scripts/comd_hold.txt

mkdir -p $OUTPUTDIR

# C16-16
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-16_fixed.pdb $OUTPUTDIR c16_16_1 $PULLFORCE $HOLDINGFILE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-16_fixed.pdb $OUTPUTDIR c16_16_2 $PULLFORCE $HOLDINGFILE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-16_fixed.pdb $OUTPUTDIR c16_16_3 $PULLFORCE $HOLDINGFILE

# C16-17
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-17_fixed.pdb $OUTPUTDIR c16_17_1 $PULLFORCE $HOLDINGFILE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-17_fixed.pdb $OUTPUTDIR c16_17_2 $PULLFORCE $HOLDINGFILE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-17_fixed.pdb $OUTPUTDIR c16_17_3 $PULLFORCE $HOLDINGFILE

# C16-18
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-18_fixed.pdb $OUTPUTDIR c16_18_1 $PULLFORCE $HOLDINGFILE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-18_fixed.pdb $OUTPUTDIR c16_18_2 $PULLFORCE $HOLDINGFILE
sbatch $SCRIPT /home/wesg/proteinpdbs/Eckert2006/eckertC16-18_fixed.pdb $OUTPUTDIR c16_18_3 $PULLFORCE $HOLDINGFILE
