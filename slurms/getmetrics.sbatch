#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=mdmetrics
#SBATCH --output=logs/mdtraj-%x_%j.out
#SBATCH --error=logs/mdtraj-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Wesley Gomersall
# Tue, 22 Jul 2025 18:11:45 -0700
# 
# Get metrics for SMD run, create csv with the following information: 
# Step,Displacement(nm),RMSD(nm),aCarbon-aCarbon Distance(nm)

if [ $# != 4 ]; then
    echo "Supply 4 arguments in the following order:"
    echo "--path --trjectory-pdb --reference-pdb --output-name"
    exit
fi

PATH=$1
TRAJECTORY=$2
REFERENCE=$3
OUTPUT=$4

SCRIPT='~/openmm_scripts/mdtraj_metrics/allstats.py'

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

cd $PATH

$SCRIPT -i $TRAJECTORY -r $REFERENCE -o $OUTPUT
