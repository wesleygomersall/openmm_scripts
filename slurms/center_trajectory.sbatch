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
# Wed, 06 Aug 2025 09:55:17 -0700
# 
# Use mdtraj to center chains 0 and 1 within periodic boundary conditions. This
# should minimize jumping across boundary while peptide (chain 1) is bound. 

if [ $# != 2 ]; then
    echo "Supply 2 arguments in the following order:"
    echo "--path --trjectory-pdb"
    exit
fi

MYDIR=$1
TRAJECTORY=$2
REFERENCE="notneeded"

SCRIPT='/home/wesg/openmm_scripts/mdtraj_metrics/clean_trajectory_pbc.py'

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

cd $MYDIR

$SCRIPT -i $TRAJECTORY -r $REFERENCE
