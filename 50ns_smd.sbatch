#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=smd
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
#SBATCH --partition=gpulong
#SBATCH --ntasks=1                  # Number of tasks per node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --time=7-00:00:00           # Time limit 7 days
#SBATCH --mem=64G                   # Memory per node
#SBATCH --gpus=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu


# TODO: 
# take input file, output dir, outputname as command line arguments

INPUTFILE=$1
OUTPARENT=$2
NAME=$3

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

OUTPUTDIR=$OUTPARENT/$(date +"%Y%m%dT%H%M%S")_$NAME

# this new script has a nonperiodic cutoff which will hopefully speed it up significantly 

python3 ~/openmm_scripts/pullcomplex_md_nonperiodic.py \
    --input $INPUTFILE \
    --output $OUTPUTDIR \
    --steps 25000000 \
    --freq 1000 \
    --add-water \
    --pull-force 5 \
    --suppress-movement
