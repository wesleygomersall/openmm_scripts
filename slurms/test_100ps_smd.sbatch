#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=smd
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
#SBATCH --partition=gpu
#SBATCH --ntasks=1                  # Number of tasks per node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --time=1-00:00:00           # Time limit 1 day
#SBATCH --mem=64G                   # Memory per node
#SBATCH --gpus=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu

# Wesley Gomersall
# Tue, 15 Jul 2025 14:13:35 -0700
# Run 0.1ns long simulation with 100 data entries, do not add explicit solvent, 
#   suppress movement of residues specified in `hold.txt` if true.

if [ $# != 5 ]; then
    echo "Supply 5 arguments in the following order:"
    echo "--input-pdb --output-dir --simulation-name --pull-force --hold(true or false)"
    exit
fi

INPUT_FILE=$1
OUT_PARENT_DIR=$2
NAME=$3
FORCE=$4
HOLD=$5

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

OUTPUT_DIR=$OUT_PARENT_DIR/$(date +"%Y%m%dT%H%M%S")_$NAME

if [ $HOLD == 'true' ]; then
    python3 ~/openmm_scripts/pullcomplex_md.py \
        --input $INPUT_FILE \
        --output $OUTPUT_DIR \
        --steps 50000 \
        --freq 100 \
        --pull-force $FORCE \
        --suppress-movement
    exit
fi

if [ $HOLD != 'true' ]; then
    python3 ~/openmm_scripts/pullcomplex_md.py \
        --input $INPUT_FILE \
        --output $OUTPUT_DIR \
        --steps 50000 \
        --freq 100 \
        --pull-force $FORCE 
    exit
fi
