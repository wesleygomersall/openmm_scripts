#!/bin/bash

#SBATCH --job-name=md
#SBATCH --output=logs/openmm-%A_%a.out
#SBATCH --error=logs/openmm-%A_%a.err
#SBATCH --partition=gpu
#SBATCH --ntasks=1                 # Number of tasks per node
#SBATCH --cpus-per-task=2         # Number of CPU cores per task
#SBATCH --time=24:00:00            # Time limit
#SBATCH --mem=64G                   # Memory per node
#SBATCH --gpus=1
#SBATCH --gres=gpu:1

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

python3 openmm_comdC16_ala.py
