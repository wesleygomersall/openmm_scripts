#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=smd
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
#SBATCH --partition=gpulong
#SBATCH --ntasks=1                  # Number of tasks per node
#SBATCH --cpus-per-task=4           # Number of CPU cores per task
#SBATCH --time=2-00:00:00           # Time limit 2 days
#SBATCH --mem=64G                   # Memory per node
#SBATCH --gpus=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

python3 ~/openmm_scripts/pullcomplex_md.py \
    --input /projects/parisahlab/wesg/openmm/comd_c16_af3multimer_pdbfixed.pdb \
    --output ~/openmm_scripts/20250626_c16_long \
    --steps 50000000 \
    --freq 1000 \
    --add-water \
    --pull-force 20
