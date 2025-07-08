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

# Wesley Gomersall
# Mon, 07 Jul 2025 19:52:05 -0700
# Run 50ns long simulation with 1000 data entries, add explicit solvent, 
#   suppress movement of residues specified in `hold.txt`.

if [ $# != 4 ]; then
    echo "Supply 4 arguments in the following order:"
    echo "--input-pdb --output-dir --simulation-name --pull-force"
    exit
fi

INPUTFILE=$1
OUTPARENT=$2
NAME=$3
FORCE=$4

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

OUTPUTDIR=$OUT_PARENT_DIR/$(date +"%Y%m%dT%H%M%S")_$NAME

python3 ~/openmm_scripts/pullcomplex_md.py \
    --input $INPUTFILE \
    --output $OUTPUTDIR \
    --steps 25000000 \
    --freq 1000 \
    --add-water \
    --pull-force $4 \
    --suppress-movement
