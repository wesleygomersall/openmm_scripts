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
# Monday, 01 December, 2025, 11:37:13 
# Run 0.5ns long simulation for troubleshoot

if [ $# != 4 ]; then
    echo "Supply 4 arguments in the following order:"
    echo "--input-pdb --output-dir --simulation-name --pull-force"
    exit
fi

INPUT_FILE=$1
OUT_PARENT_DIR=$2
NAME=$3
FORCE=$4

mkdir -p $OUT_PARENT_DIR

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

OUTPUT_DIR=$OUT_PARENT_DIR/$(date +"%Y%m%dT%H%M%S")_$NAME

python3 ~/openmm_scripts/pullcomplex_md.py \
    --input $INPUT_FILE \
    --output $OUTPUT_DIR \
    --steps 250000 \
    --freq 1000 \
    --add-water \
    --pull-force $FORCE \
    --suppress-movement 

# Don't do the rest for troubleshooting only openmm
exit 

# Get metrics from output.pdb with mdtraj
conda activate mdtraj

cd $OUTPUT_DIR

TRAJECTORY="*.pdb"
MDTRAJ_METRICS_OUTPUT='mdtrajanalysis.csv'

python3 /home/wesg/openmm_scripts/mdtraj_metrics/allstats.py \
    -i $TRAJECTORY \
    -r $INPUT_FILE \
    -o $MDTRAJ_METRICS_OUTPUT
