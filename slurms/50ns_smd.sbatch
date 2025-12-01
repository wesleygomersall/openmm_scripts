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
# Tue, 02 Sep 2025 18:57:11 -0700
# Run 50ns long simulation with 1000 data entries, add explicit solvent
# After the trajectory is created, superpose it and calculate metrics with mdtraj 

if [ $# != 4 ]; then
    echo "Supply 4 arguments in the following order:"
    echo "--input-pdb --output-dir --simulation-name --pull-force"
    exit
fi

OPENMM_INPUT=$1
OUT_PARENT_DIR=$2
NAME=$3
FORCE=$4

mkdir -p $OUT_PARENT_DIR

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

OUTPUT_DIR=$OUT_PARENT_DIR/$(date +"%Y%m%dT%H%M%S")_$NAME

python3 ~/openmm_scripts/pullcomplex_md.py \
    --input $OPENMM_INPUT \
    --output $OUTPUT_DIR \
    --steps 25000000 \
    --freq 1000 \
    --add-water \
    --pull-force $FORCE

# For now, move this up here to avoid wasting time on superpose.py and allstats.py
pigz *.pdb
exit

# openmm output
TRAJECTORY="$OUTPUT_DIR/output.pdb"

# Get metrics from output.pdb with mdtraj
conda activate mdtraj

cd $OUTPUT_DIR

python3 /home/wesg/openmm_scripts/mdtraj_metrics/superpose.py \
    -i $TRAJECTORY \
    -r $OPENMM_INPUT

SUPERPOSED="$OUTPUT_DIR/output_superposed.pdb"
MDTRAJ_METRICS_OUTPUT='mdtrajanalysis.csv'

python3 /home/wesg/openmm_scripts/mdtraj_metrics/allstats.py \
    -i $SUPERPOSED \
    -r $OPENMM_INPUT \
    -o $MDTRAJ_METRICS_OUTPUT

