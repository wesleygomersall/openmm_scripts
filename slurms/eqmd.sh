#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=mdeq
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
# Wed, 10 Sep 2025 10:32:39 -0700
# Run long simulation with 1000 data entries, add explicit solvent

if [ $# != 5 ]; then
    echo "Supply 5 arguments in the following order:"
    echo "--input-pdb --output-dir --simulation-name --timestep --steps"
    exit
fi

OPENMM_INPUT=$1
OUT_PARENT_DIR=$2
NAME=$3
TSTEP=$4
STEPS=$5

# source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
# conda activate openMM
conda init 
conda activate openmm

# module load cuda 

OUTPUT_DIR=$OUT_PARENT_DIR/$(date +"%Y%m%dT%H%M%S")_$NAME

mkdir -p $OUT_PARENT_DIR

# usage: MD_longeq.py [-h] [--input INPUT] [--output OUTPUT]
                    # [--timestep TIMESTEP] [--steps STEPS] [--freq FREQ]
                    # [--pressure PRESSURE] [--temp TEMP] [--add-water]
                    # [--variable-int]
# 
# MD trajectory simulating a protein-peptide complex over time.
# 
# options:
  # -h, --help           show this help message and exit
  # --input INPUT        Path to input pdb containing both protein and peptide.
                       # If errors exist, try using PDBfixer first.
  # --output OUTPUT      Path to output directory. Will be created if not
                       # already existing.
  # --timestep TIMESTEP  Time between steps in femtoseconds, default is 2.
  # --steps STEPS        Total steps in simulation, default is ten thousand.
  # --freq FREQ          Frequency of reporter, will record data `freq` amount
                       # of times (every `step // freq` steps).
  # --pressure PRESSURE  Pressure (atmospheres), default is 1 atm.
  # --temp TEMP          Temperature (Kelvin), default is 300K.
  # --add-water          Add this option to populate modeller with a surrounding
                       # box of water molecules.
  # --variable-int       Add this option to use variable integrator for
                       # constructing the system.


python3 ~/openmm_scripts/MD_longeq.py \
    --input $OPENMM_INPUT \
    --output $OUTPUT_DIR \
    --timestep $TSTEP \
    --steps $STEPS \
    --freq 1000 \
    --add-water \
    --variable-int
