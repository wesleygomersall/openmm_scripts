#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=openmm
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
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

# c16-15
python3 ~/openmm_scripts/pullcomplex_md.py \
    --input  ~/Eckert2006/c16_15_fixed.pdb \
    --output ~/openmm_scripts/20250624_c16_15 \
    --steps 20000 \
    --freq 100 \
    --add-water \
    --pull-force 50

# c16-16
python3 ~/openmm_scripts/pullcomplex_md.py \
    --input  ~/Eckert2006/c16_16_fixed.pdb \
    --output ~/openmm_scripts/20250624_c16_16 \
    --steps 20000 \
    --freq 100 \
    --add-water \
    --pull-force 50

# c16-17
python3 ~/openmm_scripts/pullcomplex_md.py \
    --input  ~/Eckert2006/c16_17_fixed.pdb \
    --output ~/openmm_scripts/20250624_c16_17 \
    --steps 20000 \
    --freq 100 \
    --add-water \
    --pull-force 50

# c16-18
python3 ~/openmm_scripts/pullcomplex_md.py \
    --input  ~/Eckert2006/c16_18_fixed.pdb \
    --output ~/openmm_scripts/20250624_c16_18 \
    --steps 20000 \
    --freq 100 \
    --add-water \
    --pull-force 50

# usage: pullcomplex_md.py [-h] [--input INPUT] [--output OUTPUT]
                         # [--timestep TIMESTEP] [--steps STEPS]
                         # [--pressure PRESSURE] [--temp TEMP] [--add-water]
                         # [--no-protein-move]
# 
# MD trajectory simulating constant force pulling apart a protein-peptide
# complex.
# 
# options:
  # -h, --help           show this help message and exit
  # --input INPUT        Path to input pdb containing both protein and peptide.
                       # If errors exist, try using PDBfixer first.
  # --output OUTPUT      Path to output directory. Will be created if not
                       # already existing.
  # --timestep TIMESTEP  Time between steps in femtoseconds, default is 2.
  # --steps STEPS        Total steps in simulation, default is ten thousand.
  # --pressure PRESSURE  Pressure (atmospheres), default is 1 atm.
  # --temp TEMP          Temperature (Kelvin), default is 300K.
  # --add-water          Add this option to populate modeller with a surrounding
                       # box of water molecules.
  # --no-protein-move    Add this option to hold all atoms in the protein chain
                       # in place throughout the simulation.
