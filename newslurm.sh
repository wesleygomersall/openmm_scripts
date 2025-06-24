#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=md
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
#SBATCH --partition=gpu
#SBATCH --ntasks=1                 # Number of tasks per node
#SBATCH --cpus-per-task=2         # Number of CPU cores per task
#SBATCH --time=24:00:00            # Time limit
#SBATCH --mem=64G                   # Memory per node
#SBATCH --gpus=1
#SBATCH --gres=gpu:1

# eval $(conda shell.bash hook) # for conda to work properly

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate openMM

module load cuda 

# python3 openmm_comdC16.py
python3 ~/openmm_scripts/pullcomplex_md.py \
    --input /projects/parisahlab/wesg/openmm/comd_c16_ala_af3multimer_pdbfixed.pdb \
    --output ~/openmm_scripts/testing1 \
    --steps 20000

# /projects/parisahlab/wesg/openmm/comd_c16_ala_af3multimer_pdbfixed.pdb
# other input for testing comd_c16_af3multimer_pdbfixed.pdb

# usage: pullcomplex_md.py [-h] [--input INPUT] [--output OUTPUT]
                         # [--timestep TIMESTEP] [--steps STEPS]
                         # [--pressure PRESSURE] [--temp TEMP]
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
  # --no-protein-move    Add this option to hold all atoms in the protein chain
                       # in place throughout the simulation.
