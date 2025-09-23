#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=metrics
#SBATCH --output=logs/openmm-%x_%j.out
#SBATCH --error=logs/openmm-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G

# Wesley Gomersall
# Tue, 16 Sep 2025 12:09:58 -0700
# Collect metrics from openMM trajectory pdb. 
# Superpose and metrics from mdtraj 

if [ $# != 3 ]; then
    echo "Supply 3 arguments in the following order:"
    echo "--parent-dir --trajectory-pdb --reference-pdb" # --simulation-name --pull-force"
    exit
fi

PARENT_DIR=$1
RAW_TRAJ=$2
REFERENCE=$3

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

cd $PARENT_DIR

python3 /home/wesg/openmm_scripts/mdtraj_metrics/superpose.py \
    -i $RAW_TRAJ \
    -r $REFERENCE

SUPERPOSED="$OUTPUT_DIR/output_superposed.pdb"
MDTRAJ_METRICS_OUTPUT='mdtrajanalysis.csv'

python3 /home/wesg/openmm_scripts/mdtraj_metrics/allstats.py \
    -i $SUPERPOSED \
    -r $REFERENCE \
    -o $MDTRAJ_METRICS_OUTPUT

pigz *.pdb
