#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=mdmetrics
#SBATCH --output=logs/mdtraj-%x_%j.out
#SBATCH --error=logs/mdtraj-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Wesley Gomersall
# Friday, 30 January, 2026, 14:17:39 

if [ $# != 4 ]; then
    echo "Supply 4 arguments in the following order:"
    echo "--trajectory-pdb --trajectory-dcd --reference-pdb --chain-info-file"
    exit
fi

TRAJ=$1
DCD=$2
REF=$3
CHINFO=$4

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

DISP_SCRIPT=~/openmm_scripts/mdtraj_metrics/displacement.py
RMSD_SCRIPT=~/openmm_scripts/mdtraj_metrics/rmsd.py
RMSD_DELTAS=~/openmm_scripts/mdtraj_metrics/delta_rmsd.py
RMSF_SCRIPT=~/openmm_scripts/mdtraj_metrics/rmsf.py
python3 $DISP_SCRIPT -i $TRAJ -d $DCD -r $REF -c $CHINFO
python3 $RMSD_SCRIPT -i $TRAJ -d $DCD -r $REF -c $CHINFO
python3 $RMSD_DELTAS -i $TRAJ -d $DCD -r $REF -c $CHINFO
python3 $RMSF_SCRIPT -i $TRAJ -d $DCD -r $REF -c $CHINFO

conda deactivate 
pigz *.pdb
exit
