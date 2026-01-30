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

if [ $# != 3 ]; then
    echo "Supply 3 arguments in the following order:"
    echo "--trajectory-pdb --reference-pdb --chain-info-file"
    exit
fi

TRAJ=$1
REF=$2
CHINFO=$3

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

DISP_SCRIPT=~/openmm_scripts/mdtraj_metrics/displacement.py
RMSD_SCRIPT=~/openmm_scripts/mdtraj_metrics/rmsd.py
RMSD_DELTAS=~/openmm_scripts/mdtraj_metrics/delta_rmsd.py
REM_SOLVENT=~/openmm_scripts/mdtraj_metrics/remove_solvent.py
RMSF_SCRIPT=~/openmm_scripts/mdtraj_metrics/rmsf.py
python3 $DISP_SCRIPT -i $TRAJ -r $REF -c $CHINFO
python3 $RMSD_SCRIPT -i $TRAJ -r $REF -c $CHINFO
python3 $RMSD_DELTAS -i $TRAJ -r $REF -c $CHINFO
python3 $RMSF_SCRIPT -i $TRAJ -r $REF -c $CHINFO
python3 $REM_SOLVENT -i $TRAJ -r $REF -c $CHINFO

conda deactivate 
pigz *.pdb
exit
