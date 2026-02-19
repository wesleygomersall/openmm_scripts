#!/bin/bash

#SBATCH --account=parisahlab
#SBATCH --job-name=traj_contacts
#SBATCH --output=logs/job_desc-%x_%j.out
#SBATCH --error=logs/job_desc-%x_%j.err
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=wesg@uoregon.edu

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate mdtraj

if [ "$#" -eq 3 ]; then
    python3 ~/openmm_scripts/mdtraj_metrics/get_contacts.py -i $1 -r $2 -c $3
elif [ "$#" -eq 4 ]; then
    python3 ~/openmm_scripts/mdtraj_metrics/get_contacts.py -i $1 -d $2 -r $3 -c $4
else
    echo "Supply 3 or 4 arguments in the following order:"
    echo "--trajectory-pdb (OPTIONAL --output-dcd) --reference-pdb --chain-info-file"
    exit
fi

conda deactivate 
exit
