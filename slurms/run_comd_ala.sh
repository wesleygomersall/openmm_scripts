#!/bin/bash

# Wesley Gomersall
# Mon, 07 Jul 2025 19:21:29 -0700
#
# Script to launch slurm jobs on Talapas:
# Run steered molecular dynamics openmm script

SCRIPT=/home/wesg/openmm_scripts/slurms/50ns_smd.sbatch

sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg

sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg



