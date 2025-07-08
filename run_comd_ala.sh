#!/bin/bash

# Wesley Gomersall
# Mon, 07 Jul 2025 19:21:29 -0700
#
# Script to launch slurm jobs on Talapas:
# Run steered molecular dynamics openmm script

SCRIPT=/path/to/the/script

sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg

sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg
sbatch $SCRIPT firstarg secondarg thirdarg



