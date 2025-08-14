#!/bin/bash

# Wed, 06 Aug 2025 10:21:46 -0700
# Use mdtraj to center molecules of interest within periodic boundaries

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/center_trajectory.sbatch'
T_NAME="output.pdb"

declare -a directories=(
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151700_testsink30_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151700_testsink40_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151700_testsink40_c16"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151900_testsink30_c16"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151926_testsink20_c16"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151929_testsink10_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151929_testsink20_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest2_SMD/20250807T151957_testsink10_c16"
    )

for dir in "${directories[@]}"; do
  cd $dir
  sbatch $SLURMSCRIPT $dir $T_NAME
done
