#!/bin/bash

# Wed, 06 Aug 2025 10:21:46 -0700
# Use mdtraj to center molecules of interest within periodic boundaries

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/center_trajectory.sbatch'
T_NAME="output.pdb"

declare -a directories=(
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T121807_testsink_c16_1"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T125312_testsink40_c16"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T125306_testsink30_c16"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T121807_testsink_c16_2"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T124605_testsink50_ala_1"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T124605_testsink50_ala_2"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T121807_testsink_ala_2"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T125312_testsink20_c16"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T125306_testsink30_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T125312_testsink40_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T125306_testsink20_ala"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T121807_testsink_ala_1"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T124605_testsink50_c16_2"
    "/home/wesg/openmm_scripts/20250807_sinktest_SMD/20250807T124605_testsink50_c16_1"
    )

for dir in "${directories[@]}"; do
  cd $dir
  sbatch $SLURMSCRIPT $dir $T_NAME
done
