#!/bin/bash

# Wed, 06 Aug 2025 10:21:46 -0700
# Use mdtraj to center molecules of interest within periodic boundaries

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/center_trajectory.sbatch'
T_NAME="output.pdb"

declare -a directories=(
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_ala_1"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_ala_2"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_ala_3"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16_1"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16_2"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16_3"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16g2_1"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16g2_2"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_c16g2_3"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_m8g2_1"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_m8g2_2"
    "/home/wesg/openmm_scripts/20250731_SMD/20250731T190851_m8g2_3"
    )

for dir in "${directories[@]}"; do
  cd $dir
  sbatch $SLURMSCRIPT $dir $T_NAME
done
