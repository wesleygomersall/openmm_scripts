#!/bin/bash

# Wed, 06 Aug 2025 10:04:02 -0700
# Use mdtraj to center molecules of interest within periodic boundaries

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/center_trajectory.sbatch'
T_NAME="output.pdb"

declare -a directories=(
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.33_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.33_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.33_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.66_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.66_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.66_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_12_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_12_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_12_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_13_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T124812_c16_13_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T132756_c16_13_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T134027_c16_14_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_14_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_14_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_15_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_15_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_15_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_ala_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_ala_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_ala_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_c16_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_c16_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_c16_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_pmpnn_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_pmpnn_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_pmpnn_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T202512_noforce_c16g2_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T203650_noforce_c16g2_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_c16g2_3"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_m8g2_1"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_m8g2_2"
    "/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_m8g2_3"
    )

for dir in "${directories[@]}"; do
  cd $dir
  sbatch $SLURMSCRIPT $dir $T_NAME
done
