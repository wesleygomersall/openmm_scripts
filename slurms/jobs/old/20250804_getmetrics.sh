#!/bin/bash

# Tue, 05 Aug 2025 14:29:14 -0700
# Use mdtraj to get metrics from MD trajectory pdb after the edits made to
# account for periodic boundary conditions. 

SLURMSCRIPT='/home/wesg/openmm_scripts/slurms/getmetrics.sbatch'
T_NAME='output.pdb'
OUTPUTNAME='mdtrajanalysis.csv'

# Associative array for trajectories and their reference files
declare -A directories_references

directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.33_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11.33_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.33_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11.33_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.33_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11.33_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.66_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11.66_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.66_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11.66_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_11.66_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-11.66_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_12_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_12_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_12_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-12_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124717_c16_13_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T124812_c16_13_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T132756_c16_13_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-13_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T134027_c16_14_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_14_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_14_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-14_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_15_1]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_15_2]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_c16_15_3]='/home/wesg/proteinpdbs/Eckert2006/eckertC16-15_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_ala_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_ala_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_ala_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_allalanine_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_c16_1]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_c16_2]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_noforce_c16_3]='/home/wesg/proteinpdbs/ComD_C16_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_pmpnn_1]='/home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_pmpnn_2]='/home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202146_pmpnn_3]='/home/wesg/proteinpdbs/Eckert2006/mpnn_neg_ctrl1_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T202512_noforce_c16g2_1]='/home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T203650_noforce_c16g2_2]='/home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_c16g2_3]='/home/wesg/proteinpdbs/comd_c16g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_m8g2_1]='/home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_m8g2_2]='/home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb'
directories_references[/home/wesg/openmm_scripts/20250804_SMD/20250804T224845_noforce_m8g2_3]='/home/wesg/proteinpdbs/comd_m8g2_AF3_fixed.pdb'

for key in "${!directories_references[@]}"; do 
  sbatch $SLURMSCRIPT $key $T_NAME ${directories_references[$key]} $OUTPUTNAME
done
