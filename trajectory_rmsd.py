#!/usr/bin/env python3

import mdtraj as md
import argparse

parser = argparse.ArgumentParser(description="Calculate RMSD for second chain in MD trajectory with exactly two chains.")
parser.add_argument("--ref", type=str, help="Path to pdb used as input for MD simulation. This is used as reference for RMSD.")
parser.add_argument("--input", type=str, help="Path to trajectory pdb output from MD simulation.")
parser.add_argument("--output", type=str, help="File base name to output data in current directory.")
args = parser.parse_args()

traj = md.load(args.input)
ref = md.load(args.ref)

if traj.n_chains != 2:
    raise AssertionError("Number of chains in trajectory is not 2.")
if traj.n_chains != ref.n_chains:
    print(f"{traj.n_chains} chains in trajectory.\n{ref.n_chains} chains in reference.")
    raise AssertionError("Number of chains in trajectory and reference differ.")
if traj.n_residues != ref.n_residues:
    print(f"{traj.n_residues} residues in trajectory.\n{ref.n_residues} residues in reference.")
    raise AssertionError("Number of residues in trajectory and reference differ.")
if traj.n_atoms != ref.n_atoms:
    print(f"{traj.n_atoms} atoms in trajectory.\n{ref.n_atoms} atoms in reference.")
    raise AssertionError("Number of atoms in trajectory and reference differ.")

top = traj.topology
protein = top.select("chainid == 0")
peptide = top.select("chainid == 1")

# Create a list of the IDs for all atoms of the peptide chain
peptide_atomids = [atomid for atomid in peptide]

# Calculate RMSD for each frame, using the 0th (from ref) as reference
peptide_rmsd = md.rmsd(traj, ref, 0, peptide_atomids)

# Create csv with columns: Frame, RMSD
with open(args.output.split('.csv')[0] + ".csv", 'w') as out: 
    out.write("Frame,RMSD(nm)\n")
    for i, frame_rmsd in enumerate(peptide_rmsd): 
        out.write(f"{i+1},{frame_rmsd}\n")
