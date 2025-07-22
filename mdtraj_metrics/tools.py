#!/usr/bin/env python3

import sys
import argparse
import mdtraj as md


parser = argparse.ArgumentParser(description='Calculate RMSD for second chain in MD trajectory with exactly two chains.')
parser.add_argument('-i', '--input', type=str, help='Path to trajectory pdb output from MD simulation.')
parser.add_argument('-r', '--ref', type=str, help='Path to pdb used as input for MD simulation. This is used as reference for RMSD.')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help='Output file path (default: stdout).')
args = parser.parse_args()


def check_ref_match(trajectory, ref):
    if trajectory.n_chains != 2:
        print(f"{trajectory.n_chains} chains in trajectory.")
        raise AssertionError("Number of chains in trajectory is not 2.")
    if trajectory.n_chains != ref.n_chains:
        print(f"{trajectory.n_chains} chains in trajectory.\n{ref.n_chains} chains in reference.")
        raise AssertionError("Number of chains in trajectory and reference differ.")
    if trajectory.n_residues != ref.n_residues:
        print(f"{trajectory.n_residues} residues in trajectory.\n{ref.n_residues} residues in reference.")
        raise AssertionError("Number of residues in trajectory and reference differ.")
    if trajectory.n_atoms != ref.n_atoms:
        print(f"{trajectory.n_atoms} atoms in trajectory.\n{ref.n_atoms} atoms in reference.")
        raise AssertionError("Number of atoms in trajectory and reference differ.")


def prepend_reference(trajectory, reference): 
    check_ref_match(trajectory, reference)
    newtraj1 = trajectory
    newtraj2 = reference
    newtraj1.unitcell_angles = None
    newtraj1.unitcell_lengths = None
    newtraj2.unitcell_angles = None
    newtraj2.unitcell_lengths = None
    return newtraj2 + newtraj1


if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    myref = md.load(args.ref)
    newtraj = prepend_reference(mytraj, myref)
