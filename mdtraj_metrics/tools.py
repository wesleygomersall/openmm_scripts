#!/usr/bin/env python3

import argparse
import mdtraj as md


parser = argparse.ArgumentParser(description="Calculate RMSD for second chain in MD trajectory with exactly two chains.")
parser.add_argument("--ref", type=str, help="Path to pdb used as input for MD simulation. This is used as reference for RMSD.")
parser.add_argument("--input", type=str, help="Path to trajectory pdb output from MD simulation.")
parser.add_argument("--output", type=str, help="File base name to output data in current directory.")
args = parser.parse_args()


class SmdTrajectory:
    '''
    Input file path of output trajectory from MD simulation. Trajectory should contain 
    2 chains, chain 0 is the protein, chain 1 is a peptide ligand.

    Attributes:
        traj: trajectory object with n frames of m atoms connected by topology.
        top: topology specifying the connection of atoms.
        protein: 
        peptide:
        peptide_atomids: list of integers containing atom IDs of the peptide's
                         atoms.
    '''
    def __init__(self, trajectoryfile: str):
        self.trajectory_solvent = md.load(trajectoryfile)
        self.trajectory = self.trajectory_solvent.remove_solvent()

        self.top = self.trajectory.topology

        self.protein = self.top.select("chainid == 0")
        self.peptide = self.top.select("chainid == 1")

        # Create a list of the IDs for all atoms of the peptide chain
        self.peptide_atomids = [atomid for atomid in self.peptide]


class Reference:
    '''
    Reference pdb structure for positions and topology. Properties of the
    reference must match trajectory of interest (create SmdTrajectory object
    first).

    Attributes:
        ref: reference pdb structure from md.load(file).

    Methods:
        check_ref_match: check if the number of chains, residues, and atoms
                         match between this reference structure and some
                         trajectory. Called upon __init__.
    '''
    def __init__(self, referencefile: str, smdtraj):
        self.ref = md.load(referencefile)
        self.check_ref_match(smdtraj)

    def check_ref_match(self, other):
        if other.trajectory.n_chains != 2:
            print(f"{other.trajectory.n_chains} chains in trajectory.")
            raise AssertionError("Number of chains in trajectory is not 2.")
        if other.trajectory.n_chains != self.ref.n_chains:
            print(f"{other.trajectory.n_chains} chains in trajectory.\n{self.ref.n_chains} chains in reference.")
            raise AssertionError("Number of chains in trajectory and reference differ.")
        if other.trajectory.n_residues != self.ref.n_residues:
            print(f"{other.trajectory.n_residues} residues in trajectory.\n{self.ref.n_residues} residues in reference.")
            raise AssertionError("Number of residues in trajectory and reference differ.")
        if other.trajectory.n_atoms != self.ref.n_atoms:
            print(f"{other.trajectory.n_atoms} atoms in trajectory.\n{self.ref.n_atoms} atoms in reference.")
            raise AssertionError("Number of atoms in trajectory and reference differ.")


if __name__ == "__main__":
    mytraj = SmdTrajectory(args.input)
    print(mytraj.trajectory)
    print(mytraj.top)
    print(mytraj.protein)
    print(mytraj.peptide)
    print(mytraj.peptide_atomids)
    myref = Reference(args.ref, mytraj)
