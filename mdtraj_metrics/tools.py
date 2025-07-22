#!/usr/bin/env python3

import argparse
import mdtraj as md


parser = argparse.ArgumentParser(description="Calculate RMSD for second chain in MD trajectory with exactly two chains.")
parser.add_argument("--ref", type=str, help="Path to pdb used as input for MD simulation. This is used as reference for RMSD.")
parser.add_argument("--input", type=str, help="Path to trajectory pdb output from MD simulation.")
parser.add_argument("--output", type=str, help="File base name to output data in current directory.")
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


class SmdTrajectory:
    '''
    Input file path of output trajectory from MD simulation. Trajectory should contain 
    2 chains, chain 0 is the protein, chain 1 is a peptide ligand.

    Attributes:
        traj: trajectory object with n frames of m atoms connected by topology.
        ref: reference pdb structure input to MD simulation.
        top: topology specifying the connection of atoms.
        protein: 
        peptide:
        peptide_atomids: list of integers containing atom IDs of the peptide's
                         atoms.

    Methods:
        check_ref_match: check if the number of chains, residues, and atoms
                         match between this reference structure and trajectory.
                         Called upon __init__.
    '''

    def __init__(self, trajectoryfile: str, referencefile: str):

        # store trajectories with and without solvent
        self.trajectory_solvent = md.load(trajectoryfile)
        self.trajectory = self.trajectory_solvent.remove_solvent()

        # Load reference that was the input structure for MD simulation
        # Compare and prepend this structure to the beginning of the trajectory
        self.ref = md.load_frame(referencefile, 0)
        self.check_ref_match()

        # Remove unitcell to combine trajectories
        self.trajectory.unitcell_angles = None 
        self.trajectory.unitcell_lengths = None
        self.trajectory = self.trajectory.join(self.ref)

        # Lists containing atom IDs
        self.top = self.trajectory.topology
        self.protein = self.top.select("chainid == 0")
        self.peptide = self.top.select("chainid == 1")

        # Create a list of the IDs for all atoms of the peptide chain
        self.peptide_atomids = [atomid for atomid in self.peptide]

    def check_ref_match(self):
        if self.trajectory.n_chains != 2:
            print(f"{self.trajectory.n_chains} chains in trajectory.")
            raise AssertionError("Number of chains in trajectory is not 2.")
        if self.trajectory.n_chains != self.ref.n_chains:
            print(f"{self.trajectory.n_chains} chains in trajectory.\n{self.ref.n_chains} chains in reference.")
            raise AssertionError("Number of chains in trajectory and reference differ.")
        if self.trajectory.n_residues != self.ref.n_residues:
            print(f"{self.trajectory.n_residues} residues in trajectory.\n{self.ref.n_residues} residues in reference.")
            raise AssertionError("Number of residues in trajectory and reference differ.")
        if self.trajectory.n_atoms != self.ref.n_atoms:
            print(f"{self.trajectory.n_atoms} atoms in trajectory.\n{self.ref.n_atoms} atoms in reference.")
            raise AssertionError("Number of atoms in trajectory and reference differ.")


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
        self.ref = md.load_frame(referencefile, 0)
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
    mytraj = md.load(args.input).remove_solvent()
    myref = md.load(args.ref)
    newtraj = prepend_reference(mytraj, myref)
