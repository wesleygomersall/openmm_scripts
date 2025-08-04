#!/usr/bin/env python3

import mdtraj as md 
import numpy as np
from tools import *


def peptide_displacement(trajectory, reference):
    '''
    Track the displacement between and atom from ptortein and an atom from
    peptide throughout a trajectory.
    
    Input: 
        mdtraj Trajectory object
    Output: 

    '''
    
    # Setup trajectory from tools 
    check_ref_match(trajectory, reference)
    centeredtraj = center_traj_pbc(trajectory) 

    proteinCAs = centeredtraj.topology.select("(chainid == 0) and (name == CA)")
    peptideCAs = centeredtraj.topology.select("(chainid == 1) and (name == CA)")
    atom_pairs= []
    for atom_a_idx in proteinCAs:
        for atom_b_idx in peptideCAs:
            atom_pairs.append([atom_a_idx, atom_b_idx])
    atom_pairs= np.array(atom_pairs)

    traj_disp = md.compute_displacements(centeredtraj, atom_pairs, periodic=True)
    ref_disp = md.compute_displacements(reference, atom_pairs) 
    traj_disp = np.linalg.norm(traj_disp, axis=2)
    ref_disp = np.linalg.norm(ref_disp, axis=2)
    distances = []

    # First get distance in reference
    distances.append(sum(ref_disp[0]) / len(ref_disp[0]))

    # Append distance from each frame
    for frame in range(centeredtraj.n_frames):
        # print(frame)
        # print(traj_disp[frame])
        distances.append(sum(traj_disp[frame]) / len(traj_disp[frame]))

    return distances


if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    myref = md.load(args.ref)
    
    displacement = peptide_displacement(mytraj, myref)
    print(len(displacement)) 
    print(displacement)

    for i in range(len(displacement)):
        if i == 0: args.output.write("Frame,Displacement(nm)\n")
        args.output.write(f"{i},{displacement[i]}\n")

    args.output.close()
