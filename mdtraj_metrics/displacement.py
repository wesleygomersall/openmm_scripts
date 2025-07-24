#!/usr/bin/env python3

import mdtraj as md 
import numpy as np
from tools import *

def peptide_displacement(trajectory):
    '''
    Track the displacement of peptide (chainid == 1) throughout a trajectory.
    Displacement is calculated from the first frame of trajectory. 
    
    Input: 
        mdtraj Trajectory object
    Output: 
        tuple: List of center of mass coordinates, List of displacement scalars
    '''
    com = md.compute_center_of_mass(trajectory, select="chainid == 1")
    displacement = [np.linalg.norm(np.subtract(c, com[0])) for c in com]
    return com, displacement

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    mytraj = prepend_reference(mytraj, md.load(args.ref))
    
    coords, disp = peptide_displacement(mytraj)

    for i in range(len(coords)):
        if i == 0: args.output.write("Frame,Displacement(nm),Coordinates\n")
        args.output.write(f"{i},{disp[i]},{coords[i]}\n")

    args.output.close()
