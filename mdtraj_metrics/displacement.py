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
    com = md.compute_center_of_mass(trajectory, select="chainid == 1")
    ref_com = md.compute_center_of_mass(reference, select="chainid == 1")
    displacement = [np.linalg.norm(np.subtract(c, com[0])) for c in com]
    return displacement
    
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
