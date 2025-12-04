#!/usr/bin/env python3

import mdtraj as md 
import numpy as np
from tools import *

def output_clean_pdb(trajectorypath): 
    '''
    Input: 
        trajectorypath:     path to `.pdb` output trajectory from MD with
                            periodic boundary condition and explicit solvent.
    Output: 
        writes a new file with chainid 0 and 1 centered.
    '''
    
    trajectory = md.load(trajectorypath).remove_solvent()
    centeredtraj = center_traj_pbc(trajectory)
    basename = trajectorypath.strip(".pdb")
    filename = basename + "_cleanedpbc.pdb" 
    centeredtraj.save_pdb(filename)
    return

if __name__ == "__main__":
    output_clean_pdb(args.input)
