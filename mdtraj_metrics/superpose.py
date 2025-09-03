#!/usr/bin/env python3

from tools import *
# import numpy as np
# import matplotlib.pyplot as plt

def save_superpose(trajectory, reference, fileout):
    '''
    Saves a version of the trajectory after superposing protein chain A to the
    first frame of reference.
    
    Input: 
        mdtraj Trajectory objects trajectory and reference
        fileout (str): output file path

    Output: 
        Write new trajectory {original_name}_superposed.pdb
    '''

    prot_atom_indices = trajectory.topology.select("chainid == 0")
    filename = "output_superposed.pdb"
    superposed = trajectory.superpose(reference, 0, prot_atom_indices) 
    superposed.save_pdb(fileout)
    return

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    myref = md.load(args.ref) 

    filename = args.input.strip(".pdb") + "_cleanedpbc.pdb" 
    save_superpose(mytraj, myref, filename)

    args.output.close()
