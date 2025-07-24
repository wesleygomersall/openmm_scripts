#!/usr/bin/env python3

from tools import *

def peptide_rmsd(trajectory): 
    '''
    Calculate backbone RMSD of peptide chain ("chainid == 1") for each frame of
    trajectory, using the first frame as reference. 
    
    Input: 
        mdtraj Trajectory object
    Output: 
        List: RMSD values, length = trajectory.n_frames
    '''
    peptide_bb_atomids = trajectory.topology.select("(chainid == 1) and backbone")
    peprmsd = md.rmsd(trajectory, trajectory, 0, peptide_bb_atomids)
    return peprmsd

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    mytraj = prepend_reference(mytraj, md.load(args.ref))
    
    prmsd = peptide_rmsd(mytraj)

    for i, p in enumerate(prmsd):
        if i == 0: args.output.write("Frame,RMSD(nm)\n")
        args.output.write(f"{i},{p}\n")
    
    args.output.close()
