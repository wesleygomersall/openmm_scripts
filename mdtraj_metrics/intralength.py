#!/usr/bin/env python3

import mdtraj as md 
import numpy as np
from tools import *

def acarbon_distances(trajectory): 
    '''
    ************************WIP************************
    Calculate the distance between the second and the second to last alpha carbons of the peptide chain ("chainid == 1"). 

    Input: 
        mdtraj Trajectory object
    Output: 
        List: Distance values between alpha carbons, length = trajectory.n_frames
    '''

    ac_list = trajectory.topology.select("(chainid == 1) and (name == CA)")
    index1, index2 = ac_list[1], ac_list[-2]
    
    dists = []
    coord1 = trajectory.xyz[:, index1, :]
    coord2 = trajectory.xyz[:, index2, :]
    for frame in range(trajectory.n_frames):
        dists.append(np.linalg.norm(np.subtract(coord1[frame], coord2[frame])))

    return dists

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    mytraj = prepend_reference(mytraj, md.load(args.ref))

    acac_dist = acarbon_distances(mytraj)
    
    for i in range(len(acac_dist)):
        if i == 0: args.output.write("Step,aCarbon-aCarbon Distance(nm)\n")
        args.output.write(f"{i},{acac_dist[i]}\n")

    args.output.close()
