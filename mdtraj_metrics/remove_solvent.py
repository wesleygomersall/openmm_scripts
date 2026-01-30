#!/usr/bin/env python3

import argparse
import mdtraj as md 
import collect_in

def save_without_solv(in_trajectory, outfilepath): 
    '''
    Save version of trajectory without any solvent. 
    '''
    my_traj = in_trajectory.remove_solvent()
    my_traj.save_pdb(outfilepath, force_overwrite = True) 
    return 0 

def main(): 
    data = collect_in.get_traj_inputs()
    out_path = collect_in.get_output_path(data.input_traj_filepath, 
                                          "_Solventless", ".pdb")
    save_without_solv(data.input_traj, out_path)

if __name__ == "__main__":
    main()
