#!/usr/bin/env python3

import argparse
import mdtraj as md 

def save_without_solv(trajectory_file_path): 
    '''
    Save version of trajectory without any solvent. 
    '''
    my_traj = md.load(trajectory_file_path).remove_solvent()
    outfilepath = trajectory_file_path.strip('.pdb') + "_solventRemoved.pdb"
    my_traj.save_pdb(outfilepath, force_overwrite = True) 
    return 0 

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, help="Input trajectory.pdb.")
    args = parser.parse_args()

    save_without_solv(args.input)

if __name__ == "__main__":
    main()
