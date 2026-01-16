#!/usr/bin/env python3

import argparse
import mdtraj as md 

def save_without_solv(trajectory_file_path): 
    '''
    Save version of trajectory without any solvent. 
    '''
    my_traj = md.load(trajectory_file_path).remove_solvent()
    if trajectory_file_path.endswith('.pdb'):
        outfilepath = trajectory_file_path.strip('.pdb') + "_solventRemoved.pdb"
    elif trajectory_file_path.endswith('.pdb.gz'):
        outfilepath = trajectory_file_path.strip('.pdb.gz') + "_solventRemoved.pdb"
    else: 
        raise TypeError("{trajectory_file_path} file suffix not recognized. Use .pdb or .pdb.gz file.")
    my_traj.save_pdb(outfilepath, force_overwrite = True) 
    return 0 

def main(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, help="Input trajectory.pdb.")
    args = parser.parse_args()

    save_without_solv(args.input)

if __name__ == "__main__":
    main()
