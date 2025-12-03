#!/usr/bin/env python3

import numpy as np
import pandas as pd
import mdtraj as md 
import collect_in

def displacement(trajectory, reference, chains): 
    '''
    Track the displacement between chains center of masses throughout the
    states of a trajectory.
    '''
    df = pd.DataFrame()
    for chain_name, residues in chains: 
        if len(residues) > 1: 
            selection_string = f"index {residues[0] - 1} to {residues[-1] - 1}"
        else: 
            selection_string = f"index {residues[0] - 1}"

        com = md.compute_center_of_mass(trajectory, select=selection_string)
        ref_com = md.compute_center_of_mass(reference, select=selection_string)
        displacement = [np.linalg.norm(np.subtract(c, ref_com)) for c in com]

        col_name = chain_name + "_displacement(nm)"
        df[col_name] = displacement

    return df
    
def main():
    data = collect_in.get_traj_inputs()
    displacements = displacement(data.input_traj, data.reference, data.chains)
    outfilepath = data.input_traj_filepath.strip('.pdb') + "_displacement.csv"
    displacements.to_csv(outfilepath)

if __name__ == "__main__":
    main()
