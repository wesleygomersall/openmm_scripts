#!/usr/bin/env python3

import pandas as pd
import mdtraj as md 
import collect_in

def rmsd(trajectory, reference, chains, bb_only = False): 
    '''
    Return a data frame with RMSD as separate columns corresponding to the atom
    groupings in input variable chains from collect_in.Traject.chains.
    Default to backbone RMSD of the chain(s). 
    '''
    df = pd.DataFrame()
    for chain_name, residues in chains: 
        if len(residues) > 1: 
            selection_string = f"index {residues[0] - 1} to {residues[-1] - 1}"
        else: 
            selection_string = f"index {residues[0] - 1}"

        if bb_only: 
            selection_string = selection_string + " and backbone" 

        chain_selection = trajectory.topology.select(selection_string)
        chain_rmsd = md.rmsd(trajectory, reference, 0, chain_selection)
        chain_rmsd = chain_rmsd * 10 #(output in A)

        col_name = chain_name + "_RMSD(A)"
        df[col_name] = chain_rmsd

    return df

def main(): 
    data = collect_in.get_traj_inputs()
    rmsds = rmsd(data.input_traj, data.reference, data.chains)
    outfilepath = data.input_traj_filepath.strip('.pdb') + "_RMSDs.csv"
    rmsds.to_csv(outfilepath)

if __name__ == "__main__":
    main()
