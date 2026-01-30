#!/usr/bin/env python3

import pandas as pd
import mdtraj as md 
import collect_in

def delta_rmsd(trajectory, reference, chains, bb_only = False): 
    '''
    Return a data frame with RMSD as separate columns corresponding to the atom
    groupings in input variable chains from collect_in.Traject.chains.
    Delta RMSD is the RMSD between each subsequent frame. The first frame of
    the simulation utilizes the reference pdb input (there is no previous
    frame).
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

        chain_rmsd = []
        for framenum, frame in enumerate(trajectory):
            if framenum == 0: # first frame 
                frame_drmsd = md.rmsd(frame, 
                                      reference, 
                                      0, 
                                      chain_selection) 
            else: 
                frame_drmsd = md.rmsd(frame, 
                                      trajectory,
                                      (framenum - 1), 
                                      chain_selection)
            chain_rmsd.append(frame_drmsd[0]) 

        # (output in A)
        rmsd_A = [item * 10.0 for item in chain_rmsd] 
        col_name = chain_name + "_RMSD(A)"
        df[col_name] = rmsd_A

    return df

def main(): 
    data = collect_in.get_traj_inputs()
    rmsds = delta_rmsd(data.input_traj, data.reference, data.chains)
    rmsds.to_csv(collect_in.get_output_path(data.input_traj_filepath, 
                                            "_deltaRMSDs", 
                                            ".csv"))

if __name__ == "__main__":
    main()
