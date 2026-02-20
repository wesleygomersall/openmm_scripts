#!/usr/bin/env python3

import pandas as pd
import mdtraj as md 
import collect_in

def delta_rmsd(trajectory, reference, chains, bb_only = False): 
    '''
    Return a data frame with delta RMSD as separate columns corresponding to
    the chains in trajectory. Delta RMSD is the RMSD between each subsequent
    frame. The first frame of the simulation has delta RMSD of 0. (there is no
    previous frame). Default to backbone RMSD of the chain(s). 
    '''
    df = pd.DataFrame()

    for chain_num, chain in enumerate(trajectory.topology.chains):
        indices = trajectory.topology.select(f"chainid {chain_num}")

        chain_drmsd = []
        for framenum, frame in enumerate(trajectory):
            if framenum == 0: 
                chain_drmsd.append(0) # no previous frame
            else: 
                chain_drmsd.append(md.rmsd(frame, 
                                      trajectory,
                                      (framenum - 1), 
                                      indices)[0]) 

        drmsd_A = [item * 10.0 for item in chain_drmsd] 
        col_name = f"chain{chain_num}_RMSD(A)"
        df[col_name] = drmsd_A

    return df

def main(): 
    data = collect_in.get_traj_inputs()
    rmsds = delta_rmsd(data.input_traj, data.reference, data.chains)
    rmsds.to_csv(collect_in.get_output_path(data.input_traj_filepath, 
                                            "_deltaRMSDs", 
                                            ".csv"))

if __name__ == "__main__":
    main()
