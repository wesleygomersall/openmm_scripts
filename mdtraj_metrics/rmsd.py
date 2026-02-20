#!/usr/bin/env python3

import pandas as pd
import numpy as np
import mdtraj as md 
import collect_in

def rmsd(trajectory, reference, chains, bb_only = True): 
    '''
    Return a data frame with RMSD as separate columns corresponding to each chain.
    Default to backbone RMSD of the chain(s). 
    '''
    df = pd.DataFrame()
    for chain_num, chain in enumerate(trajectory.topology.chains):
        indices = trajectory.topology.select(f"chainid {chain_num}")

        if bb_only: 
            indices = np.intersect1d(indices, trajectory.topology.select('backbone'))

        chain_rmsd = md.rmsd(trajectory, trajectory, 0, indices)
        chain_rmsd = chain_rmsd * 10 #(output in A)

        col_name = f"chain{chain_num}_RMSD(A)"
        df[col_name] = chain_rmsd

    return df

def main(): 
    data = collect_in.get_traj_inputs()
    rmsds = rmsd(data.input_traj, data.reference, data.chains)
    outfilepath = collect_in.get_output_path(data.input_traj_filepath,
                                             "_RMSDs", ".csv")
    rmsds.to_csv(outfilepath)

if __name__ == "__main__":
    main()
