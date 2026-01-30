#!/usr/bin/env python3

import numpy as np
import pandas as pd
import mdtraj as md 
import collect_in

def displacement(trajectory, reference, chains): 
    # compute all Ca-Ca distances between chains
    pairs = trajectory.topology.select_pairs(selection1='chainid 0 and name CA', 
                                             selection2='chainid 1 and name CA')

    displacement = md.compute_distances(trajectory, pairs)
    row_means = np.mean(displacement, axis=1)
    df = pd.DataFrame(row_means, columns=['Ligand Displacement'])
    return df
    
def main():
    data = collect_in.get_traj_inputs()
    displacements = displacement(data.input_traj, data.reference, data.chains)
    outfilepath = collect_in.get_output_path(data.input_traj_filepath, 
                                             "_displacement", 
                                             ".csv")
    displacements.to_csv(outfilepath)

if __name__ == "__main__":
    main()
