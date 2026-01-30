#!/usr/bin/env python3

import pandas as pd
import numpy as np
import mdtraj as md 
import collect_in

def calc_rmsf(trajectory, ref):
    '''
    Calculate per-residue RMSF for each chain in a trajectory, first frame as
    reference. 
    
    Input: 
        mdtraj Trajectory object 
        reference pdb 
    Output: 
        multiple (N=trajectory.getNumChains) RMSF lists
    '''

    rmsf_nm = md.rmsf(trajectory, ref, 0) 
    rmsf = [item * 10.0 for item in rmsf_nm] # nm -> A

    aa_rmsf = []
    for res in trajectory.topology.residues:
        indices=[atom.index for atom in res.atoms]
        resi_rmsf_values= [rmsf[i] for i in indices]
        aa_rmsf.append(np.mean(resi_rmsf_values))

    counter = 0 
    rmsf_by_chain = []
    for chain in trajectory.topology.chains: 
        templist = []
        for residue in chain.residues: 
            templist.append(aa_rmsf[counter])
            counter += 1
        rmsf_by_chain.append(templist) 

    return rmsf_by_chain

def main(): 
    data = collect_in.get_traj_inputs()
    rmsfs = calc_rmsf(data.input_traj, data.reference) # return list of lists

    for chainnum, chainrmsd in enumerate(rmsfs):
        df = pd.DataFrame({'RSMF(A)':chainrmsd})
        df.to_csv(collect_in.get_output_path(data.input_traj_filepath, 
                                             f"_rmsf_ch{chainnum}", 
                                             ".csv")
                  )

if __name__ == "__main__":
    main()
