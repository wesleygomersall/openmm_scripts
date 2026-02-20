#!/usr/bin/env python3

import pandas as pd
import mdtraj as md 
import collect_in
import itertools
import numpy as np
import progressbar

def contacts(trajectory, cutoff: float = 0.35, countatoms: bool = True):
    '''
    Inputs: 
    trajectory:     MDTraj Trajectory object. Must have at least two chains. 
                    Pairs will be found between the first two chains of the 
                    trajectory.
    cutoff:         Maximum distance in Angrstrom to return pairs.
                    Default 3.5 A.
    Output:
                    Pandas dataframe of atom/residue pairs and their distances. 
    '''

    atom_contact_df = pd.DataFrame(columns=['Frame', 
                                            'chainA_atomID', 
                                            'chainA_name', 
                                            'chainA_resn', 
                                            'chainA_resi', 
                                            'chainB_atomID', 
                                            'chainB_name', 
                                            'chainB_resn', 
                                            'chainB_resi', 
                                            'distance', 
                                            ])

    chA_res_indices = [residue.index for residue in trajectory.topology.chain(0).residues]
    chB_res_indices = [residue.index for residue in trajectory.topology.chain(1).residues]

    all_pairs = list(itertools.product(chA_res_indices, chB_res_indices))
    distances, _ = md.compute_contacts(trajectory, all_pairs) 
    # distancesnp.ndarray, shape=(n_frames, n_pairs), dtype=np.float32

    bar = progressbar.ProgressBar(maxval=trajectory.n_frames)
    bar.start()

    for framenum, frame in enumerate(trajectory):
        bar.update(framenum + 1) # Progress updates

        frame_contacts = atom_contact_df.iloc[0:0].copy()

        for index, d in enumerate(distances[framenum,]):
            if d < cutoff:

                # For pairs within cutoff, loop through atoms in chainA res
                for atom in trajectory.topology.select(f"resid {all_pairs[index][0]}"):
                    chA_atom = trajectory.topology.atom(atom)
                    a_index = chA_atom.index
                    a_name = chA_atom.name
                    a_elem = chA_atom.element.name
                    a_res = chA_atom.residue
                    a_resnum = chA_atom.residue.index

                    if a_elem == 'hydrogen': continue 

                    haystack = [atomid for atomid in trajectory.topology.select(f"resid {all_pairs[index][1]}")] 
                    neighbors = md.compute_neighbors(frame, cutoff, np.array([atom]), haystack, periodic=True)


                    for i in neighbors[0]: 
                        if i != atom and trajectory.topology.atom(i).element.name != 'hydrogen': 
                            chB_atom = trajectory.topology.atom(i)
                            frame_contacts = pd.concat([frame_contacts if not frame_contacts.empty else None, 
                                                        pd.DataFrame([{'Frame':framenum, 
                                                                        'chainA_atomID': a_index, 
                                                                        'chainA_name': a_name,
                                                                        'chainA_elem': a_elem,
                                                                        'chainA_resn': a_res, 
                                                                        'chainA_resi': a_resnum, 
                                                                        'chainB_atomID': chB_atom.index, 
                                                                        'chainB_name': chB_atom.name, 
                                                                        'chainB_elem': chB_atom.element.name,
                                                                        'chainB_resn': chB_atom.residue, 
                                                                        'chainB_resi': chB_atom.residue.index, 
                                                                        'distance': d}])], 
                                                        ignore_index=True)
        atom_contact_df = pd.concat([atom_contact_df if not atom_contact_df.empty else None, 
                                     frame_contacts], 
                                    ignore_index=True)

    bar.finish()
    return atom_contact_df

def atoms_to_res_contacts(atom_contact_df):

    res_contact_df = pd.DataFrame(columns=['Frame', 
                                           'chainA_resn', 
                                           'chainB_resn', 
                                           'distance', 
                                           ])
    
    maxframe = atom_contact_df['Frame'].max()
    for i in range(maxframe + 1):
        temp_df = atom_contact_df.loc[atom_contact_df['Frame'] == i]
        temp_df['res_pair'] = list(zip(temp_df['chainA_resn'], temp_df['chainB_resn']))
        unique_pairs = set(temp_df['res_pair'])

        for pair in unique_pairs:
            min_dist = atom_contact_df.loc[atom_contact_df['res_pair'] == pair]['distance'].min()

            frame_contacts = pd.concat([res_contact_df if not res_contact_df.empty else None, 
                                        pd.DataFrame([{'Frame': i, 
                                                       'chainA_resn': pair[0], 
                                                       'chainB_resn': pair[1], 
                                                       'distance': min_dist, }])], 
                                       ignore_index=True)

    return res_contact_df

def main():
    data = collect_in.get_traj_inputs()
    atom_df = contacts(data.input_traj, 0.35)

    atom_df.to_csv(collect_in.get_output_path(data.input_traj_filepath, 
                                             "_contacts_byatom", 
                                             ".csv"))

    res_df = atoms_to_res_contacts(atom_df)
    res_df.to_csv(collect_in.get_output_path(data.input_traj_filepath, 
                                            "_contacts_byres", 
                                            ".csv"))

if __name__ == "__main__":
    main()
