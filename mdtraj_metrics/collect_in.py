#!/usr/bin/env python3

import argparse
import mdtraj as md 
import os
import sys

class Traject:
    '''
    self.input_traj (mdtraj trajectory object) 
    self.input_traj_filepath (/path/to/trajectory.pdb) 
    self.input_dcd_filepath (/path/to/trajectory.dcd)
    self.reference (mdtraj single frame pose)
    self.chains (list of tuples: (chain name, atoms in chain) )
    '''

    def __init__(self, in_trajectory, reference_file, chain_file, dcd_path: None):

        self.input_traj_filepath = in_trajectory
        self.input_dcd_filepath = dcd_path

        # Get input_traj, make sure exists and is more than one frame 
        if in_trajectory is None or not os.path.exists(in_trajectory):
            print(f"File {in_trajectory} was not found. Exiting now.")
            sys.exit(1)
        if dcd_path is None or not os.path.exists(dcd_path):
            # Just load the pdb as trajectory
            self.input_traj = md.load(in_trajectory).remove_solvent()
        else:
            self.input_traj = md.load_dcd(dcd_path, top=in_trajectory)
        assert self.input_traj.n_frames > 1

        # Reference is first frame of input trajectory if no reference file is provided
        self.reference = self.input_traj[0] if reference_file == "" else md.load(reference_file).remove_solvent()

        # chain info from chain file. default: all particles/residues in one chain
        self.chains = []
        if chain_file is None or not os.path.exists(chain_file):
            print(f"Chain info file was not found. Output will not be separated by chains. This may cause error with identifying anchors for peroidic boundary.")
            self.chains.append(("Chain0", range(1, self.input_traj.n_atoms + 1)))
        else: # read chinfo file into chains list
            with open(chain_file, 'r') as fin: 
                name = "initial value, shouldn't be saved"
                linenumber = 1
                while True: 
                    line = fin.readline()
                    if not line: break
                    if linenumber == 1 and line[0] != '@': 
                        raise ValueError("Unexpected first line of file")
                    if line == '': continue
                    if line[0] == '@': 
                        name = line.strip('@ \n') 
                    else: 
                        atom_range = range(int(line.split()[0]), int(line.split()[-1]) + 1)
                        if name is None: pass
                        else: self.chains.append((name, atom_range))
                        name = "reset name, shouldn't be saved!"
                    linenumber += 1
        self.check_chain_overlap()

        # For periodic boundary
        # TODO: Currently can break if no chain info file specified
        anchors = []
        for index in self.chains[0][1]: # anchor is first chain in chain info file
            anchors.append(self.input_traj.topology.atom(index))
        self.input_traj.image_molecules(inplace = True, 
                                        anchor_molecules = [set(anchors)])

    def check_chain_overlap(self):
        # ensure no overlap between chains
        allselected = set()
        for c in self.chains: 
            thischain = set(c[1]) 
            if not allselected.isdisjoint(thischain):
                raise AssertionError("Chains overlap!")
            allselected = allselected.union(thischain)

def get_traj_inputs(): 
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, help="Input trajectory.pdb.")
    parser.add_argument("--dcd", "-d", type=str, default = None, help="Input trajectory.dcd.")
    parser.add_argument("--ref", "-r", type=str, default = "", help="Reference input.pdb to be used as t0 state. If none is specified, then first state of input trajectory.pdb is used as t0 state.")
    parser.add_argument("--chains", "-c", type=str, default = "", help="Chain Info File. See README.md")
    args = parser.parse_args()

    my_traj = Traject(args.input, args.ref, args.chains, args.dcd)

    return my_traj

def get_output_path(in_path: str, 
                    new_label: str,
                    new_ext: str,
                    possible_exts: list = ['.pdb', '.pdb.gz']):

    for ext in possible_exts: 
        if in_path.endswith(ext):
            return in_path.strip(ext) + new_label + new_ext
    raise TypeError("{in_path} file suffix not recognized. Use .pdb or .pdb.gz file.")
    return 1
