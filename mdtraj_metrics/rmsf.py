#!/usr/bin/env python3

from tools import *
import numpy as np
import matplotlib as plt

def rmsf(trajectory):
    '''
    Calculate per-residue RMSF for each chain in a trajectory, first frame as
    reference. 
    
    Input: 
        mdtraj Trajectory object 
    Output: 
        multiple (N=trajectory.getNumChains) RMSF lists
    '''
    rmsf = md.rmsf(trajectory, trajectory, 0) 
    rmsf = 10 * rmsf # nm -> A

    aa_rmsf = []
    for res in trajectory.topology.residues:
        indices=[atom.index for atom in res.atoms]
        resi_rmsf_values=rmsf[indices]
        aa_rmsf.append(np.mean(resi_rmsf_values))

    counter = 0 
    rmsf_by_chain = []
    for chain in trajectory.topology.chains: 
        templist = []
        for residue in chain.residues: 
            templist.append(aa_rmsf[counter])
            counter += 1
        rmsf_by_chain.append(templist) 

    # testing, remove prints() later
    for chainrmsf in rmsf_by_chain: 
        print(chainrmsf)
        for aarmsf in chainrmsf:
            print(aarmsf)

    return rmsf_by_chain

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    mytraj = prepend_reference(mytraj, md.load(args.ref))
    
    my_rmsfs = rmsf(mytraj)

    chain_no = 0 
    for chain_rmsf in my_rmsfs:
        filename = f"rmsf_chain{chain_no}.csv"
        with open(filename, 'w') as fout:
            fout.write("Residue,RMSF(A)\n")
            for i, rmsf in enumerate(chain_rmsf): 
                fout.write(f"{i},{rmsf}\n")

        # plot 
        plt.clf() # clear figure from any previous plotting
        plt.plot(range(1, len(chain_rmsf) + 1), chain_rmsf)
        plt.xlabel('Residue')
        plt.ylabel('RMSF (A)')
        plt.title('Per residue RMSF')
        plt.grid(False)
        output_path = f"RMSF_chain{chain_no}.pdf" 
        plt.savefig(output_path, format="pdf")

        chain_no += 1

    args.output.close()
