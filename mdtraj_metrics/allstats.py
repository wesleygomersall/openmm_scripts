#!/usr/bin/env python3

from tools import *
from rmsd import peptide_rmsd
from displacement import peptide_displacement
from intralength import acarbon_distances

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    mytraj = prepend_reference(mytraj, md.load(args.ref))

    coords, disp = peptide_displacement(mytraj)
    prmsd = peptide_rmsd(mytraj)
    acdist = acarbon_distances(mytraj)

    for i in range(mytraj.n_frames):
        if i == 0: args.output.write("Step,Displacement(nm),RMSD(nm),aCarbon-aCarbon Distance(nm)\n")
        args.output.write(f"{i},{disp[i]},{prmsd[i]},{acdist[i]}\n")
    
    args.output.close()
