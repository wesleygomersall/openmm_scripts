#!/usr/bin/env python3

from tools import *
from rmsd import peptide_rmsd
from displacement import peptide_displacement
from intralength import acarbon_distances

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    myref = md.load(args.ref)
    disp = peptide_displacement(mytraj, myref)

    prepended = prepend_reference(mytraj, myref)
    prmsd = peptide_rmsd(prepended)
    acdist = acarbon_distances(prepended)
    
    assert len(disp) == len(prmsd) == len(acdist)

    for i in range(prepended.n_frames):
        if i == 0: args.output.write("Frame,Displacement(nm),RMSD(A),aCarbon-aCarbon Distance(nm)\n")
        args.output.write(f"{i},{disp[i]},{prmsd[i]},{acdist[i]}\n")
    
    args.output.close()
