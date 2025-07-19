#!/usr/bin/env python3

from tools import *

if __name__ == "__main__":
    mytraj = SmdTrajectory(args.input)
    myref = Reference(args.ref, mytraj)

    # Calculate RMSD for each frame, using the 0th (from ref) as reference
    peptide_rmsd = md.rmsd(mytraj.trajectory, myref.ref, 0, mytraj.peptide_atomids)

    # Create csv with columns: Frame, RMSD
    with open(args.output.split('.csv')[0] + ".csv", 'w') as out: 
        out.write("Frame,RMSD(nm)\n")
        for i, frame_rmsd in enumerate(peptide_rmsd): 
            out.write(f"{i+1},{frame_rmsd}\n")
