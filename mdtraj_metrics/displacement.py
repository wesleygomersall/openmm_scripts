#!/usr/bin/env python3

from tools import *

if __name__ == "__main__":
    mytraj = SmdTrajectory(args.input)
    myref = Reference(args.ref, mytraj)

    peptide_displacement = []

    # Create csv with columns: Frame, Displacement
    with open(args.output.split('.csv')[0] + ".csv", 'w') as out: 
        out.write("Frame,Displacement(nm)\n")
        for i, d in enumerate(peptide_displacement): 
            out.write(f"{i+1},{d}\n")
