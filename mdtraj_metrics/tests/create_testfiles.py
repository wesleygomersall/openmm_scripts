#!/usr/bin/env python3

from tools import * 

if __name__ == "__main__":
    mytraj = md.load(args.input).remove_solvent()
    reference = md.load(args.ref)

    truncated_test = mytraj[:10]

    truncated_test.save_pdb('test_trajectory.pdb', force_overwrite=True)
    reference.save_pdb('test_reference.pdb', force_overwrite=True)
