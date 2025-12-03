#!/usr/bin/env python3

import mdtraj as md
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", "-i", type=str, help="")
    parser.add_argument("--ref", "-r", type=str, help="")
    args = parser.parse_args()

    mytraj = md.load(args.input).remove_solvent()
    reference = md.load(args.ref)

    truncated_test = mytraj[:10]

    truncated_test.save_pdb('test_trajectory.pdb', force_overwrite=True)
    reference.save_pdb('test_reference.pdb', force_overwrite=True)
