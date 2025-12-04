# mdtraj_metrics

This is intended to be a bit more universal analysis of trajectories. 

The basic idea of this one will be to calculate similar metrics to the original
version, `mdtraj_metrics`, but which will work with a universal
`trajectory.pdb` input and the following parameters: 
- (OPTIONAL) Reference `input.pdb` which shall be considered as the starting
  point for the system's movement. If this is not provided, then it is assumed
  that the first state or frame of the input trajectory should be used as the
  reference starting point for these calculations. 
- Chain Info File. Ouput statistics are given on a per-chain basis. Chain
  information is given in a file with format given in example file below. 
  `chinfo` is the default name that the scripts will look for if
  none is given, otherwise the user may specify the file path. This file must
  be in the following format: A line containing a `@Chainidentifier`, followed
  by a line containing integers separated by white space. These integers
  correspond to the residue numbers of the residues/bases of that chain. If
  there is a skip between integers in this list then the skip is ignored.
  Therefore `1 2 5` is equivalent to `1 2 3 4 5`. Any overlap between chains is
  not permitted. In the absence of any file, then the entire system is analyzed
  as a single chain given the degault name "Chain0".

```
$ cat chinfo
@Chain1: # Arbitrary chain name, will be used in the output file path.
1 2 3 4 5 6 7 ... 200
@Chain2: 
201 202 203 204 205
```

There is an example in tests/test_chaininfo. 

Currently not compatible are backbone RMSD and selections including non-protein
atoms (tRNA in my case).

## OLD/ 

This directory contains the old mdtraj metrics scripts. Most of these scripts
rely on prepending the reference pdb to the trajectory and subsequently
calculating metrics using the 0th frame as the reference point. There are also
scripts which are intended for use only with trajectories which have periodic
boundary conditions (not currently implemented in openMM scripts).
