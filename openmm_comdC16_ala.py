#!usr/bin/env python3

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

import argparse

# parameters tuned for comparison
steps = 1000000
store = steps / 100
k = 300

# pdb = PDBFile('af3_c16.pdb')
pdb = PDBFile('comd_c16_ala_af3multimer_pdbfixed.pdb')
forcefield = ForceField('amber14/protein.ff14SB.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds) # Noora suggested 1.0*picoseconds step size
# integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, 0.001) # Variable step sizes

# Suppress movement of the protein but not the peptide
# Mass =0 suppressed atom movement

protein = range(3865)
peptide = range(3866, 4027)
# Use CustomCentroidBondForce to add forces between centers of mass
# Or use CustomExternalForce? 
repulse = CustomCentroidBondForce(2, "-1*k*distance(g1, g2)")
repulse.addPerBondParameter("k")
repulse.addGroup(protein)
repulse.addGroup(peptide)
repulse.addBond([0, 1], [k])

system.addForce(repulse)

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('testing_output_c16_ala.pdb', store))
simulation.reporters.append(StateDataReporter(stdout, store, step=True, potentialEnergy=True, temperature=True))
simulation.step(steps)
