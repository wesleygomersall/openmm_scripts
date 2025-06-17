#!/usr/bin/env python3

from openmm.app import *
from openmm import *
from openmm.unit import picosecond, femtoseconds, atmosphere, kelvin, nanometer
# from openmm.unit import *
# from sys import stdout

import argparse

# previously imported from params
steps = 1_000_000
store = steps / 100
k = 300

parser = argparse.ArgumentParser(description="MD trajectory simulating constant force pulling apart a protein-peptide complex.")
parser.add_argument("--input", type=str, help="Path to input pdb containing both protein and peptide. If errors exist, try using PDBfixer first.")
parser.add_argument("--output", type=str, help="Path to output directory. Will be created if not already existing.")
parser.add_argument("--timestep", type=float, default=2.0, help="Time between steps in femtoseconds, default is 2.")
parser.add_argument("--steps", type=int, default=10_000, help="Total steps in simulation, default is ten thousand.")
parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (atmospheres), default is 1 atm.")
parser.add_argument("--temp", type=float, default=300.0, help="Temperature (Kelvin), default is 300K.")
args = parser.parse_args()

if not os.path.exists(args.output):
    os.mkdir(args.output)

tstep = args.timestep * femtoseconds
# Noora suggested 1.0*picoseconds step size
pressure = args.pressure * atmosphere
temperature = args.temp * kelvin
pdb = PDBFile(args.input)

# forcefield = ForceField('amber14/protein.ff14SB.xml')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)

# residue_all = defaultdict(list)
# for i, c in enumerate(modeller.topology.chains()):
    # residue_all[i].extend([r for r in c.residues()])








# modeller.addHydrogens(forcefield)

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(temperature, 1/picosecond, tstep) 
# integrator = VariableLangevinIntegrator(300*kelvin, 1/picosecond, 0.001) # Variable step sizes


# Suppress movement of the protein but not the peptide
# Mass =0 suppressed atom movement







"""

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

"""
