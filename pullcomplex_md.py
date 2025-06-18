#!/usr/bin/env python3

from openmm.app import *
from openmm import *
from openmm.unit import *
# from sys import stdout

import argparse

# previously imported from params
# steps = 1_000_000
# store = steps / 100
# k = 300 # for force calculations

parser = argparse.ArgumentParser(description="MD trajectory simulating constant force pulling apart a protein-peptide complex.")
parser.add_argument("--input", type=str, help="Path to input pdb containing both protein and peptide. If errors exist, try using PDBfixer first.")
parser.add_argument("--output", type=str, help="Path to output directory. Will be created if not already existing.")
parser.add_argument("--timestep", type=float, default=2.0, help="Time between steps in femtoseconds, default is 2.")
parser.add_argument("--steps", type=int, default=10_000, help="Total steps in simulation, default is ten thousand.")
parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (atmospheres), default is 1 atm.")
parser.add_argument("--temp", type=float, default=300.0, help="Temperature (Kelvin), default is 300K.")
parser.add_argument("--no-protein-move", action="store_true", help="Add this option to hold all atoms in the protein chain in place throughout the simulation.")
args = parser.parse_args()


def custom_force(atoms, force_constant):
    '''
    Defines custom force pulling on a group of atoms
    Input(s): 
        atoms (str):            List of atoms ids to apply force to.
        f_constant (float):     Constant force * kcal/mole/A

    Output(s): 
        openmm.CustomCentroidBondForce object
    '''
    force_constant = force_constant * kilocalories_per_mole / angstroms

    # pull = CustomCentroidBondForce(2, "-1*force_constant*distance(g1, g2)")
    pull = CustomCentroidBondForce(1, "force_constant * (x1*fx + y1*fy + z1*fz")

    # pull.addPerBondParameter("force_constant", force_constant)
    pull.addGlobalParameter("force_constant", force_constant)
    # pull.addPerBondParameter("fx", fx)
    pull.addGlobalParameter("fx", fx)
    # pull.addPerBondParameter("fy", fy)
    pull.addPerBondParameter("fy", fy)
    # pull.addPerBondParameter("fz", fz)
    pull.addGlobalParameter("fz", fz)

    # pull.addGroup(g1)
    # pull.addGroup(g2)
    pull.addGroup(atoms)

    return pull

if not os.path.exists(args.output):
    os.mkdir(args.output)

tstep = args.timestep * femtoseconds
# Noora suggested 1.0*picoseconds step size
pressure = args.pressure * atmosphere
temperature = args.temp * kelvin
pdb = PDBFile(args.input)
print(f"Reading input pdb {args.input}")

forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)

if len(list(modeller.topology.chains())) != 2:
    raise Exception("Topology does not have exactly two chains.")
    
for i, c in enumerate(modeller.topology.chains()):
    if i == 0: protein = [res for res in c.residues()]
    else: peptide = [res for res in c.residues()]

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)

if args.no_protein_move: # optionally change mass of these atoms to 0 to suppress movement
    print("Setting mass of the first chain (protein chain) to zero to suppress motion.")
    for residue in protein:
        for atom in residue.atoms():
            system.setParticleMass(int(atom.id), 0)

modeller.addHydrogens(forcefield)

# Add water 
modeller.addSolvent(forcefield, model='tip3p', 
                    boxSize=None, boxVectors=None, padding=None, numAdded=None, 
                    positiveIon='Na+', negativeIon='Cl-', 
                    ionicStrength=Quantity(value=0, unit=molar), 
                    neutralize=True)

# Add force
peptide_atoms = [atom.index for residue in peptide for atom in residue.atoms()]
pullforce = custom_force(peptide_atoms, 200)
system.addForce(pullforce)

print(f"Temperature: {temperature}")
print(f"Time step: {tstep}")
integrator = LangevinIntegrator(temperature, 1/picosecond, tstep) 

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# get initial potential energy before the simulation begins
# t0_potential = simulation.context.getState(getEnergy=True.getPotentialEnergy()
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('testing_output_c16_ala.pdb', store))
simulation.reporters.append(StateDataReporter(stdout, store, step=True, potentialEnergy=True, temperature=True))
simulation.step(steps)


