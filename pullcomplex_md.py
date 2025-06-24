#!/usr/bin/env python3

from openmm.app import *
from openmm import *
from openmm.unit import *

from sys import stdout
import argparse
import numpy as np

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
parser.add_argument("--add-water", action="store_true", help="Add this option to populate modeller with a surrounding box of water molecules.")
parser.add_argument("--no-protein-move", action="store_true", help="Add this option to hold all atoms in the protein chain in place throughout the simulation.")
args = parser.parse_args()

def center_of_mass(pos, system, atoms): 
    '''
    Find the center of mass for a group of atoms in topology.
    This should not be used for groups of atoms with altered masses (mass set
    to zero to suppress movement). 

    Input(s): 
        pos:                    Positions within topology
        system:                 System object from openMM
        atoms (integer array):  Array of atom indices to compute center of mass 

    Output(s): 
        Center of mass (numpy array of euclidian coordinates).
    '''
    masses = np.array([system.getParticleMass(i.index) for i in atoms])
    center_of_mass = sum([m * pos[i.index] for i, m in zip(atoms, masses)]) / masses.sum()
    return center_of_mass

def custom_force(atoms, force_constant):
    '''
    Define custom force pulling on a group of atoms

    Input(s): 
        atoms (str):            List of atoms ids to apply force to.
        f_constant (float):     Constant force * kcal/mole/A

    Output(s): 
        openmm.CustomCentroidBondForce object
    '''
    force_constant *= kilocalories_per_mole / angstroms

    fx = 0.0
    fy = 0.0
    fz = -1.0 # apply force in the z direction only

    # pull = CustomCentroidBondForce(2, "-1*force_constant*distance(g1, g2)")
    pull = CustomCentroidBondForce(1, "force_constant * (x1*fx + y1*fy + z1*fz)")

    # pull.addPerBondParameter("force_constant", force_constant)
    pull.addGlobalParameter("force_constant", force_constant)
    # pull.addPerBondParameter("fx", fx)
    pull.addGlobalParameter("fx", fx)
    # pull.addPerBondParameter("fy", fy)
    pull.addGlobalParameter("fy", fy)
    # pull.addPerBondParameter("fz", fz)
    pull.addGlobalParameter("fz", fz)

    # pull.addGroup(g1)
    # pull.addGroup(g2)
    pull.addGroup(atoms)
    pull.addBond([0])

    return pull

tstep = args.timestep * femtoseconds
# Noora suggested 1.0*picoseconds step size
pressure = args.pressure * atmosphere
temperature = args.temp * kelvin
freq = 1000 # store data `freq` times depending on total steps 
store = args.steps // freq
output_log_path = args.output + "/log.txt"
output_pdb_path = args.output + "/output.pdb"
output_energy_stats = args.output + "/stats.csv"
output_displacement = args.output + "/displacement.csv"

log = open(output_log_path, 'w') # write run info to this log rather than stdout

if not os.path.exists(args.output):
    log.write(f"making directory: {args.output}\n")
    os.mkdir(args.output)

log.write(f"Reading input pdb {args.input}\n")
pdb = PDBFile(args.input)

forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)

if len(list(modeller.topology.chains())) != 2:
    # take care of possibly unbound protein and peptide list of residues
    raise Exception("Topology does not have exactly two chains.")
    
for i, c in enumerate(modeller.topology.chains()):
    if i == 0: protein = [res for res in c.residues()]
    else: peptide = [res for res in c.residues()]

# optionally change mass of these atoms to 0 to suppress movement
if args.no_protein_move: 
    log.write("Setting mass of the first chain (protein chain) to zero to suppress motion.\n")
    for residue in protein:
        for atom in residue.atoms():
            system.setParticleMass(int(atom.id), 0)

log.write("Adding hydrogen atoms.\n")
modeller.addHydrogens(forcefield)

# Optionally add water as solvent
if args.add_water: 
    log.write("Adding solvent: water\n")
    modeller.addSolvent(forcefield, model='tip3p', 
                        positiveIon='Na+', negativeIon='Cl-',
                        padding= 1 * nanometer, 
                        neutralize=True)

# create system needs to be after adding solvent and hydrogens
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=app.NoCutoff, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds)

# Add custom force pulling on the peptide to the system
peptide_atoms = [atom.index for residue in peptide for atom in residue.atoms()]
pullforce = custom_force(peptide_atoms, 25)
system.addForce(pullforce)

integrator = LangevinIntegrator(temperature, 1/picosecond, tstep) 

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

log.write("Calculating initial center of mass location.\n")
t0_peptide_center = center_of_mass(
        modeller.positions,
        system,
        [a for i in peptide for a in i.atoms()])

log.write("Calculating initial potential energy\n")
t0_sim_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()

log.write(f"Temperature: {temperature}\n")
log.write(f"Time step: {tstep}\n")
log.write(f"Initial peptide center of mass: {t0_peptide_center}\n")

log.write(f"Initial potential energy: {t0_sim_energy}\n")
simulation.minimizeEnergy()
minimized_sim_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
log.write(f"Potential energy after minimization: {minimized_sim_energy}\n")

log.write(f"Running sumulation, storing data every {store} iterations.\n")

# Andrew imports a custom reporter from md_helper.py 
# but it looks like he prints displacement to stdout for SMD.

simulation.reporters.append(PDBReporter(output_pdb_path, store))
simulation.reporters.append(StateDataReporter(output_energy_stats, 
                                              store, 
                                              step=True, 
                                              time=True,
                                              potentialEnergy=True, 
                                              kineticEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              elapsedTime=True))

with open(output_displacement, "w") as dispout: 
    dispout.write("Step,Distance(nm),PeptideCenterCoordinate\n")
    for step in range(args.steps):
        simulation.step(1)
        peptide_center = center_of_mass(
                simulation.context.getState(getPositions=True).getPositions(asNumpy=True),
                simulation.context.getSystem(),
                [a for i in peptide for a in i.atoms()])
        # distance = peptide_center - t0_peptide_center
        distance = np.linalg.norm(peptide_center - t0_peptide_center)
        if (step + 1) % store == 0: 
            dispout.write(f"{step+1},{distance},{peptide_center}\n")

# Before using the above loop I used this to advance the steps 
# simulation.step(args.steps)
