#!/usr/bin/env python3

from openmm.app import *
from openmm import *
from openmm.unit import *

from sys import stdout
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="MD trajectory simulating constant force pulling apart a protein-peptide complex.")
parser.add_argument("--input", type=str, help="Path to input pdb containing both protein and peptide. If errors exist, try using PDBfixer first.")
parser.add_argument("--output", type=str, help="Path to output directory. Will be created if not already existing.")
parser.add_argument("--timestep", type=float, default=2.0, help="Time between steps in femtoseconds, default is 2.")
parser.add_argument("--steps", type=int, default=10_000, help="Total steps in simulation, default is ten thousand.")
parser.add_argument("--freq", type=int, default=100, help="Frequency of reporter, will record data `freq` amount of times (every `step // freq` steps).")
parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (atmospheres), default is 1 atm.")
parser.add_argument("--temp", type=float, default=300.0, help="Temperature (Kelvin), default is 300K.")
parser.add_argument("--pull-force", type=float, default=25.0, help="Pull force constant which will be applied to the peptide.")
parser.add_argument("--add-water", action="store_true", help="Add this option to populate modeller with a surrounding box of water molecules.")
parser.add_argument("--suppress-movement", type=str, default='', help="Add this option to hold atoms in the protein chain in place throughout the simulation. If specified, program will look for file specified which contains residue numbers to hold stationary. This file should contain integer IDs separated by lines and/or spaces.") 
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

def custom_force(direction, atoms, force_constant):
    '''
    Define custom force pulling on a group of atoms

    Input(s): 
        direction (list):       (x, y, z) vector for direction of force.
        atoms (str):            List of atoms ids to apply force to.
        f_constant (float):     Constant force * kcal/mole/A

    Output(s): 
        openmm.CustomCentroidBondForce object
    '''
    force_constant *= kilocalories_per_mole / angstroms

    fx = direction[0]
    fy = direction[1]
    fz = direction[2]

    pull = CustomCentroidBondForce(1, "-1 * force_constant * (x1*fx + y1*fy + z1*fz)")
    pull.addGlobalParameter("force_constant", force_constant)
    pull.addGlobalParameter("fx", fx)
    pull.addGlobalParameter("fy", fy)
    pull.addGlobalParameter("fz", fz)
    pull.addGroup(atoms)
    pull.addBond([0])

    return pull

def container_restraint(atoms, radius: float):
    '''
    Define force to contain peptide chain to a spherical container. 

    Input(s): 
        atoms (list):           List of atoms ids to apply force to.
        radius (float):         Radius of container in nm. 

    Output(s): 
        openmm.CustomExternalForce object
    '''
    radius *= nanometer
    force_expression = "container_force*max(0, r - radius)^2; r = sqrt(x^2 + y^2 + z^2) "
    container = CustomExternalForce(force_expression)
    container.addGlobalParameter("radius", radius) 
    container.addGlobalParameter("container_force", 100.0 * kilocalories_per_mole / angstroms) 
    for p in atoms: 
        container.addParticle(p, [])
    return container

# def COM_restraint(atoms, radius):
    # radius *= nanometer
    # force_expression = "COM_restraint_force * max(0, d-rad)^2; d = distance(g1, g2) " 
    # restraint = CustomCentroidBondForce(force_expression) 
    # restraint.addGlobalParameter("rad", radius)
    # restraint.addGlobalParameter("COM_restraint_force", 100.0 * kilocalories_per_mole / angstroms) 

if not os.path.exists(args.output):
    os.mkdir(args.output)

tstep = args.timestep * femtoseconds
pressure = args.pressure * atmosphere
temperature = args.temp * kelvin
store = args.steps // args.freq
output_log_path = args.output + "/log.txt"
output_pdb_path = args.output + "/output.pdb"
output_energy = args.output + "/energy.csv"

log = open(output_log_path, 'w') # write run info to this log rather than stdout

log.write(f"Reading input pdb {args.input}\n")
pdb = PDBFile(args.input)

forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)

if len(list(modeller.topology.chains())) != 2: 
    raise Exception("Topology does not have exactly two chains.")
    
for i, c in enumerate(modeller.topology.chains()):
    if i == 0: protein = [res for res in c.residues()]
    if i == 1: peptide = [res for res in c.residues()]

log.write("Adding hydrogen atoms.\n")
modeller.addHydrogens(forcefield)

# Optionally explicit water solvent
if args.add_water: 
    log.write("Add water solvent\n")
    modeller.addSolvent(forcefield, model='tip3p', 
                        positiveIon='Na+', negativeIon='Cl-',
                        padding= 1 * nanometer, 
                        neutralize=True)

# create system needs to be after adding solvent and hydrogens
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=app.CutoffNonPeriodic, 
                                 nonbondedCutoff=1*nanometer, 
                                 removeCMMotion=True,
                                 constraints=HBonds) 

if args.suppress_movement != '': 
    # Optionally change masses of these residues to zero to suppress their motion
    # Previously had issue when adding too many residues due to contraints=HBonds,
    # though this was not modifying CA mass exclusively

    hold_residues = []
    with open(args.suppress_movement, 'r') as holdfile: 
        for line in holdfile: 
            splitline = line.strip().split()
            hold_residues.extend([int(index) for index in splitline])
    log.write(f"Holding residue alpha carbons: {hold_residues}\n")
    for resnum, residue in enumerate(protein):
        if (resnum + 1) in hold_residues:
            for atom in residue.atoms():
                if atom.name == "CA": # Only modify alpha carbon mass
                    system.setParticleMass(int(atom.id),0)

t0_protein_center = center_of_mass(modeller.positions, system,
        [a for i in protein for a in i.atoms()])
t0_peptide_center = center_of_mass(modeller.positions, system,
        [a for i in peptide for a in i.atoms()])
log.write(f"Initial protein center of mass: {t0_protein_center}\n")
log.write(f"Initial peptide center of mass: {t0_peptide_center}\n")

# Add custom force pulling on the peptide to the system
peptide_atoms = [atom.index for residue in peptide for atom in residue.atoms()]
pull_direction = t0_peptide_center - t0_protein_center
pull_direction = pull_direction / np.linalg.norm(pull_direction) # make unit vector
pullforce = custom_force(pull_direction, peptide_atoms, args.pull_force)
system.addForce(pullforce)
log.write(f"Pull force {args.pull_force} kcal/mole/Angstrom added in {pull_direction} direction.")

# Enforce a sphere boundary of radius 10nm
sphere_container = container_restraint(range(system.getNumParticles()), 10)
system.addForce(sphere_container)
log.write(f"Spherical container of radius 10 nm added.") 

integrator = LangevinIntegrator(temperature, 1/picosecond, tstep) 

simulation = Simulation(modeller.topology, system, integrator, openmm.Platform.getPlatformByName('CUDA'))
simulation.context.setPositions(modeller.positions)

log.write(f"Temperature (K): {temperature}\n")
log.write(f"Pressure (atm): {pressure}\n")
t0_sim_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
log.write(f"Initial potential energy: {t0_sim_energy}\n")

simulation.minimizeEnergy()
minimized_sim_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
log.write(f"Potential energy after minimization: {minimized_sim_energy}\n")

log.write(f"{args.steps} steps\n")
log.write(f"Time step: {tstep}\n")
log.write(f"energy.csv rows = {store}\n")

simulation.reporters.append(PDBReporter(output_pdb_path, args.steps // 100)) # only store 100 frames in trajectory pdb
simulation.reporters.append(StateDataReporter(output_energy, 
                                              store, 
                                              step=True, 
                                              time=True,
                                              potentialEnergy=True, 
                                              kineticEnergy=True,
                                              totalEnergy=True,
                                              temperature=True,
                                              elapsedTime=True,
                                              progress=True,
                                              remainingTime=True,
                                              totalSteps=args.steps))

simulation.step(args.steps)
