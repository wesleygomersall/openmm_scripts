#!/usr/bin/env python3

from openmm.app import *
from openmm import *
from openmm.unit import *

from sys import stdout
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="MD trajectory simulating constant force pulling apart a protein-ligand complex.")
parser.add_argument("--input", type=str, help="Path to input pdb containing both protein and ligand. If errors exist, try using PDBfixer first.")
parser.add_argument("--output", type=str, help="Path to output directory. Will be created if not already existing.")
parser.add_argument("--timestep", type=float, default=2.0, help="Time between steps in femtoseconds, default is 2.")
parser.add_argument("--steps", type=int, default=10_000, help="Total steps in simulation, default is ten thousand.")
parser.add_argument("--freq", type=int, default=100, help="Frequency of reporter, will record data `freq` amount of times (every `step // freq` steps).")
parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (atmospheres), default is 1 atm.")
parser.add_argument("--temp", type=float, default=300.0, help="Temperature (Kelvin), default is 300K.")
parser.add_argument("--pull-force", type=float, default=25.0, help="Pull force constant which will be applied to the ligand.")
parser.add_argument("--add-water", action="store_true", help="Add this option to populate modeller with a surrounding box of water molecules.")
parser.add_argument("--suppress-specific", type=str, default='', help="Movement of some protein alpha carbons will be suppressed for the duration of the simulation. If this not specified, will default to all alpha carbons 1.5 nm from ligand. If specified, program will look for file specified which contains residue numbers to hold stationary. This file should contain integer IDs separated by lines and/or spaces. If this file is empty, all atoms will be allowed to move for the duration of the simulation.") 
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

if not os.path.exists(args.output):
    os.mkdir(args.output)

tstep = args.timestep * femtoseconds
pressure = args.pressure * atmosphere
temperature = args.temp * kelvin
store = args.steps // args.freq
output_log_path = args.output + "/log.txt"
output_pdb_path = args.output + "/output.pdb"
output_dcd_path = args.output + "/output.dcd"
output_energy = args.output + "/energy.csv"

log = open(output_log_path, 'w') # write run info to this log rather than stdout

log.write(f"Reading input pdb {args.input}\n")
pdb = PDBFile(args.input)

forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)

if len(list(modeller.topology.chains())) != 2: 
    raise Exception("Topology does not have exactly two chains.")
    
all_residues = []
for i, c in enumerate(modeller.topology.chains()):
    all_residues.append([res for res in c.residues()])
if len(all_residues[0]) < len(all_residues[1]): 
    protein = all_residues[1]
    ligand = all_residues[0]
elif len(all_residues[0]) > len(all_residues[1]):
    protein = all_residues[0]
    ligand = all_residues[1]

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
                                 nonbondedMethod=app.PME, # barostat requires
                                 nonbondedCutoff=1*nanometer, 
                                 removeCMMotion=False) 

t0_protein_center = center_of_mass(modeller.positions, system,
        [a for i in protein for a in i.atoms()])
t0_ligand_center = center_of_mass(modeller.positions, system,
        [a for i in ligand for a in i.atoms()])
log.write(f"Initial protein center of mass: {t0_protein_center}\n")
log.write(f"Initial ligand center of mass: {t0_ligand_center}\n")

# Add movement suppression and pull force after adjust period
if args.suppress_specific == '': 
    '''
    Cannot use constraints (contraints=HBonds in
    forcefield.createSystem) involving massless particles - Hydrogens bonded
    to alpha carbons) 
    '''
    # beyond this distance from ligand, alpha carbons will be massless
    suppress_distance = 1.5 * nanometer 

    hold_residues = []

    positions = modeller.getPositions()

    for resnum1, residue1 in enumerate(protein):
        for atom1 in residue1.atoms():
            if atom1.name == "CA": 

                distances = [] # alpha carbon distance from all atoms in ligand
                for residue2 in ligand: 
                    for atom2 in residue2.atoms(): 
                        distances.append(norm(positions[int(atom1.id)] - positions[int(atom2.id)]))

                if min(distances) >= suppress_distance: 
                    system.setParticleMass(int(atom1.id),0)
                    hold_residues.append(resnum1)

    log.write(f"Holding residue alpha carbons: {hold_residues}\n")
if args.suppress_specific != '': 
    hold_residues = []
    with open(args.suppress_movement, 'r') as holdfile: 
        for line in holdfile: 
            splitline = line.strip().split()
            hold_residues.extend([int(item) for item in splitline])

    if hold_residues != []: 
        for resnum, residue in enumerate(protein):
            if (resnum + 1) in hold_residues:
                for atom in residue.atoms():
                    if atom.name == "CA": # Only modify alpha carbon mass
                        system.setParticleMass(int(atom.id),0)

        log.write(f"Holding residue alpha carbons: {hold_residues}\n")
    else: 
        log.write(f"Holding residue alpha carbons: None!\n")

pullstart_protein_center = center_of_mass(modeller.positions, system,
        [a for i in protein for a in i.atoms()])
pullstart_ligand_center = center_of_mass(modeller.positions, system,
        [a for i in ligand for a in i.atoms()])
log.write(f"Protein center of mass at pull start: {pullstart_protein_center}\n")
log.write(f"Ligand center of mass at pull start: {pullstart_ligand_center}\n")

# Add custom force pulling on the ligand atoms to the system
pull_direction = pullstart_ligand_center - pullstart_protein_center
pull_direction = pull_direction / np.linalg.norm(pull_direction) # make unit vector
log.write(f"Initializing external force. Start value 0 kcal/mol/A.\n")

restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
system.addForce(restraint)
restraint.addGlobalParameter('k', 0 * kilocalories_per_mole / angstroms)
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')

# Apply pull to only to the first alpha carbon of ligand 
for atom in ligand[0].atoms():
    if atom.name == 'CA': # modify for RNA ligand 
        restraint.addParticle(atom.index, pdb.positions[atom.index] + 3 * pull_direction * nanometer)

# Must create system after adding all of the forces.
integrator = LangevinIntegrator(temperature, 1/picosecond, tstep) 

simulation = Simulation(modeller.topology, system, integrator, openmm.Platform.getPlatformByName('CUDA'))
simulation.context.setPositions(modeller.positions)

log.write(f"Temperature (K): {temperature}\n")
log.write(f"Pressure (atm): {pressure}\n")
t0_sim_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
log.write(f"Initial potential energy: {t0_sim_energy}\n")

simulation.minimizeEnergy()

# NVT (2 fs * 10000 = 20 picoseconds)
simulation.step(10000)
log.write(f"NVT equillibration for 10000 time steps of {args.timestep} fs complete.\n")

# NPT (2 fs * 10000 = 20 picoseconds)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin, 25))
simulation.context.reinitialize(preserveState=True)
simulation.step(10000)
log.write(f"NPT equillibration for 10000 time steps of {args.timestep} fs complete.\n")

minimized_sim_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
log.write(f"Potential energy after minimization and equillibration: {minimized_sim_energy}\n")

with open(output_pdb_path, 'w') as f:
    PDBFile.writeFile(simulation.topology, 
                      simulation.context.getState(getPositions=True).getPositions(),
                      f)

simulation.reporters.append(DCDReporter(output_dcd_path, args.steps // 1000)) # store 1000 frames in trajectory pdb
# simulation.reporters.append(PDBReporter(output_pdb_path, args.steps // 1000)) # Used to use pdb reporter

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
                                              remainingTime=True, totalSteps=args.steps)) # Advance simulation 500k steps (to 0.5 ns for default 2fs step size) 2 fs * 500000 steps = 500 picoseconds simulation.step(500000)



# Advance simulation 500k steps (to 0.5 ns for default 2fs step size)
# 2 fs * 500000 steps = 500 picoseconds 
simulation.step(500000)

# Add nonzero pull force with the specified strength
simulation.context.setParameter('k', args.pull_force * kilocalories_per_mole / angstroms)
log.write(f"Pull force {args.pull_force} kcal/mole/Angstrom added in {pull_direction} direction.\n")

log.write(f"Will run to total steps, time/step of:\n")
log.write(f"{args.steps} steps\n")
log.write(f"Time step: {tstep}\n")
log.write(f"energy.csv rows = {store}\n")

# Advance simulation (args.steps minus X steps / calc steps from time from above)
if args.steps < 500000: args.steps = 500000 + 20000 # minimum 20k time frames with force applied.
simulation.step(args.steps - 500000)
