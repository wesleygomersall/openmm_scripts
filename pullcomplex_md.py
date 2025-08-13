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
parser.add_argument("--suppress-movement", action="store_true", help="Add this option to hold atoms in the protein chain in place throughout the simulation. If specified, program will look for file called 'hold.txt' containing residue numbers to hold stationary. 'hold.txt' should contain integer IDs separated by lines and/or spaces.") 
parser.add_argument("--periodic", action="store_true", help="Add this option to use PME for nonbonded interaction, enforcing a periodic box on the trajectory.") 
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

    pull = CustomCentroidBondForce(1, "force_constant * (x1*fx + y1*fy + z1*fz)")
    pull.addGlobalParameter("force_constant", force_constant)
    pull.addGlobalParameter("fx", fx)
    pull.addGlobalParameter("fy", fy)
    pull.addGlobalParameter("fz", fz)
    pull.addGroup(atoms)
    pull.addBond([0])

    return pull

def restraining_force(direction, atoms, anchor_point, distance, force_constant):
    '''
    Define force which is applied only if there is sufficient space between the two groups.  

    Input(s): 
        direction (list):       (x, y, z) vector for direction of force.
        atoms (str):            List of atoms ids to apply force to.
        anchor_point (str):     (x, y, z) coordinate for the point to restrain around. 
        distance (list):        Distance in nm before restraining force is applied.
        f_constant (float):     Constant force * kcal/mole/A

    Output(s): 
        openmm.CustomCentroidBondForce object
    '''
    force_constant *= kilocalories_per_mole / angstroms
    # distance *= nanometer

    fx = direction[0]
    fy = direction[1]
    fz = direction[2]

    restrain = CustomCentroidBondForce(1, "step(sqrt((x1-t0_x)^2 + (y1-t0_y)^2 + (z1-t0_z)^2) - distance) * -1 * force_constant * (x1*fx + y1*fy + z1*fz)")
    restrain.addGlobalParameter("force_constant", force_constant)
    restrain.addGlobalParameter("fx", fx)
    restrain.addGlobalParameter("fy", fy)
    restrain.addGlobalParameter("fz", fz)
    restrain.addGlobalParameter("t0_x", anchor_point[0] * nanometer)
    restrain.addGlobalParameter("t0_y", anchor_point[1] * nanometer)
    restrain.addGlobalParameter("t0_z", anchor_point[2] * nanometer)
    restrain.addGlobalParameter("distance", distance)
    restrain.addGroup(atoms)
    restrain.addBond([0])

    return restrain

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
    # take care of possibly unbound protein and peptide list of residues
    raise Exception("Topology does not have exactly two chains.")
    
for i, c in enumerate(modeller.topology.chains()):
    if i == 0: protein = [res for res in c.residues()]
    if i == 1: peptide = [res for res in c.residues()]

# For tracking peptide stability, get the second and the second-to-last alpha carbon ids
CA_1_id, CA_2_id = 0, 0
for resnum, res in enumerate(peptide): 
    if resnum == 1:
        for atom in res.atoms(): 
            if atom.name == 'CA': 
                CA_1_id = atom.index
    if resnum == len(peptide) - 2:
        for atom in res.atoms(): 
            if atom.name == 'CA': 
                CA_2_id = atom.index

log.write("Adding hydrogen atoms.\n")
modeller.addHydrogens(forcefield)

# Optionally explicit water solvent
if args.add_water: 
    log.write("Adding solvent: water\n")
    modeller.addSolvent(forcefield, model='tip3p', 
                        positiveIon='Na+', negativeIon='Cl-',
                        padding= 1 * nanometer, 
                        neutralize=True)

# create system needs to be after adding solvent and hydrogens
if args.periodic:
    log.write("nonbondedMethod: PME\n")
    system = forcefield.createSystem(modeller.topology, 
                                    nonbondedMethod=PME, 
                                    nonbondedCutoff=1*nanometer, 
                                    removeCMMotion=False,
                                    constraints=HBonds) # remove this line if error with constraints involving zero-mass atoms
if not args.periodic: 
    log.write("nonbondedMethod: CutoffNonPeriodic\n")
    system = forcefield.createSystem(modeller.topology, 
                                    nonbondedMethod=app.CutoffNonPeriodic,
                                    nonbondedCutoff=1*nanometer, 
                                    removeCMMotion=False,
                                    constraints=HBonds) # remove this line if error with constraints involving zero-mass atoms

if args.suppress_movement: 
    # Optionally change mass of some residues' atoms to 0 to suppress movement
    # Warning: 
    #   This was causing N-terminal nitrogen atom of peptide to also remain locked 
    #   in place. Mass of that atom had not been set to zero.  

    hold_residues = []
    with open('hold.txt', 'r') as holdfile:
        for line in holdfile:
            splitline = line.strip().split()
            hold_residues.extend([int(index) for index in splitline])
    log.write("Setting subset of atoms' mass to zero to suppress motion.\n")
    for resnum, residue in enumerate(protein):
        if (resnum + 1) in hold_residues: 
            for atom in residue.atoms():
                print(atom)
                system.setParticleMass(int(atom.id), 0)
    # Debug suppress movement
    # for r in protein: 
        # for a in r.atoms(): print(a, system.getParticleMass(int(a.id)))
    # for r in peptide: 
        # for a in r.atoms(): print(a, system.getParticleMass(int(a.id)))

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
restraint = restraining_force(pull_direction, peptide_atoms, t0_peptide_center, 8, args.pull_force)
system.addForce(pullforce)
system.addForce(restraint)
log.write(f"Pull force {args.pull_force} kcal/mole/Angstrom added in {pull_direction} direction.")

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

log.write(f"Running simulation for {args.steps} steps.\n")
log.write(f"Time step: {tstep}\n")
log.write(f"Running simulation, storing data every {store} iterations.\n")

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
