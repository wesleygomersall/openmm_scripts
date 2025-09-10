#!/usr/bin/env python3

# This script was adapted from pullcomplex_md.py 2025-09-09
# Changed from SMD: Remove pulling force and add option for variable integrator

from openmm import LangevinIntegrator, VariableLangevinIntegrator, CustomExternalForce, Platform
from openmm.app import ForceField, PDBFile, Modeller, Simulation, StateDataReporter, PDBReporter, CutoffNonPeriodic, HBonds
# from openmm.unit import nanometer, femtoseconds, kelvin, atmosphere, kilocalories_per_mole, angstroms
import openmm.unit as unit

# from sys import stdout
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser(description="MD trajectory simulating a protein-peptide complex over time.")
parser.add_argument("--input", type=str, help="Path to input pdb containing both protein and peptide. If errors exist, try using PDBfixer first.")
parser.add_argument("--output", type=str, help="Path to output directory. Will be created if not already existing.")
parser.add_argument("--timestep", type=float, default=2.0, help="Time between steps in femtoseconds, default is 2.")
parser.add_argument("--steps", type=int, default=10_000, help="Total steps in simulation, default is ten thousand.")
parser.add_argument("--freq", type=int, default=100, help="Frequency of reporter, will record data `freq` amount of times (every `step // freq` steps).")
parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (atmospheres), default is 1 atm.")
parser.add_argument("--temp", type=float, default=300.0, help="Temperature (Kelvin), default is 300K.")
parser.add_argument("--add-water", action="store_true", help="Add this option to populate modeller with a surrounding box of water molecules.")
parser.add_argument("--variable-int", action="store_true", help="Add this option to use variable integrator for constructing the system.")
args = parser.parse_args()


def container_restraint(atoms, radius: float):
    '''
    Define force to contain peptide chain to a spherical container. 

    Input(s): 
        atoms (list):           List of atoms ids to apply force to.
        radius (float):         Radius of container in nm. 

    Output(s): 
        openmm.CustomExternalForce object
    '''
    radius *= unit.nanometer
    force_expression = "container_force*max(0, r - radius)^2; r = sqrt(x^2 + y^2 + z^2) "
    container = CustomExternalForce(force_expression)
    container.addGlobalParameter("radius", radius) 
    container.addGlobalParameter("container_force", 100.0 * unit.kilocalories_per_mole / unit.angstroms) 
    for p in atoms: 
        container.addParticle(p, [])
    return container


if not os.path.exists(args.output):
    os.mkdir(args.output)

tstep = args.timestep * unit.femtoseconds
total_time = tstep * args.steps
pressure = args.pressure * unit.atmosphere
temperature = args.temp * unit.kelvin
store = int(args.steps) // int(args.freq)
pdbstore = int(args.steps) // 100 # only store 100 frames in trajectory pdb
if pdbstore == 0: pdbstore += 1 # if frames == 100 this will happen  
if store == 0: store += 1 # if frames == 100 this will happen  
output_log_path = args.output + "/log.txt"
output_pdb_path = args.output + "/output.pdb"
output_energy = args.output + "/energy.csv"

## TESTING
# print(int(args.steps))


gpu = False
# if gpu: plat = 'CUDA'
# else: plat = 'CPU'

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

if args.add_water: 
    log.write("Add water solvent\n")
    modeller.addSolvent(forcefield, model='tip3p', 
                        positiveIon='Na+', negativeIon='Cl-',
                        padding= 1 * unit.nanometer, 
                        neutralize=True)

# create system needs to be after adding solvent and hydrogens
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=CutoffNonPeriodic, 
                                 nonbondedCutoff=1*unit.nanometer, 
                                 removeCMMotion=False,
                                 constraints=HBonds) 


# Enforce a sphere boundary of radius 20nm
sphere_container = container_restraint(range(system.getNumParticles()), 20)
system.addForce(sphere_container)
log.write(f"Spherical container of radius 20 nm added.") 

if args.variable_int:
    error_tolerance = 0.001
    integrator = VariableLangevinIntegrator(temperature, 1/unit.picosecond, error_tolerance) 
    log.write(f"VariableLangevinIntegrator error tolerance: {error_tolerance}.") 
else: 
    integrator = LangevinIntegrator(temperature, 1/unit.picosecond, tstep) 

simulation = Simulation(modeller.topology, system, integrator, Platform.getPlatformByName('CUDA' if gpu else 'CPU'))
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

simulation.reporters.append(PDBReporter(output_pdb_path, pdbstore)) 
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


if args.variable_int:
    integrator.stepTo(total_time)
else: 
    simulation.step(int(args.steps))
