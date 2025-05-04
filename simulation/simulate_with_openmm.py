#%%
from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, kelvin, kilojoules_per_mole, atmosphere, angstrom
import os
from time import time
from tqdm import tqdm
import argparse
import sys

parser = argparse.ArgumentParser(description='Run a NPT simulation from a processed protein file.')
parser.add_argument('protein', help='Name of the processed protein structure file (.pdb) for the simulation job, e.g. myProtein_processed.pdb')
parser.add_argument('topol', help='Name of the Gromacs topology file (.top) for the simulation job, e.g. topol.top')
parser.add_argument('-ns','--nanoseconds',default=5,type=float,help='Time in ns you wish to simulate.')
parser.add_argument('-r','--restrain',action='store_true',help='Restrain heavy atoms of protein.')
parser.add_argument('--random_seed',default=42,type=int,help='Random seed for the simulation.')
parser.add_argument('-o','--output',type=str,help='Output trajectory file name (.dcd), will default to {protein_name}_traj.dcd')
parser.add_argument('--noCUDA',action='store_true',help='set to avoid using CUDA.')
args = parser.parse_args()

# example usage: python simulate_with_openmm.py myProtein_processed -ns 5 -r -o traj.dcd

# Read the processed protein structure file (i.e. solvated and neutralized)
protein_file = args.protein
if not protein_file.endswith(".pdb"):
    protein_file += ".pdb"

# extract the protein name from the file name (e.g. myProtein_processed.pdb -> myProtein_processed)
protein_name = os.path.splitext(os.path.basename(protein_file))[0]
if protein_name.endswith("_processed"):
    protein_name = protein_name[:-10]

pdb = PDBFile(protein_file)

# find box vectors in the PDB file
box_vectors = None
with open(protein_file, 'r') as f:
    for line in f:
        if line.startswith('CRYST1'):
            fields = line.split()
            a = float(fields[1]) * 0.1 # convert to nm
            b = float(fields[2]) * 0.1
            c = float(fields[3]) * 0.1
            box_vectors = (Vec3(a, 0, 0), Vec3(0, b, 0), Vec3(0, 0, c)) # only orthogonal box vectors are supported (or else we need to update the water triplet minimum image convention)
            # check if the box is orthogonal and throw error if not; values are usually 90 90 90
            if not (abs(float(fields[4]) - 90) < 1e-3 and abs(float(fields[5]) - 90) < 1e-3 and abs(float(fields[6]) - 90) < 1e-3):
                raise ValueError("Box vectors are not orthogonal. Please provide a box with orthogonal vectors.\n(only orthogonal box vectors are supported for water triplet analysis, unless we update the minimum image convention in waterlib.c)")
            break
if box_vectors is None:
    raise ValueError("Can't find box size in PBD file (looking for a line starting with CRYST1)")


# load gromacs topology file
top = GromacsTopFile(args.topol, periodicBoxVectors=box_vectors)


pressure = 1*atmosphere  # Store pressure
print("Pressure:", pressure)

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometer
constraints = HBonds
system = top.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, constraints=constraints)

# Pressure & Barostat
temperature = 300*kelvin
print("Temperature:", temperature)
barostatInterval = 25
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

# Integration Options
dt = 0.004*picosecond  # 4 fs timestep
friction = 2/picosecond
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setRandomNumberSeed(args.random_seed)

# Setup Platform for GPU

# Check if CUDA is available
if args.noCUDA:
    print("Ignoring CUDA...")
    try: 
        platform = Platform.getPlatformByName('OpenCL')
        print("Falling back to OpenCL platform.")
    except Exception as e:
        print("Falling back to CPU (will be extremely slow).")
        platform = Platform.getPlatformByName('CPU')
else:
    try:
        platform = Platform.getPlatformByName('CUDA')
        print("Using CUDA platform.")
    except Exception as e:
        print("Error: CUDA platform not avaialble.")
        print("Hint: Enable CUDA or use the --noCUDA flag to fall back to OpenCL (slower) or CPU (extremely slow).")
        sys.exit(1)

# Set reporter frequency
report_frequency_ps = 1  # Every 1 ps
steps_per_report = int(report_frequency_ps / (dt/picosecond))
steps_per_checkpoint = int(steps_per_report)*100

# Equilibration phase (before production run)
equilibration_time_ps = 100  # 100 ps equilibration time
equilibration_steps = int(equilibration_time_ps / (dt/picosecond))

# Set total production simulation time in ns and calculate the number of steps
total_simulation_time = args.nanoseconds  # ns
steps = int(total_simulation_time * 1e3 / (dt/picosecond))  # Convert total time to ps and divide by timestep

# Add restraints to heavy atoms in the protein if -r flag is set
if args.restrain:
    # Define a force for restraining atoms
    force = CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", 1000.0 * kilojoules_per_mole/nanometer**2)  # force constant
    force.addPerParticleParameter("x0")  # x coordinate of the restrained atom
    force.addPerParticleParameter("y0")  # y coordinate of the restrained atom
    force.addPerParticleParameter("z0")  # z coordinate of the restrained atom

    # Loop through all residues in the topology
    for residue in top.topology.residues():
        # print(residue.name)
        # Check if the residue is not water or typical ion names (you can expand this list if needed)
        # alternatively uncomment this to instead check if residue.name is in a list of amino acid names (may not work for odd residue names)
        # if residue.name in ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
        #                     'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
        #                     'ACE','NME','HID','HIE','HIP','LYN','ORN','CYX','CME','CYM',
        #                     'ASH','GLH''HYP','NH2','NHE','NH3']:
        if residue.name not in ["SOL","HOH","HO4","Na+","NA","Cl-","CL","K+","K","Mg2+","MG","Ca2+","CA","ZN","URE"]:
            # Loop through all atoms in the residue
            for atom in residue.atoms():
                # Check if the atom is not a hydrogen
                if atom.element.symbol != 'H':
                    force.addParticle(atom.index, pdb.positions[atom.index].value_in_unit(nanometer))
                    print(f"Adding restraint to atom: {atom.name} in residue: {residue.name}")

    # Add the restraining force to the system
    system.addForce(force)

# Setup the Simulation
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
if args.output is None:
    traj_name = f'{protein_name}_traj.dcd'
else:
    traj_name = args.output

# Load from the checkpoint if it exists
checkpoint_file = f'{protein_name}_checkpoint.chk'
energies_file = f'{protein_name}_energies.log'
endState_file = f'{protein_name}_endState'

if os.path.exists(checkpoint_file):
    tik = time()
    print("Found checkpoint file. Resuming simulation from the checkpoint.")
    # Load from the checkpoint
    with open(checkpoint_file, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

    # subtract equilibration (equilibration_time_ps (=100 ps) and equilibration_steps)
    simulation.context.setTime(simulation.context.getTime() - equilibration_time_ps* picosecond)
    simulation.context.setStepCount(simulation.context.getStepCount() - equilibration_steps)

    # Adjust the number of steps to simulate based on the total desired steps and the current step count
    steps_remaining = steps - simulation.currentStep
    if steps_remaining < 0:
        steps_remaining = 0

    # Add the reporters
    simulation.reporters.append(DCDReporter(traj_name, steps_per_report, append=True))
    simulation.reporters.append(StateDataReporter(energies_file, steps_per_report, step=True, time=True,
                                                  potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                                  temperature=True, volume=True, separator='\t', append=True))
    simulation.reporters.append(CheckpointReporter(checkpoint_file, steps_per_checkpoint))
    remove_duplicates = True

    # Continue the simulation with progress bar
    print(f"Continuing simulation for {steps_remaining * dt.value_in_unit(picosecond) * 1e-3:.1f} / {total_simulation_time:.1f} ns...")
    print("Showing progress bar in 0.1 ns increments:")
    # Set up simulation in 0.1 ns increments
    steps_per_point1ns = int(0.1 * 1e3 / dt.value_in_unit(picosecond))  # 0.1 ns = 100 ps
    num_points = steps_remaining // steps_per_point1ns
    leftover_steps = steps_remaining % steps_per_point1ns

    for _ in tqdm(range(num_points), desc="Simulating", unit="0.1 ns"):
        simulation.step(steps_per_point1ns)
    if leftover_steps > 0:
        simulation.step(leftover_steps)

else:
    # If checkpoint doesn't exist, perform energy minimization
    print('Performing energy minimization...')
    simulation.minimizeEnergy()
    print(f'Performing {total_simulation_time:.1f} ns simulation...\nShowing progress bar in 0.1 ns increments:')
    tik = time()
    simulation.step(equilibration_steps)
    # reset time & step count so production starts at t=0
    simulation.context.setTime(0.0)
    simulation.context.setStepCount(0)

    # Set up reporters
    simulation.reporters.append(DCDReporter(traj_name, steps_per_report))
    simulation.reporters.append(StateDataReporter(energies_file, steps_per_report, step=True, time=True,
                                                  potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                                  temperature=True, volume=True, separator='\t'))
    simulation.reporters.append(CheckpointReporter(checkpoint_file, steps_per_checkpoint))
    simulation.context.setVelocitiesToTemperature(temperature)
    remove_duplicates = False

    # Production Run with progress bar
    steps_per_point1ns = int(0.1 * 1e3 / dt.value_in_unit(picosecond))  # 0.1 ns = 100 ps
    num_points = steps // steps_per_point1ns
    leftover_steps = steps % steps_per_point1ns

    for _ in tqdm(range(num_points), desc="Simulating", unit="0.1 ns"):
        simulation.step(steps_per_point1ns)
    if leftover_steps > 0:
        simulation.step(leftover_steps)

# Save the final state
simulation.saveState(endState_file)

print(f'Done! Saved trajectory ({traj_name}), state data ({energies_file}), and checkpoint files if you want to keep simulating later ({checkpoint_file}, {endState_file}).')

tok = time()
print(f'Total wall clock time: {(tok - tik)/60:.1f} minutes')
print(f"We recommend visualizing your trajectory (e.g. with Pymol or ChimeraX or VMD),\nand to evaluating the energies/temperature/box volumes in the {energies_file} file, to check if the simulation is healthy.")

# clean up the duplicate frames from the trajectory and energies log if a checkpoint was used
if remove_duplicates:
    print(f'\nRemoving duplicate frames in {args.output} and duplicate entries in {energies_file}...')
    import subprocess
    subprocess.run(["python", "remove_checkpointed_duplicates.py",protein_file,energies_file, traj_name])

# %%
