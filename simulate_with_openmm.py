#%%
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Run a NPT simulation from a processed protein file.')
parser.add_argument('protein', help='Name of the processed protein structure file (.gro) for the simulation job, e.g. myProtein_processed[.gro]')
parser.add_argument('-ns','--nanoseconds',default=5,type=float,help='Time in ns you wish to simulate.')
parser.add_argument('-r','--restrain',action='store_true',help='Restrain heavy atoms of protein.')
parser.add_argument('-o','--output',default='traj.dcd',type=str,help='Output trajectory file name (.dcd)')
args = parser.parse_args()

# Read the processed protein structure file (i.e. solvated and neutralized)
protein_file = args.protein
if not protein_file.endswith(".gro"):
    protein_file += ".gro"

# Load Gromacs Files
gro = GromacsGroFile(protein_file)
top = GromacsTopFile('topol.top', periodicBoxVectors=gro.getPeriodicBoxVectors())

# Store box size and pressure
box_vectors = gro.getPeriodicBoxVectors()
print("Box Vectors:", box_vectors)

pressure = 1*atmospheres  # Store pressure
print("Pressure:", pressure)

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.0*nanometers
constraints = HBonds
system = top.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, constraints=constraints)

# Pressure & Barostat
temperature = 300*kelvin
print("Temperature:", temperature)
barostatInterval = 25
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

# Integration Options
dt = 0.004*picoseconds  # 4 fs timestep
friction = 2/picosecond
integrator = LangevinMiddleIntegrator(temperature, friction, dt)

# Setup Platform for GPU
platform = Platform.getPlatformByName('CUDA')

# Set reporter frequency
report_frequency_ps = 1  # Every 1 ps
steps_per_report = int(report_frequency_ps / (dt/picoseconds))
steps_per_checkpoint = int(steps_per_report)*100

# Equilibration phase (before production run)
equilibration_time_ps = 100  # 100 ps equilibration time
equilibration_steps = int(equilibration_time_ps / (dt/picoseconds))

# Set total production simulation time in ns and calculate the number of steps
total_simulation_time = args.nanoseconds  # ns
steps = int(total_simulation_time * 1e3 / (dt/picoseconds))  # Convert total time to ps and divide by timestep

# Add restraints to heavy atoms in the protein if -r flag is set
# HAVENT TESTED THIS YET
if args.restrain:
    # Define a force for restraining atoms
    force = CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", 1000.0 * kilojoules_per_mole/nanometers**2)  # force constant
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
                    force.addParticle(atom.index, gro.positions[atom.index].value_in_unit(nanometers))
                    print(f"Adding restraint to atom: {atom.name} in residue: {residue.name}")

    # Add the restraining force to the system
    system.addForce(force)

# Setup the Simulation
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)
traj_name = args.output

# Load from the checkpoint if it exists
checkpoint_file = 'checkpoint.chk'
if os.path.exists(checkpoint_file):
    print("Found checkpoint file. Resuming simulation from the checkpoint.")
    # Load from the checkpoint
    with open(checkpoint_file, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())
    # Adjust the number of steps to simulate based on the total desired steps and the current step count
    steps_remaining = steps - simulation.currentStep
    if steps_remaining < 0:
        steps_remaining = 0

    # Add the reporters
    simulation.reporters.append(DCDReporter(traj_name, steps_per_report,append=True))
    simulation.reporters.append(StateDataReporter('energies.log', steps_per_report, step=True,time=True,
                                                  potentialEnergy=True, kineticEnergy=True,totalEnergy=True,
                                                  temperature=True, volume=True, separator='\t',append=True))
    simulation.reporters.append(CheckpointReporter('checkpoint.chk', steps_per_checkpoint))
    # NOTE: the energies.log and traj.dcd files are 99% likely to have duplicates.
    # We can remove these duplicates with `python remove_checkpointed_duplicates.py`
    remove_duplicates = True
    
    # Continue the simulation
    simulation.step(steps_remaining)

else:
    # If checkpoint doesn't exist, perform energy minimization
    print('Performing energy minimization...')
    simulation.minimizeEnergy()
    print(f'Simulating for {total_simulation_time} ns...')
    # Equilibrate briefly before saving properties
    simulation.step(equilibration_steps)

    # Set up reporters to report coordinates, energies, and checkpoints
    simulation.reporters.append(DCDReporter(traj_name, steps_per_report))
    simulation.reporters.append(StateDataReporter('energies.log', steps_per_report, step=True, time=True, 
                                                potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
                                                temperature=True, volume=True, separator='\t'))
    simulation.reporters.append(CheckpointReporter('checkpoint.chk', steps_per_checkpoint))
    simulation.context.setVelocitiesToTemperature(temperature)
    remove_duplicates = False # no need to clean up any repeated frames from checkpointing

    # Production Run
    simulation.step(steps)

# Save the final state
simulation.saveState("endState")

print(f'Done! Saved trajectory ({args.output}), state data (energies.log), and checkpoint files if you want to keep simulating later (checkpoint.chk, endState).')

# clean up the duplicate frames from the trajectory and energies log if a checkpoint was used
if remove_duplicates:
    print(f'\nRemoving duplicate frames in {args.output} and duplicate entries in energies.log...')
    import subprocess
    subprocess.run(["python", "remove_checkpointed_duplicates.py",protein_file])

# %%
