from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import os

# Read structure file using sys.argv
if len(sys.argv) < 2:
    print("Usage: simulate_with_openmm.py <structure_file[.gro]>")
    sys.exit(1)

gro_file_name = sys.argv[1]
if not gro_file_name.endswith(".gro"):
    gro_file_name += ".gro"

# Load Gromacs Files
gro = GromacsGroFile(gro_file_name)
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

# Set total simulation time in ns and calculate the number of steps
total_simulation_time = 6  # ns
steps = int(total_simulation_time * 1e3 / (dt/picoseconds))  # Convert total time to ps and divide by timestep

# Setup the Simulation
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)

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
    simulation.reporters.append(DCDReporter('traj.dcd', steps_per_report,append=True))
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

    # Set up reporters to report coordinates, energies, and checkpoints
    simulation.reporters.append(DCDReporter('traj.dcd', steps_per_report))
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

print('Done! Saved trajectory (traj.dcd), state data (energies.log), and checkpoint files if you want to keep simulating later (checkpoint.chk, endState).')

# clean up the duplicate frames from the trajectory and energies log if a checkpoint was used
if remove_duplicates:
    print('\nRemoving duplicate frames in traj.dcd and duplicate entries in energies.log...')
    import subprocess
    subprocess.run(["python", "remove_checkpointed_duplicates.py",gro_file_name])
