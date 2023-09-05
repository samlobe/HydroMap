from openmm.app import *
from openmm import *
from openmm.unit import *
import sys

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

# Setup the Simulation
simulation = Simulation(top.topology, system, integrator, platform)
simulation.context.setPositions(gro.positions)

# Energy Minimization
print('Performing energy minimization...')
simulation.minimizeEnergy()

# 3. Set reporter frequency
report_frequency_ps = 10  # Every 10 ps
steps_per_report = int(report_frequency_ps / (dt/picoseconds))

# Set total simulation time in ns and calculate the number of steps
total_simulation_time = 0.5  # ns
steps = int(total_simulation_time * 1e3 / (dt/picoseconds))  # Convert total time to ps and divide by timestep

# Production Run
print(f'Simulating for {total_simulation_time} ns...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(DCDReporter('traj.dcd', steps_per_report))
simulation.reporters.append(StateDataReporter('energies.log', steps_per_report, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, separator='\t'))
simulation.reporters.append(CheckpointReporter('checkpoint.chk', 5*steps_per_report))  # checkpoint is 5x less frequent

simulation.step(steps)

# Save the final state
simulation.saveState("endState")

print('Done! Saved trajectory (traj.dcd), state data (energies.log), and checkpoint files if you want to keep simulating later (checkpoint.chk, endState).')
