#%%
from openmm.app import *
from openmm import *
from openmm.unit import nanometer, picosecond, kelvin, kilojoules_per_mole, atmosphere, angstrom
from pathlib import Path
import os
from time import time
from tqdm import tqdm
import argparse
import sys

parser = argparse.ArgumentParser(description='Run a NPT simulation from a processed protein file.')
parser.add_argument('protein', help='Name of the processed protein structure file (.pdb) for the simulation job, e.g. myProtein_processed.pdb')
parser.add_argument(
    'topol',
    help='Serialized OpenMM System XML for the simulation job.'
)
parser.add_argument('-ns','--nanoseconds',default=1,type=float,help='Time in ns you wish to simulate.')
parser.add_argument(
    '-r','--restrain',
    nargs='?',            # optional value
    const='protein',      # if provided with no value, use "protein"
    default=None,         # if not provided, no restraints
    help='Selection for atoms to restrain (using MDAnalysis selection language). If the restraint flag is given with no selection, it will default to restraining all protein heavy atoms.'
         'Examples: --restrain "resid 295-311 and name CA", --restrain "backbone"'
)
parser.add_argument(
    '--restraint_k', type=float, default=1000.0,
    help='Restraint force constant (kJ/mol/nm^2). Default: 1000'
)
parser.add_argument('--random_seed',default=42,type=int,help='Random seed for the simulation.')
parser.add_argument('--velocity_seed',default=None,type=int,help='Random seed used when initializing velocities. Defaults to --random_seed.')
parser.add_argument('--barostat_seed',default=None,type=int,help='Random seed used by MonteCarloBarostat. Defaults to --random_seed.')
parser.add_argument('-o','--output',type=str,help='Output trajectory file name (.dcd), will default to {protein_name}_traj.dcd')
parser.add_argument('--checkpoint_policy',choices=['error','resume','overwrite'],default='error',help='How to handle existing checkpoint/output files: error (default), resume, or overwrite.')
parser.add_argument('--deterministic',action='store_true',help='Enable strict deterministic platform settings (when available).')
parser.add_argument('--cuda_precision',choices=['single','mixed','double'],default='mixed',help='CUDA precision mode (default: mixed).')
parser.add_argument('--equilibration_ps',type=float,default=100.0,help='Equilibration time before production, in ps (default: 100).')
parser.add_argument('--timestep_ps', type=float, default=0.003, help='Integrator timestep in ps (default: 0.003, i.e. 3 fs).')
parser.add_argument('--report_interval_ps', type=float, default=0.5, help='Trajectory/report interval in ps (default: 0.5).')
parser.add_argument('--noCUDA',action='store_true',help='set to avoid using CUDA.')
parser.add_argument('--debug',action='store_true',help='print additional debug information.')
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
    raise ValueError("Could not find box vectors in the PDB file (expected a CRYST1 record).")

if not args.topol.endswith(".xml"):
    raise ValueError("Prepared system input must be a serialized OpenMM System XML file.")
with open(args.topol, "r", encoding="utf-8") as handle:
    system = XmlSerializer.deserialize(handle.read())
system.setDefaultPeriodicBoxVectors(*box_vectors)
simulation_topology = pdb.topology


pressure = 1*atmosphere  # Store pressure
print("Pressure:", pressure)

# Pressure & Barostat
temperature = 300*kelvin
print("Temperature:", temperature)
barostatInterval = 25
integrator_seed = args.random_seed
velocity_seed = args.velocity_seed if args.velocity_seed is not None else integrator_seed
barostat_seed = args.barostat_seed if args.barostat_seed is not None else integrator_seed
print(f"Seeds: integrator={integrator_seed}, velocities={velocity_seed}, barostat={barostat_seed}")
barostat = MonteCarloBarostat(pressure, temperature, barostatInterval)
barostat.setRandomNumberSeed(barostat_seed)
system.addForce(barostat)

# Integration Options
if args.timestep_ps <= 0:
    raise ValueError("--timestep_ps must be > 0.")
if args.report_interval_ps <= 0:
    raise ValueError("--report_interval_ps must be > 0.")
if args.report_interval_ps < args.timestep_ps:
    raise ValueError("--report_interval_ps must be >= --timestep_ps.")

dt = args.timestep_ps * picosecond
friction = 2/picosecond
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setRandomNumberSeed(integrator_seed)

# Setup Platform for GPU

# Check if CUDA is available
if args.noCUDA:
    print("CUDA disabled; using a non-CUDA platform.")
    if args.deterministic:
        print("Using CPU platform for deterministic non-CUDA execution.")
        platform = Platform.getPlatformByName('CPU')
    else:
        try: 
            platform = Platform.getPlatformByName('OpenCL')
            print("Using OpenCL platform.")
        except Exception as e:
            print("Using CPU platform. Performance will be slow.")
            platform = Platform.getPlatformByName('CPU')
else:
    try:
        platform = Platform.getPlatformByName('CUDA')
        print("Using CUDA platform.")
    except Exception as e:
        print("Error: CUDA platform not available.")
        print("Hint: Enable CUDA or use --noCUDA to run with OpenCL or CPU.")
        sys.exit(1)

platform_properties = {}
try:
    supported_platform_properties = set(platform.getPropertyNames())
except Exception:
    supported_platform_properties = set()
if platform.getName() == "CUDA":
    if "Precision" in supported_platform_properties:
        platform_properties["Precision"] = args.cuda_precision
    elif "CudaPrecision" in supported_platform_properties:
        platform_properties["CudaPrecision"] = args.cuda_precision
    if args.deterministic:
        if "DeterministicForces" not in supported_platform_properties:
            raise RuntimeError("CUDA platform does not expose DeterministicForces; cannot run with --deterministic.")
        platform_properties["DeterministicForces"] = "true"
elif args.deterministic and platform.getName() == "OpenCL":
    # OpenCL deterministic behavior is backend-specific and less standardized.
    print("Warning: --deterministic requested with OpenCL; strict deterministic force setting is CUDA-specific.")
elif args.deterministic and platform.getName() == "CPU":
    print("Note: --deterministic requested on CPU. DeterministicForces setting is not required on this platform.")

# Set reporter frequency
report_frequency_ps = args.report_interval_ps
steps_per_report = max(1, int(round(report_frequency_ps / (dt/picosecond))))
steps_per_checkpoint = int(steps_per_report)*100
actual_report_ps = steps_per_report * (dt / picosecond)
print(f"Trajectory/report interval target: {report_frequency_ps} ps (actual {actual_report_ps:.6f} ps)")

# Equilibration phase (before production run)
equilibration_time_ps = args.equilibration_ps
equilibration_steps = int(equilibration_time_ps / (dt/picosecond))

# Set total production simulation time in ns and calculate the number of steps
total_simulation_time = args.nanoseconds  # ns
steps = int(total_simulation_time * 1e3 / (dt/picosecond))  # Convert total time to ps and divide by timestep

# ---- Restraints (via MDAnalysis selection) ----
if args.restrain is not None:
    try:
        import MDAnalysis as mda
    except ImportError:
        raise ImportError("MDAnalysis is required for --restrain. Install with `conda install MDAnalysis`.")

    selection_str = args.restrain if isinstance(args.restrain, str) else 'protein'
    # Build an MDAnalysis universe from the same PDB you loaded above
    u = mda.Universe(protein_file)

    # Sanity checks to ensure atom order/counts align (important because we index into pdb.positions/system)
    n_pdb = len(pdb.positions)
    n_u = u.atoms.n_atoms
    n_sys = system.getNumParticles()
    if n_u != n_pdb or n_sys != n_pdb:
        raise ValueError(
            f"Atom count mismatch: MDAnalysis={n_u}, PDBFile={n_pdb}, System={n_sys}. "
            "Ensure your prepared system input and .pdb describe the same atom ordering."
        )

    sel = u.select_atoms(selection_str)
    if sel.n_atoms == 0:
        raise ValueError(f"--restrain selection matched 0 atoms: '{selection_str}'")

    # Create the harmonic positional restraint about initial coordinates
    force = CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", args.restraint_k * kilojoules_per_mole/nanometer**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    # Add all selected atoms
    for idx in sel.atoms.indices:  # 0-based; matches PDB order
        xyz = pdb.positions[int(idx)].value_in_unit(nanometer)
        force.addParticle(int(idx), xyz)

    system.addForce(force)
    print(f"Added harmonic restraints to {sel.n_atoms} atoms ")
    if args.debug:
        print(f"(selection: '{selection_str}', k={args.restraint_k} kJ/mol/nm^2).")
# -----------------------------------------------


# Setup the Simulation
if platform_properties:
    print(f"Platform properties: {platform_properties}")
    simulation = Simulation(simulation_topology, system, integrator, platform, platform_properties)
else:
    simulation = Simulation(simulation_topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
if args.output is None:
    traj_name = f'{protein_name}_traj.dcd'
else:
    traj_name = args.output

# Load from the checkpoint if it exists
checkpoint_file = f'{protein_name}_checkpoint.chk'
energies_file = f'{protein_name}_energies.log'
endState_file = f'{protein_name}_endState'
checkpoint_exists = os.path.exists(checkpoint_file)

if checkpoint_exists and args.checkpoint_policy == "error":
    raise FileExistsError(
        f"Found existing checkpoint '{checkpoint_file}' with checkpoint policy 'error'. "
        "Use --checkpoint_policy resume to continue or --checkpoint_policy overwrite to restart."
    )

if checkpoint_exists and args.checkpoint_policy == "overwrite":
    print(f"Found checkpoint and outputs for {protein_name}; removing them because checkpoint policy is 'overwrite'.")
    for path in [checkpoint_file, energies_file, endState_file, traj_name]:
        if os.path.exists(path):
            os.remove(path)
    checkpoint_exists = False

if checkpoint_exists and args.checkpoint_policy == "resume":
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
    simulation.context.setVelocitiesToTemperature(temperature, velocity_seed)
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

print(
    f"Saved trajectory ({traj_name}), state data ({energies_file}), "
    f"and checkpoint/state files ({checkpoint_file}, {endState_file})."
)

tok = time()
print(f'Total wall clock time: {(tok - tik)/60:.1f} minutes')
print(
    "Inspect the trajectory and review the energies, temperature, and box volume "
    f"in {energies_file} to confirm the simulation behaved as expected."
)

# clean up the duplicate frames from the trajectory and energies log if a checkpoint was used
if remove_duplicates:
    print(f'\nRemoving duplicate frames in {traj_name} and duplicate entries in {energies_file}...')
    import subprocess
    subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "remove_checkpointed_duplicates.py"), protein_file, energies_file, traj_name])

# %%
