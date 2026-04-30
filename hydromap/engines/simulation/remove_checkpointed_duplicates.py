import MDAnalysis as mda
import numpy as np
import os
import sys

structure_file = sys.argv[1]
energies_file = sys.argv[2]
traj_file = sys.argv[3]

# Step 1: Process the energies.log
times = [] # in ps

with open(energies_file, "r") as file:
    lines = file.readlines()
    for line in lines:
        if line.startswith("#"):  # Assuming lines starting with # are comments or headers
            continue
        # Assuming time step is the first column in the log file, adjust if not
        time = round(float(line.split()[1]))
        times.append(time)

# Backup the original file
backup_energies_file = energies_file[:-4] + "_with_duplicates.log"
os.rename(energies_file, backup_energies_file)

# Step 2: Remove duplicate time steps and overwrite the energies.log
unique_times, indices_to_keep = np.unique(times, return_index=True)
with open(energies_file, "w") as file:
    file.write(lines[0])  # Write the header
    for i in indices_to_keep:
        file.write(lines[i+1])

# Step 3 & 4: Process the traj.dcd with MDAnalysis
topology = structure_file
u = mda.Universe(topology, traj_file)

# Backup the original file
backup_traj_file = traj_file[:-4] + "_with_duplicates.dcd"
os.rename(traj_file, backup_traj_file)

with mda.Writer(traj_file, n_atoms=u.atoms.n_atoms) as W:
    for i in range(len(u.trajectory)):
        if i in indices_to_keep:
            u.trajectory[i]  # Set the trajectory frame to the given index
            W.write(u.atoms)

print("Done removing duplicate coordinates and energy reporters!")
