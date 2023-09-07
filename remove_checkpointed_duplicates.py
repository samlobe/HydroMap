import MDAnalysis as mda
import numpy as np
import os
import sys

topology_file = sys.argv[1]

# Step 1: Process the energies.log
times = [] # in ps

with open("energies.log", "r") as file:
    lines = file.readlines()
    for line in lines:
        if line.startswith("#"):  # Assuming lines starting with # are comments or headers
            continue
        # Assuming time step is the first column in the log file, adjust if not
        time = round(float(line.split()[1]))
        times.append(time)

# Backup the original file
os.rename("energies.log", "energies_with_duplicates.log")

# Step 2: Remove duplicate time steps and overwrite the energies.log
unique_times, indices_to_keep = np.unique(times, return_index=True)
print(indices_to_keep)
with open("energies.log", "w") as file:
    file.write(lines[0])  # Write the header
    for i in indices_to_keep:
        file.write(lines[i+1])

# Step 3 & 4: Process the traj.dcd with MDAnalysis
topology = topology_file
u = mda.Universe(topology, "traj.dcd")

# Backup the original file
os.rename("traj.dcd", "traj_with_duplicates.dcd")

with mda.Writer("traj.dcd", n_atoms=u.atoms.n_atoms) as W:
    for i in range(len(u.trajectory)):
        if i in indices_to_keep:
            u.trajectory[i]  # Set the trajectory frame to the given index
            W.write(u.atoms)

print("Done removing duplicate coordinates and energy reporters!")
