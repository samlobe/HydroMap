import numpy as np
import MDAnalysis as mda
import water_properties as wp
from tqdm import tqdm
import sys
import os

resnum = int(sys.argv[1]) # argv1 in `python backbone_triplet.py <argv1>` is set to the residue number variable
u = mda.Universe('hydrophobin_rect.tpr','hydrophobin_rect.xtc')

my_residue = f'resid {resnum} and not element H' # selecting heavy atoms

try:
    res_group = u.select_atoms(my_residue)
    if len(res_group) == 0:
        raise ValueError
except ValueError:
    print(f"No atoms were selected from residue {resnum}. Please check your selection criteria.")
    sys.exit(1)
resname = u.select_atoms(my_residue).resnames[0]

print(f'Looking at {len(res_group)} atoms from this residue:')
res = res_group[0].resname;  resid = res_group[0].resid
print(f'   {res}{resid}')

waters = 'resname SOL and element O' # looking at just water oxygens

BoxDims = u.dimensions[:3] # Å; needed for periodic boundary conditions
lowCut = 0
highCut = 3.5 # Å; cutoff distance to establish the neighbors of each hydration water

# Create a list of lists of 3-body angles.
# One sublist per configuration.
# Each sublist contains one 3-body angle per water in the first shell.
angles_list = [[] for i in range(len(u.trajectory))]

# Checkpoint filename
checkpoint_filename = f'checkpoint_res{resnum}_angles.txt'

# Load from checkpoint if exists
start_frame = 0
if os.path.exists(checkpoint_filename):
    with open(checkpoint_filename, 'r') as file:
        for line in file:
            frame, angles = line.split(":")
            angles_str = angles.strip()[1:-1]
            if angles_str:
                angles = [float(angle) for angle in angles_str.split(',')]
            else: angles = []
            angles_list[int(frame)] = angles
        start_frame = int(frame) + 1

for i,ts in tqdm(enumerate(u.trajectory[start_frame:])):
    # select hydration waters
    shell_waters = u.select_atoms(f'({waters}) and around 4.25 ({my_residue})') # 4.25Å is ~2 water layers from the residue
    subPos = shell_waters.positions # hydration water positions
    Pos = u.select_atoms(waters).positions # all water positions
    
    three_body_angs, _ = wp.getCosAngs(subPos,Pos,BoxDims,lowCut,highCut)
    #print(three_body_angs)
    angles_list[i+start_frame] = list(np.around(three_body_angs,1))
    
    # Save checkpoint every 1000 frames
    if (i+1) % 1000 == 0:
        with open(checkpoint_filename, 'w') as txtfile:
            for j in range(i+start_frame+1):
                txtfile.write(f"{j}:{str(angles_list[j])}\n")

# Save final result
with open(f'res{resnum}_angles.txt','w') as txtfile:
    # each line contains the 3-body angles for one configuration
    # in each line, there is one element per water in the first shell
    for line in angles_list:
        txtfile.write(str(line)[1:-1].replace(",","").replace("\n","  ")+'\n') # string formatting

