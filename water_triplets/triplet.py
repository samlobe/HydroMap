#%%
import numpy as np
import MDAnalysis as mda
import water_properties as wp
from tqdm import tqdm
import sys
import os
import argparse
import shutil

# Setting up argparse
parser = argparse.ArgumentParser(description='Process some protein-related inputs.\nExample usage: `python triplet.py ../myProtein_processed.gro ../traj.dcd -res 42`')

# Add arguments
parser.add_argument('protein',type=str,help="Processed protein structure file (e.g. '../myProtein_processed.gro')")
parser.add_argument('trajectory', type=str, help="Trajectory file name (e.g. '../traj.dcd')")
parser.add_argument('-res', '--resid', type=int, help='resid argument')
parser.add_argument('-ch', '--chain', help='segid (i.e. chainID) argument')
parser.add_argument('-t', '--time', type=float, default=5.0, help='Last x ns. Default is 5 ns.')
parser.add_argument('--groupsFile', type=str, help='File containing MDAnalysis selection strings, one per line.')
parser.add_argument('--groupNum', type=int, help='Line number in groupFile to use as the selection string.')
parser.add_argument('--hydrationCutoff', type=float, default=4.25, help='Cutoff distance (in Å) to define the hydration waters. Default is 4.25 Å.)')

args = parser.parse_args()

# Assign the arguments
protein_processed = args.protein
resid = args.resid
segid = args.chain
last_x_ns = args.time
hydrationCutoff = args.hydrationCutoff

# Ensure that the protein file ends with '.gro'
if not protein_processed.endswith('.gro'):
    protein_processed += '.gro'

# Ensure that the protein file exists
if not os.path.exists(protein_processed):
    print(f"Error: Can't find {protein_processed}")
    sys.exit(1)

# Ensure that the trajectory file exists
if not os.path.exists(args.trajectory):
    print(f"Error: Can't find {args.trajectory}")
    sys.exit(1)

structure_path = protein_processed  # Looking at structure file (usually in parent directory)
traj_path = args.trajectory

# Get the protein name from the protein file name
if protein_processed.endswith('_processed.gro'):
    protein_name = protein_processed[:-14] # excluding the '_processed.gro' part
else:
    protein_name = protein_processed[:-4] # excluding the '.gro' part

# check if the protein name has '/' in it (e.g. '../myProtein')
# if so read the part after the last '/' (e.g. 'myProtein')
if '/' in protein_name:
    protein_name = protein_name.split('/')[-1]

### SELECT YOUR GROUP OF INTEREST WITH MDANALYSIS SELECTION LANGUAGE \
if args.groupsFile and args.groupNum:
    # Load selection strings from the provided file
    try:
        with open(args.groupsFile, 'r') as f:
            group_strings = f.readlines()
        my_residue = group_strings[args.groupNum - 1].strip()  # 1-based index for user-friendliness
    except (IndexError, FileNotFoundError, IOError) as e:
        print(f"Error loading custom groups' selection strings from file: {e}")
        sys.exit(1)
elif args.resid:
    if args.chain:  # If a chain was provided
        my_residue = f'resid {resid} and segid {segid} and not name H*'  # selecting heavy atoms
    else:  # No chain provided
        my_residue = f'resid {resid} and not name H*'  # selecting heavy atoms without chain
else:
    print("Error: Please either provide the -res (and -ch flags) or the --groupFile and --groupNum flags.")
    sys.exit(1)

### LOAD THE MD TRAJECTORY
u = mda.Universe(structure_path,traj_path)
total_frames = len(u.trajectory)
timestep = u.trajectory.dt
# convert ns to ps and then divide by timestep to get number of frames
frames_to_load = int((last_x_ns * 1e3) / timestep)

### ERROR CHECKING YOUR SELECTION
try:
    res_group = u.select_atoms(my_residue)
    if len(res_group) == 0:
        raise ValueError
except ValueError:
    print(f"No atoms were selected. Please check your selection criteria.")
    sys.exit(1)

# throw an error if the user selected any atoms that are waters
resnames = [atom.resname for atom in u.select_atoms(my_residue)]
if 'SOL' in resnames:
    print(f"Error: one or more of your selected atoms were waters.")
    print("You probably meant to select protein atoms. Please check your selection criteria.")
    sys.exit(1)

resname = resnames[0] # get the residue name from the first atoms (assuming all atoms are from the same residue)
print(f'Looking at {len(res_group)} atoms from this residue:')
res = res_group[0].resname;  resid = res_group[0].resid
print(f'   {res}{resid}')

# Continue setting up analysis
waters = 'resname SOL and name O*' # looking at just water oxygens
BoxDims = u.dimensions[:3] # Å; needed for periodic boundary conditions
lowCut = 0
highCut = 3.5 # Å; cutoff distance to establish the neighbors of each hydration water

# Create a list of lists of 3-body angles.
# One sublist per configuration (i.e. per frame of trajectory)
# Each sublist contains the 3-body angle for the waters in the first shell around your group of interest.
angles_list = [[] for i in range(len(u.trajectory))]

# Create a checkpoint file to save progress (every 1000 frames)
if args.groupsFile and args.groupNum:
    checkpoint_filename = f'checkpoint{protein_name}_group{args.groupNum}_angles.txt'
elif args.chain:
    checkpoint_filename = f'checkpoint_{protein_name}_res{resid}_chain{segid}_angles.txt'
else:
    checkpoint_filename = f'checkpoint_{protein_name}_res{resid}_angles.txt'


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

### CALCULATE THE WATER TRIPLET ANGLES
start_frame_for_last_x_ns = max(0,total_frames - frames_to_load) # ensure it's not negative
for i,ts in tqdm(enumerate(u.trajectory[start_frame_for_last_x_ns:])):
    # SELECT HYDRATION WATERS
    shell_waters = u.select_atoms(f'({waters}) and around {hydrationCutoff} ({my_residue})') # 4.25Å is ~2 water layers from the residue
    subPos = shell_waters.positions # hydration water positions
    Pos = u.select_atoms(waters).positions # all water positions
    
    # MEAUSURE TRIPLET ANGLES
    triplet_angs, _ = wp.getTripletAngs(subPos,Pos,BoxDims,lowCut,highCut)
    #print(three_body_angs) # for debugging
    angles_list[i+start_frame] = list(np.around(triplet_angs,1))
    
    # Save checkpoint every 1000 frames
    if (i+1) % 1000 == 0:
        with open(checkpoint_filename, 'w') as txtfile:
            for j in range(i+start_frame+1):
                txtfile.write(f"{j}:{str(angles_list[j])}\n")

### SAVE FINAL RESULT
if args.groupsFile and args.groupNum:
    output_filename = f'{protein_name}_group{args.groupNum}_angles.txt'
elif args.chain != None:
    output_filename = f'{protein_name}_res{resid}_chain{segid}_angles.txt'
else:
    output_filename = f'{protein_name}_res{resid}_angles.txt'

with open(output_filename,'w') as txtfile:
    # each line contains the 3-body angles for one configuration (i.e. frame of trajectory)
    # in each line, there are triplet angles for each of the hydration waters
    for line in angles_list:
        txtfile.write(str(line)[1:-1].replace(",","").replace("\n","  ")+'\n') # string formatting

# Delete checkpoint file
if os.path.exists(checkpoint_filename):
    os.remove(checkpoint_filename)

# create directory 'angles' if it doesn't exist
if not os.path.exists('angles'):
    os.makedirs('angles')

# move output file to angles
shutil.move(output_filename,f'angles/{output_filename}')

# %%
