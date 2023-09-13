#%%
import MDAnalysis as mda
import os
import subprocess
import argparse
import numpy as np
from glob import glob

# Check for the compiled fortran code in the current directory
matching_files = glob("waterlib*.so") + glob("waterlib*.pyd")
if not matching_files:
    raise FileNotFoundError("Could not find a compiled Fortran code in the current directory. (Looking for these patterns: waterlib*.so or waterlib*.pyd) \n"
                            "Please compile the Fortran code first. Example compilation using f2py (included with anaconda):`f2py -c -m waterlib waterlib.f90`\n"
                            "Please consult the tutorial (on GitHub) and leave an issue on the GitHub page if you have trouble compiling waterlib.f90's fortran code.\n")

# Setting up arguments parser
parser = argparse.ArgumentParser(description='Set groups in your protein to analyze their triplets.\nSubmit analysis jobs in batches.\nSelects all residues (assuming 1 chain) by default.')
# adding arguments
parser.add_argument('protein', help="unprocessed protein file (e.g. '../myProtein.pdb')")
parser.add_argument('trajectory', help="trajectory file (e.g. '../traj.dcd')")
parser.add_argument('--multiChain', action='store_true', help="protein has multiple chains")
parser.add_argument('--groupsFile', help='File containing MDAnalysis selection strings, one per line.')
parser.add_argument('-t','--timeLimit', type=int, default=10, help='time limit (minutes) for analysis')
args = parser.parse_args()

# make sure user uses (1) --multiChain or (2) --groupsFile or (3) neither
if args.multiChain and args.groupsFile:
    raise ValueError("You can't use both --multiChain and --groupsFile. Please choose one.")

# reading protein file name
protein = args.protein
if not protein.endswith('.pdb'):
    protein += '.pdb'
pdb_path = protein  # Uses the protein name from command line argument
protein_processed = f'{protein[:-4]}_processed.gro' # name of the processed protein (ignoring extension)

# check if the processed protein file exists
if not os.path.exists(f'{protein_processed}'):
    raise ValueError(f"Can't find the processed protein file: {protein_processed}")

# extract protein name from the protein file name
protein_name = protein[:-4] # excluding the '.pdb' part
# check if the protein name has '/' in it (e.g. '../myProtein')
# if so read the part after the last '/' (e.g. 'myProtein')
if '/' in protein_name:
    protein_name = protein_name.split('/')[-1]


trajectory = args.trajectory
# check if the trajectory file exists
if not os.path.exists(args.trajectory):
    raise ValueError(f"Can't find the trajectory file: {trajectory}")

# convert the time limit to hh:mm:ss format
timeLimit = f"{args.timeLimit//60:02d}:{args.timeLimit%60:02d}:00"

# read in the protein file
u = mda.Universe(pdb_path)

if args.groupsFile: # reading the custom groups in the groups file if given
    with open(args.groupsFile, 'r') as f:
        groups = [line for line in f if not line.strip().startswith('#')]
    total_groups = len(groups)
else: # if no groups file, then use the resids (and segids)
    resids = u.residues.resids
    segids = u.residues.segids # used when there are multiple chains
    total_groups = len(resids)

jobs_per_script = 40 # how many processors are on the computing node you're using (e.g. 40 for UCSB's Pod cluster)

header = f"""#!/bin/bash
#SBATCH -J {protein_name}_{{batch_num}}    # Job name
#SBATCH -o outLog.%j                            # Name of stdout output file
#SBATCH -e errLog.%j                            # Name of stderr error file
#SBATCH -p batch                                # Queue (partition) name for UCSB Pod cluster
#SBATCH -N 1                                    # Total # of nodes
#SBATCH -n {{jobs_in_script}}                   # Total # of mpi tasks
#SBATCH --cpus-per-task=1
#SBATCH -t {timeLimit}                         # Run time (hh:mm:ss)
# #SBATCH --mail-user=<yourEmail> # uncomment these two lines and include email if desired
# #SBATCH --mail-type=END,FAIL    

module list
pwd
date
conda activate hydrophobicity

"""

# function to make scripts to analyze each residue (assuming 1 chain)
def create_each_residue_script(batch_num, jobs_in_script, some_resids):
    with open(f'tripletsBatch{batch_num}.sh', 'w') as f:
        f.write(header.format(batch_num=batch_num, jobs_in_script=jobs_in_script))
        # loop through the groups
        for resid in some_resids:
            f.write(f"srun --ntasks 1 --exclusive -c 1 python triplet.py {protein_processed} {trajectory} -res {resid} &\n")
        f.write("wait\n")

 # function to make scripts to analyze each residue when there are multiple chains
def create_each_residue_multiChain_script(batch_num, jobs_in_script, some_resids, some_segids):
    with open(f'tripletsBatch{batch_num}.sh', 'w') as f:
        f.write(header.format(batch_num=batch_num, jobs_in_script=jobs_in_script))
        # loop through the groups
        for resid,segid in zip(some_resids,some_segids):
            f.write(f"srun --ntasks 1 --exclusive -c 1 python triplet.py {protein_processed} {trajectory} -res {resid} -ch {segid} &\n")
        f.write("wait\n")

 # function to make scripts to analyze custom groups
def create_custom_groups_script(batch_num, jobs_in_script, some_groupNums):
    with open(f'tripletsBatch{batch_num}.sh', 'w') as f:
        f.write(header.format(batch_num=batch_num, jobs_in_script=jobs_in_script))
        # loop through the groups
        for groupNum in some_groupNums: # groupNum is the index of the group in the groups file, using 1-indexing
            f.write(f"srun --ntasks 1 --exclusive -c 1 python triplet.py {protein_processed} {trajectory} --groupsFile {args.groupsFile} --groupNum {groupNum} &\n")
        f.write("wait\n")

### STEP THROUGH ALL THE GROUPS, ASSIGN TO BATCHES, AND CREATE BATCH ANALYSIS SCRIPTS ###
batch_num = 1
for group_index in range(0, total_groups, jobs_per_script):
    end_group_index = min(group_index+jobs_per_script, total_groups)
    jobs_in_script = end_group_index - group_index
    if args.groupsFile:
        create_custom_groups_script(batch_num, jobs_in_script, np.arange(group_index+1, end_group_index+1))
    elif args.multiChain:
        create_each_residue_multiChain_script(batch_num, jobs_in_script, resids[group_index:end_group_index],
                                              segids[group_index:end_group_index])
    else:
        create_each_residue_script(batch_num, jobs_in_script, resids[group_index:end_group_index])
    batch_num += 1

print(f"Generated {batch_num-1} bash scripts.")

# Execute the bash scripts using sbatch
for j in range(1, batch_num):
    subprocess.run(["sbatch", f"tripletsBatch{j}.sh"])
