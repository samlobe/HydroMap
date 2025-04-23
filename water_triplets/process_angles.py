#%%
import numpy as np
import pandas as pd
import MDAnalysis as mda
import argparse
import os
import sys
from collections import OrderedDict
from tqdm import tqdm

# Setting up argparse
parser = argparse.ArgumentParser(description='Read the triplet angles of your protein and output triplet distributions.')
# Add arguments
parser.add_argument('protein', help="unprocessed protein file ( e.g. <protein>[.pdb] )")
parser.add_argument('--anglesDir', type=str, default='angles', help="directory containing the angles files")
parser.add_argument('--multiChain', action='store_true', help="protein has multiple chains")
parser.add_argument('--groupsFile', type=str, help='to reference a file containing your custom groups (i.e. MDAnalysis selection strings, one group per line).')
parser.add_argument('--oneAnglesFile', type=str, help='to process data from one angles file (e.g. angles/<protein>_<selection>_angles.txt)')
parser.add_argument('-o','--outputFile', type=str, help="Specify a custom name for the output CSV file")

# Read arguments
args = parser.parse_args()
anglesFile = args.oneAnglesFile

# Count how many of the exclusive options are set
num_set = sum([args.multiChain, args.groupsFile is not None, args.oneAnglesFile is not None])
# Check if more than one of the exclusive arguments are set
if num_set > 1:
    parser.error("You can only use one of --multiChain, --groupsFile, or --oneAnglesFile, not multiple.")

# Assign arguments
protein = args.protein
if not protein.endswith('.pdb'):
    protein += '.pdb'
pdb_path = protein  # Uses the protein name from command line argument

# if the protein file doesn't exist, exit
if not os.path.exists(pdb_path):
    print(f"Error: Can't find {pdb_path}")
    sys.exit(1)

if 'processed' in pdb_path:
    print(f"Error: I see 'processed' in the protein file name. Please provide the unprocessed PDB file as input (i.e. no waters/ions).")
    sys.exit(1)

# extract protein name from the protein file name
protein_name = protein[:-4] # excluding the '.pdb' part
# check if the protein name ends with '_withH', if so remove it
if protein_name.endswith('_withH'):
    protein_name = protein_name[:-5]
# check if the protein name has '/' in it (e.g. '../myProtein')
# if so read the part after the last '/' (e.g. 'myProtein')
if '/' in protein_name:
    protein_name = protein_name.split('/')[-1]

# read in the protein file
u = mda.Universe(pdb_path)

# initalize histogramming of the water triplet distribution
bin_width = 5 # degrees
min = 40
max = 180
bins = np.arange(min,max+bin_width,bin_width)
def histo_line(data):
    histo_height, bin_edges = np.histogram(data, bins=bins, density=True)
    bin_middle = np.diff(bin_edges)/2
    bin_middle = bin_middle + bin_edges[:len(bin_middle)]
    return bin_middle, histo_height

# function to read the angles from the .txt files created by triplets.py
def read_angles(filepath):
    with open(filepath,'r') as f:
        angles = []
        for i,line in enumerate(f):
            if not line:
                print(f'Line {i} is empty.')
                angles.append([])
                continue
            else:
                angles.append([float(x) for x in line.split()])
                
    all_angles = np.array([item for sublist in angles for item in sublist])
    avg_measures = len(all_angles) / len(angles) # avg number of angle measurements per frame
    return all_angles, avg_measures 
#%%
# Mapping of 3-letter code to 1-letter code for amino acids
aa_mapping = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'HID': 'H', 'HIE': 'H', 'HIP': 'H', 'ACE': 'X',
    'NME': 'X', 'UNK': 'X'
}

def read_pdb_get_sequence(pdb_file):
    sequence = []
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
        residues = OrderedDict()
        for line in lines:
            if line.startswith("ATOM"):
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                # throw error if res_name = 'SOL'
                if res_name == 'SOL':
                    raise ValueError(f"Error: one or more of your selected atoms were waters. The input should be the unprocessed pdb (i.e. just protein without water or ions)")
                residues[res_num] = aa_mapping[res_name]
        sequence = list(residues.values())
    return ''.join(sequence)

# Read sequence from PDB
sequence = read_pdb_get_sequence(pdb_path)

# Write sequence to FASTA
with open(f'{protein_name}.fasta', 'w') as out:
    out.write(f'>{protein_name}\n')
    out.write(sequence + '\n')

#%% get residue/group data
script_dir = os.path.dirname(__file__) # absolute dir the script is in
group_distros = []
group_names = []
avg_num_angles_list = []
group_counter = 0

if anglesFile: # reading the angles from the angles file if given
    print(f"Processing your angles file {anglesFile} ...")
    # select after the protein name and before the '_angles.txt' part
    group_name = anglesFile.split('/')[-1][len(protein_name)+1:-11] # just one group name
    # replace the _ with a space in the group_name
    group_selection = group_name.replace('_',' ') # MDAnalysis selection string
    groups_selection = [group_selection] # just one entry
    group_names.append(group_name)
    group_angles,avg_num_angles = read_angles(anglesFile)
    bin_mids, histo = histo_line(group_angles)
    group_distros.append(histo)
    avg_num_angles_list.append(avg_num_angles)

elif args.groupsFile: # reading the custom groups in the groups file if given
    print("Processing your custom groups' data...")
    with open(args.groupsFile, 'r') as f:
        groups_selection = [line for line in f if not line.strip().startswith('#')]
    total_groups = len(groups_selection)
    for group_num in tqdm(np.arange(1,total_groups+1)):
        group_file = f'{script_dir}/{args.anglesDir}/{protein_name}_group{group_num}_angles.txt'
        if not os.path.exists(group_file):
            print(f'Group {group_num} data is missing.')
            continue
    
        group_names.append(f'group{group_num}')
        group_angles,avg_num_angles = read_angles(group_file)
        bin_mids, histo = histo_line(group_angles)
        group_distros.append(histo)
        avg_num_angles_list.append(avg_num_angles)
        group_counter += 1

# if no custom groups file, then use the resids (and segids)
else:
    print("Processing each residue's data ...")
    resids = u.residues.resids
    segids = u.residues.segids # used when there are multiple chains
    total_groups = len(resids)
    groups_selection = [] # for MDAnalysis selection strings
    for i, (resid, segid) in enumerate(zip(tqdm(resids), segids)):
        if args.multiChain:
            group_file = f'{script_dir}/{args.anglesDir}/{protein_name}_res{resid}_chain{segid}_angles.txt'
            if not os.path.exists(group_file):
                print(f'Data missing for: residue {resid} chain {segid}')
                continue
            seq_id = resid - u.residues.resids[0]
            group_names.append(f'{sequence[seq_id]}{resid}_chain{segid}')
            groups_selection.append(f'resid {resid} and segid {segid}') # not excluding hydrogens
        else:
            group_file = f'{script_dir}/{args.anglesDir}/{protein_name}_res{resid}_angles.txt'
            if not os.path.exists(group_file):
                print(f'Data missing for: residue {resid}')
                continue
            seq_id = resid - u.residues.resids[0]
            group_names.append(f'{sequence[seq_id]}{resid}')
            groups_selection.append(f'resid {resid}') # not excluding hydrogens

        group_angles,avg_num_angles = read_angles(group_file)
        bin_mids, histo = histo_line(group_angles)
        group_distros.append(histo)
        avg_num_angles_list.append(avg_num_angles)
        group_counter += 1

# Turn into DataFrame
avg_num_angles_list = np.array(avg_num_angles_list).reshape(-1, 1)
groups_selection = np.array(groups_selection).reshape(-1, 1)
groups_data = np.hstack((np.array(group_distros), avg_num_angles_list, groups_selection))
# create column names list
columns_after_bin_mids = np.array(['avg_residue_angles', 'MDAnalysis_selection_strings'])
columns = np.append(bin_mids,columns_after_bin_mids)
groups_df = pd.DataFrame(data=groups_data, index=group_names, columns=columns)
output_filename = args.outputFile if args.outputFile else f'{protein_name}_triplet_data.csv'
groups_df.to_csv(output_filename)
print(f'Successfully processed angles and outputted {output_filename}')

#%%
