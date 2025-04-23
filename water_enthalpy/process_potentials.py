#!/usr/bin/env python3
"""
Process potentials computed by potential.py and output average energies per group.

Example
-------
python process_potentials.py myProtein.pdb
python process_potentials.py myProtein.pdb --multiChain
python process_potentials.py myProtein.pdb --groupsFile myGroups.txt
python process_potentials.py myProtein.pdb --onePotentialsFile energies/myProtein_res45_energies.csv
"""

import numpy as np
import pandas as pd
import argparse
import os
import sys
from collections import OrderedDict
import MDAnalysis as mda
from tqdm import tqdm

# ----------------------------------------------------------------------
#  Argument parser
# ----------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Process potentials and output group averages.")
parser.add_argument("protein", help="Unprocessed protein file (e.g. myProtein.pdb)")
parser.add_argument("--energiesDir", type=str, default="energies", help="Directory containing the energies files")
parser.add_argument("--multiChain", action="store_true", help="Protein has multiple chains")
parser.add_argument("--groupsFile", type=str, help="File listing custom groups (MDAnalysis selection strings)")
parser.add_argument("--onePotentialsFile", type=str, help="Process just a single potentials CSV file")
parser.add_argument("-o", "--outputFile", type=str, help="Custom name for output CSV file")
args = parser.parse_args()

potentialsFile = args.onePotentialsFile

# Check exclusivity of arguments
num_set = sum([args.multiChain, args.groupsFile is not None, args.onePotentialsFile is not None])
if num_set > 1:
    parser.error("Use only one of --multiChain, --groupsFile, or --onePotentialsFile.")

# ----------------------------------------------------------------------
#  Setup
# ----------------------------------------------------------------------
protein = args.protein
if not protein.endswith('.pdb'):
    protein += '.pdb'

if not os.path.exists(protein):
    print(f"Error: Can't find {protein}")
    sys.exit(1)

protein_name = os.path.splitext(os.path.basename(protein))[0]
if protein_name.endswith('_withH'):
    protein_name = protein_name[:-5]

# for relative paths
script_dir = os.path.dirname(__file__)

# ----------------------------------------------------------------------
#  Sequence reader (only needed for labeling residues properly)
# ----------------------------------------------------------------------
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
    residues = OrderedDict()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                res_num = int(line[22:26].strip())
                res_name = line[17:20].strip()
                if res_name == 'SOL':
                    raise ValueError(f"Error: Your input file seems to include waters. Please use an unprocessed protein-only PDB.")
                residues[res_num] = aa_mapping.get(res_name, 'X')
    sequence = list(residues.values())
    return ''.join(sequence)

# Read sequence
sequence = read_pdb_get_sequence(protein)

# Write sequence to FASTA
with open(f'{protein_name}.fasta', 'w') as out:
    out.write(f'>{protein_name}\n')
    out.write(sequence + '\n')

# ----------------------------------------------------------------------
#  Processing potentials
# ----------------------------------------------------------------------
group_means = []
group_names = []

if potentialsFile:
    print(f"Processing {potentialsFile} ...")
    group_name = potentialsFile.split('/')[-1][len(protein_name)+1:-13]  # strip "proteinname_" and "_energies.csv"
    df = pd.read_csv(potentialsFile)
    avg_coulomb = df['coulomb'].mean()
    avg_lj = df['lj'].mean()
    avg_total = df['total'].mean()
    group_means.append((avg_coulomb, avg_lj, avg_total))
    group_names.append(group_name)

elif args.groupsFile:
    print("Processing your custom groups' potentials...")
    with open(args.groupsFile, 'r') as f:
        groups_selection = [line for line in f if not line.strip().startswith('#')]
    total_groups = len(groups_selection)
    for group_num in tqdm(np.arange(1, total_groups+1)):
        potentials_path = f'{script_dir}/{args.energiesDir}/{protein_name}_group{group_num}_energies.csv'
        if not os.path.exists(potentials_path):
            print(f"Group {group_num} data missing: {potentials_path}")
            continue
        df = pd.read_csv(potentials_path)
        avg_coulomb = df['coulomb'].mean()
        avg_lj = df['lj'].mean()
        avg_total = df['total'].mean()
        group_means.append((avg_coulomb, avg_lj, avg_total))
        group_names.append(f'group{group_num}')

else:
    print("Processing each residue's potentials...")
    u = mda.Universe(protein)
    resids = u.residues.resids
    segids = u.residues.segids
    for resid, segid in zip(tqdm(resids), segids):
        if args.multiChain:
            potentials_path = f'{script_dir}/{args.energiesDir}/{protein_name}_res{resid}_chain{segid}_energies.csv'
            if not os.path.exists(potentials_path):
                print(f"Data missing for: residue {resid} chain {segid}")
                continue
            seq_id = resid - u.residues.resids[0]
            group_names.append(f'{sequence[seq_id]}{resid}_chain{segid}')
        else:
            potentials_path = f'{script_dir}/{args.energiesDir}/{protein_name}_res{resid}_energies.csv'
            if not os.path.exists(potentials_path):
                print(f"Data missing for: residue {resid}")
                continue
            seq_id = resid - u.residues.resids[0]
            group_names.append(f'{sequence[seq_id]}{resid}')

        df = pd.read_csv(potentials_path)
        avg_coulomb = df['coulomb'].mean()
        avg_lj = df['lj'].mean()
        avg_total = df['total'].mean()
        avg_nwaters = df['n_waters'].mean()
        group_means.append((avg_coulomb, avg_lj, avg_total, avg_nwaters))

# ----------------------------------------------------------------------
#  Assemble dataframe
# ----------------------------------------------------------------------
columns = ['avg_coulomb (kJ/mol)', 'avg_lj (kJ/mol)', 'avg_total (kJ/mol)', 'avg_n_waters']
groups_df = pd.DataFrame(data=group_means, index=group_names, columns=columns)

output_filename = args.outputFile if args.outputFile else f'{protein_name}_potential_data.csv'
groups_df.to_csv(output_filename)
print(f'Successfully processed potentials and outputted {output_filename}')
