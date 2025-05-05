#!/usr/bin/env python3
"""
color_by_property.py
────────────────────
Colour a PDB structure by one or more per-group properties
stored in a CSV that contains an MDAnalysis selection string
column called 'MDAnalysis_selection_strings'.

Example
-------
python color_by_property.py myProtein.pdb myProtein_results.csv
# writes:  myProtein_PC1_colored.pdb
#          myProtein_PC2_colored.pdb
#          myProtein_PC3_colored.pdb
#          myProtein_Fdewet_colored.pdb

Only colour by Fdewet_pred:
python color_by_property.py myProtein.pdb myProtein_results.csv \
        --properties Fdewet_pred
"""
import argparse, os, sys
import numpy as np, pandas as pd, MDAnalysis as mda

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('pdb',   help='input PDB file')
    p.add_argument('csv',   help='results CSV with MDAnalysis_selection_strings column')
    p.add_argument('--properties', nargs='+',
                   default=['PC1','PC2','PC3','Fdewet_pred','water_potential'],
                   help='property column(s) to map onto B-factors '
                        '(default: PC1 PC2 PC3 Fdewet_pred)')
    p.add_argument('--minWaters', type=float, default=5,
                   help='Only colour groups with avg_n_waters >= this value')
    p.add_argument('-o','--outdir', default='results',
                   help='directory for coloured PDBs (default: results)')
    p.add_argument('--pad', type=float, default=-999.0,
                   help='B-factor for atoms not in any selection (default: -999)')
    return p.parse_args()

def main():
    a = parse_args()
    if not os.path.exists(a.pdb): sys.exit(f'Cannot find PDB {a.pdb}')
    if not os.path.exists(a.csv): sys.exit(f'Cannot find CSV {a.csv}')

    df = pd.read_csv(a.csv)
    if 'MDAnalysis_selection_strings' not in df.columns:
        sys.exit("CSV must contain 'MDAnalysis_selection_strings' column.")

    # check avg_n_waters presence if needed
    if a.minWaters is not None and 'avg_n_waters' not in df.columns:
        sys.exit("CSV must contain 'avg_n_waters' column if --minWaters is used.")

    # filter by minimum avg_n_waters
    if a.minWaters is not None:
        before = len(df)
        df = df[df['avg_n_waters'] >= a.minWaters]
        print(f"Filtered groups: kept {len(df)}/{before} entries with >= {a.minWaters} avg hydration waters (within 4.25A of selection).")


    # if water_potential in properties, rename to U_pw
    if 'water_potential' in a.properties:
        a.properties = [prop if prop != 'water_potential' else 'U_pw' for prop in a.properties]

    # ensure requested properties exist
    missing = [prop for prop in a.properties if prop not in df.columns]
    if missing:
        sys.exit(f'Column(s) not found in CSV: {", ".join(missing)}')

    os.makedirs(a.outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(a.pdb))[0]

    # load structure once
    u = mda.Universe(a.pdb)

    for prop in a.properties:
        # reset all B-factors
        u.atoms.tempfactors = a.pad

        # assign property to requested selections
        for _,row in df.iterrows():
            sel = row['MDAnalysis_selection_strings']
            val = row[prop] if prop != 'water_potential' else row['U_pw'] # accepting both names

            atoms = u.select_atoms(sel)
            if len(atoms) == 0:
                print(f'Warning – selection "{sel}" matched 0 atoms.')
                continue
            atoms.tempfactors = np.round(val, 2)
        
        # rename "Fdewet_pred" to "Fdewet", and U_pw to water_potential
        label = prop.replace('Fdewet_pred', 'Fdewet').replace('U_pw', 'water_potential')
    
        out_pdb = os.path.join(a.outdir, f'{base}_{label}_colored.pdb')

        u.atoms.write(out_pdb)
        print('wrote', out_pdb)

if __name__ == '__main__':
    main()
