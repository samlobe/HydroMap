#!/usr/bin/env python3
"""
process_and_predict.py
────────────────────────────────
Processes water-triplet angle files + per-frame potential CSVs,
computes PCs via convert_triplets, builds the exact 3-feature
table your exported model expects, scales it with the saved
StandardScaler, and outputs per-group Fdewet predictions.

Outputs these CSVs to the output directory (e.g. ./results):
  angles_and_potentials_<protein>.csv
  <protein>_results.csv

Example usage:
python process_and_predict.py SARS2.pdb --anglesDir angles --potentialsDir
 energies --model ../models/Fdewet.joblib
"""
# %%
import argparse, os, sys, re, warnings, pickle, joblib
from collections import OrderedDict
import numpy as np, pandas as pd, MDAnalysis as mda
from tqdm import tqdm
import convert_triplets

BIN_WIDTH   = 5
ANGLE_MIN   = 40
ANGLE_MAX   = 180
BINS        = np.arange(ANGLE_MIN, ANGLE_MAX + BIN_WIDTH, BIN_WIDTH)
BIN_MIDS    = [f'{(BINS[i]+BINS[i+1])/2:.1f}' for i in range(len(BINS)-1)]
AA_MAP      = {k:v for k,v in
               [('ALA','A'),('CYS','C'),('ASP','D'),('GLU','E'),('PHE','F'),
                ('GLY','G'),('HIS','H'),('ILE','I'),('LYS','K'),('LEU','L'),
                ('MET','M'),('ASN','N'),('PRO','P'),('GLN','Q'),('ARG','R'),
                ('SER','S'),('THR','T'),('VAL','V'),('TRP','W'),('TYR','Y'),
                ('HID','H'),('HIE','H'),('HIP','H'),('ACE','X'),('NME','X'),
                ('UNK','X')]}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('protein', help='unprocessed protein-only PDB')
    p.add_argument('--anglesDir',     required=True)
    p.add_argument('--potentialsDir', required=True)
    p.add_argument('--model',         required=True,
                   help='Joblib/Pickle file with dict: {model, scaler, feature_order}')
    p.add_argument('--groupsFile', help='custom grouping file (one selection per line)')
    p.add_argument('--multiChain', action='store_true',
                   help="ensures each chain's residues are analyzed separately")
    p.add_argument('--mask', help='NumPy .npy frame-mask (bool/0-1)') # haven't tested yet
    p.add_argument('-o','--outdir', default='results')
    return p.parse_args()

def extract_resnum(s):
    m=re.search(r'(\d+)$',str(s)); return int(m.group(1)) if m else None

def histo_line(angles):
    if len(angles)==0 or np.all(np.isnan(angles)): return np.zeros(len(BIN_MIDS))
    h,_=np.histogram(angles,bins=BINS,density=True)
    return np.nan_to_num(h)

def flat_angles(path):
    per_frame=[]
    with open(path) as fh:
        for ln in fh:
            per_frame.append([float(x) for x in ln.split()] if ln.strip() else [])
    flat=np.array([a for sub in per_frame for a in sub])
    return flat, len(flat)/len(per_frame)

def pdb_sequence(path):
    residues=OrderedDict()
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('ATOM'):
                rn=int(ln[22:26]); res=ln[17:20].strip()
                if res=='SOL': raise ValueError('Input PDB contains water.')
                residues[rn]=AA_MAP.get(res,'X')
    return ''.join(residues.values()), list(residues.keys())

def main():
    a=parse_args(); os.makedirs(a.outdir, exist_ok=True)
    pdb=a.protein if a.protein.endswith('.pdb') else a.protein+'.pdb'
    if not os.path.exists(pdb): sys.exit(f'Cannot find {pdb}')
    prot=os.path.splitext(os.path.basename(pdb))[0].removesuffix('_withH')
    mask=np.load(a.mask).astype(bool) if a.mask else None

    # mutually exclusive grouping flags
    if sum([bool(a.groupsFile), a.multiChain]) > 1:
        sys.exit('Choose at most one of --groupsFile or --multiChain.')

    u=mda.Universe(pdb)
    seq,resnums=pdb_sequence(pdb)

    # determine selection strings
    if a.groupsFile:
        with open(a.groupsFile) as fh:
            sel=[ln.strip() for ln in fh if ln.strip() and not ln.startswith('#')]
    elif a.multiChain:
        sel=[f'resid {r} and segid {s}'
             for r,s in zip(u.residues.resids,u.residues.segids)]
    else:
        sel=[f'resid {r}' for r in u.residues.resids]

    # angles
    hrows = []
    for ix,s in enumerate(tqdm(sel,desc='Angles')):
        tag=f'group{ix+1}' if a.groupsFile else None
        path=( f'{a.anglesDir}/{prot}_{tag}_angles.txt' if tag else
               (f'{a.anglesDir}/{prot}_res{resnums[ix]}_chain{u.residues[ix].segid}_angles.txt'
                if a.multiChain else
                f'{a.anglesDir}/{prot}_res{resnums[ix]}_angles.txt') )
        if not os.path.exists(path): warnings.warn(f'{path} missing'); continue
        flat,avgN=flat_angles(path)
        if mask is not None: flat=flat[mask[:len(flat)]]
        hrows.append(histo_line(flat))

    hist_df=pd.DataFrame(hrows,index=sel,columns=BIN_MIDS)
    hist_df['MDAnalysis_selection_strings']=sel

    # add 10° bins needed by model
    for rng in [(40,50),(140,150)]:
        c1,f1=f'{rng[0]+2.5:.1f}', f'{rng[0]+7.5:.1f}'
        hist_df[f'{rng[0]}-{rng[1]}_tri']=hist_df[c1]+hist_df[f1]

    # PCs via convert_triplets
    pc_df=convert_triplets.get_PCs(hist_df,'a99SBdisp')  # PC1–3

    angles_df=hist_df.drop(columns=['MDAnalysis_selection_strings'])

    # potentials
    tots,nwats,idx=[],[],[]
    for ix,s in enumerate(tqdm(sel,desc='Potentials')):
        tag=f'group{ix+1}' if a.groupsFile else None
        pot=( f'{a.potentialsDir}/{prot}_{tag}_energies.csv' if tag else
              (f'{a.potentialsDir}/{prot}_res{resnums[ix]}_chain{u.residues[ix].segid}_energies.csv'
               if a.multiChain else
               f'{a.potentialsDir}/{prot}_res{resnums[ix]}_energies.csv') )
        if not os.path.exists(pot): warnings.warn(f'{pot} missing'); continue
        df=pd.read_csv(pot)
        if mask is not None: df=df[mask[:len(df)]]
        tots.append(df['total'].mean()); nwats.append(df['n_waters'].mean()); idx.append(s)

    pot_df=pd.DataFrame({'U_pw':tots,'avg_n_waters':nwats},index=idx)
    pot_df['MDAnalysis_selection_strings']=idx

    # Compile model inputs (needs two 10° + total_pot) 
    features=pd.merge(hist_df,pc_df[['PC1','PC2','PC3']],
                      left_index=True,right_index=True)
    features=pd.merge(features,
                      pot_df.drop(columns=['MDAnalysis_selection_strings']),
                      left_index=True,right_index=True)
    # output features df to angles_and_potentials_<protein>.csv
    features.to_csv(f'{a.outdir}/angles_and_potentials_{prot}.csv', index=False)
    
    # match the names the model may expect
    features['PC1_tri'] = features['PC1']

    # load model
    bundle=joblib.load(a.model)
    if isinstance(bundle, dict):
        mdl    = bundle['model']
        scaler = bundle['scaler']
        order  = bundle['feature_order']
    else:
        mdl    = bundle
        scaler = None
        order  = getattr(mdl,'feature_names_in_', features.columns)

    # ensure all features are present
    for f in order:
        if f not in features.columns:
            raise ValueError(f'Missing feature {f} in input data.')

    # ensure column order & scale
    X = features[order].copy()
    if scaler is not None:  X = scaler.transform(X)

    preds = mdl.predict(X)

    # compile results
    res=pd.DataFrame({'MDAnalysis_selection_strings':features['MDAnalysis_selection_strings'],
                      'PC1':features['PC1'],
                      'PC2':features['PC2'],
                      'PC3':features['PC3'],
                      'Fdewet_pred':preds,
                      'U_pw':features['U_pw'],
                      'avg_n_waters':features['avg_n_waters'],})
    
    # round PC1, PC2, PC3, and Fdewet_pred to 3 decimal places
    res['PC1'] = res['PC1'].round(3)
    res['PC2'] = res['PC2'].round(3)
    res['PC3'] = res['PC3'].round(3)
    res['Fdewet_pred'] = res['Fdewet_pred'].round(3)
    res['avg_n_waters'] = res['avg_n_waters'].round(1)
    res['U_pw'] = res['U_pw'].round(3)

    res.to_csv(f'{a.outdir}/{prot}_results.csv', index=False)
    print('Finished. Results written to', f'{a.outdir}/{prot}_results.csv')

if __name__ == '__main__':
    main()
