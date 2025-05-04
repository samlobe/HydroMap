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
"""
# %%
import argparse, os, sys, re, warnings, pickle, joblib
from collections import OrderedDict
import numpy as np, pandas as pd, MDAnalysis as mda
from tqdm import tqdm
try:
    from get_PCs import get_PCs
except ImportError:
    sys.exit("ERROR: Cannot import get_PCs.py. Make sure it is in the same directory as the process_and_predict.py script (consider using a symbolic link: `ln -s /path/to/get_PCs.py .`).")

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
    p.add_argument('-ff','--forcefield', default='a99SBdisp',
                   help='forcefield for convert_triplets and read bulk (default: a99SBdisp)')
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

def main():
    a=parse_args(); os.makedirs(a.outdir, exist_ok=True)

    model_name = os.path.basename(a.model)
    if model_name == 'Fdewet.joblib' and a.forcefield != 'a99SBdisp':
        sys.exit('Fdewet.joblib model only trained with a99SBdisp forcefield.')
    
    if model_name == 'Fdewet_isolated_aa_multi_forcefield.joblib' and a.forcefield not in ['a99SBdisp', 'a03ws','C36m']:
        sys.exit('Only these forcefields are currently supported for single aa model: a99SBdisp, a03ws, C36m.')

    pdb=a.protein if a.protein.endswith('.pdb') else a.protein+'.pdb'
    if not os.path.exists(pdb): sys.exit(f'Cannot find {pdb}')
    prot=os.path.splitext(os.path.basename(pdb))[0].removesuffix('_withH')
    mask=np.load(a.mask).astype(bool) if a.mask else None

    if sum([bool(a.groupsFile), a.multiChain]) > 1:
        sys.exit('Choose at most one of --groupsFile or --multiChain.')

    u=mda.Universe(pdb)

    if a.groupsFile:
        with open(a.groupsFile) as fh:
            groups=[ln.strip() for ln in fh if ln.strip() and not ln.startswith('#')]
    elif a.multiChain:
        groups=[f'resid {r} and segid {s}' for r,s in zip(u.residues.resids,u.residues.segids)]
    else:
        groups=[f'resid {r}' for r in u.residues.resids]

    # angles
    hrows = []
    for i,group in enumerate(tqdm(groups,desc='Angles')):
        residue = u.residues[i]
        resid = residue.resid
        segid = residue.segid.strip()

        tag=f'group{i+1}' if a.groupsFile else None
        path=( f'{a.anglesDir}/{prot}_{tag}_angles.txt' if tag else
               (f'{a.anglesDir}/{prot}_res{resid}_chain{segid}_angles.txt'
                if a.multiChain else
                f'{a.anglesDir}/{prot}_res{resid}_angles.txt') )
        if not os.path.exists(path): warnings.warn(f'{path} missing'); continue
        flat,_=flat_angles(path)
        if mask is not None: flat=flat[mask[:len(flat)]]
        hrows.append(histo_line(flat))

    hist_df=pd.DataFrame(hrows,index=groups,columns=BIN_MIDS)
    hist_df['MDAnalysis_selection_strings']=groups

    # add 10-degree bins needed by model
    for rng in [(40,50),(140,150)]:
        c1,f1=f'{rng[0]+2.5:.1f}', f'{rng[0]+7.5:.1f}'
        hist_df[f'{rng[0]}-{rng[1]}_tri']=hist_df[c1]+hist_df[f1]

    # get PCs from triplet angles
    if not os.path.exists('data/bulk_water_triplets.csv'):
        sys.exit('Missing bulk_water_triplets.csv (try: ln -s <data_dir>)')
    if not os.path.exists('data/principalComps.csv'):
        sys.exit('Missing principalComps.csv (try: ln -s <data_dir>)')

    df_bulk = pd.read_csv('data/bulk_water_triplets.csv',index_col=0)
    df_bulk.columns = df_bulk.columns.astype(float)

    PCs = pd.read_csv('data/principalComps.csv',index_col=0)
    PCs.columns = PCs.columns.astype(float)

    pc_df = get_PCs(hist_df, a.forcefield, df_bulk, PCs)
    angles_df=hist_df.drop(columns=['MDAnalysis_selection_strings'])

    # potentials
    tots,nwats,idx=[],[],[]
    for i,group in enumerate(tqdm(groups,desc='Potentials')):
        residue = u.residues[i]
        resid = residue.resid
        segid = residue.segid.strip()

        tag=f'group{i+1}' if a.groupsFile else None
        pot=( f'{a.potentialsDir}/{prot}_{tag}_potentials.csv' if tag else
              (f'{a.potentialsDir}/{prot}_res{resid}_chain{segid}_potentials.csv'
               if a.multiChain else
               f'{a.potentialsDir}/{prot}_res{resid}_potentials.csv') )
        if not os.path.exists(pot): warnings.warn(f'{pot} missing'); continue
        df=pd.read_csv(pot)
        if mask is not None: df=df[mask[:len(df)]]
        tots.append(df['total'].mean()); nwats.append(df['n_waters'].mean()); idx.append(group)

    pot_df=pd.DataFrame({'U_pw':tots,'avg_n_waters':nwats},index=idx)
    pot_df['MDAnalysis_selection_strings']=idx

    features=pd.merge(hist_df,pc_df[['PC1','PC2','PC3']],
                      left_index=True,right_index=True)
    features=pd.merge(features,
                      pot_df.drop(columns=['MDAnalysis_selection_strings']),
                      left_index=True,right_index=True)
    features.to_csv(f'{a.outdir}/angles_and_potentials_{prot}.csv', index=False)

    if model_name == 'Fdewet_isolated_aa_multi_forcefield.joblib':
        features['PC1_tri'] = features['PC1']
        if a.forcefield == 'a99SBdisp':
            features['bulk_80-90_tri'] = 0.023396
        elif a.forcefield == 'a03ws':
            features['bulk_80-90_tri'] = 0.023285
        elif a.forcefield == 'C36m':
            features['bulk_80-90_tri'] = 0.021381
        else:
            raise ValueError(f'Unknown forcefield {a.forcefield}.')

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

    for f in order:
        if f not in features.columns:
            raise ValueError(f'Missing feature {f} in input data.')

    # ensure all features are present
    X = features[order].copy()
    if scaler is not None:  X = scaler.transform(X)

    print("Predicting Fdewet from protein-water potential and triplet angle features (assuming you used the default 4.25A cutoff).")
    preds = mdl.predict(X)

    # output results
    res=pd.DataFrame({'MDAnalysis_selection_strings':features['MDAnalysis_selection_strings'],
                      'PC1':features['PC1'],
                      'PC2':features['PC2'],
                      'PC3':features['PC3'],
                      'Fdewet_pred':preds,
                      'U_pw':features['U_pw'],
                      'avg_n_waters':features['avg_n_waters'],})
    
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
