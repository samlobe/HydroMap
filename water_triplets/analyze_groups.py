#%%
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

# import the data
df_raw = pd.read_csv('residue_data.csv',index_col=0)

# filter out the rows where avg_residue_angles is <1
for res in df_raw.index:
    if df_raw['avg_residue_angles'][res] < 5:
        df_raw[res] = np.nan
        print(f'{res} has <5 measurement for frame so we will ignore.')
        # turn row into NaNs
        df_raw.loc[res,:] = np.nan

# load the 45-50 degree column (labelled 47.5)
hydrophilic_sig = df_raw['47.5'] * 5 * 100 # 5 is the bin width, x100 to get percentage

# sum the 4 columns (100-120°) to get the hydrophilic signal
hydrophobic_sig = df_raw.loc[:,'102.5':'117.5'].sum(axis=1) * 5 * 100 # 5 is the bin width, x100 to get percentage

# save these into a new dataframe
df = pd.DataFrame({'hydrophobic':hydrophobic_sig,'hydrophilic':hydrophilic_sig})

# calculate predicted dewetting free energy based on regression in regress.py in ~/Desktop/Research/forcefield_comp/
m1 = -0.148; m2 = 1.538; b = 3.899
df['F_dewetting'] = m1 * df['hydrophobic'] + m2 * df['hydrophilic'] + b
df['F_dewetting'] = df['F_dewetting'] * 300 * 0.008314 # convert to kJ/mol
 
# delete rows with NaNs
df = df.dropna()

# create list of labels from index names, but ignore the part after the underscore
labels = [res.split('_')[0] for res in df.index]

# graph the hydrophobic and hydrophilic signals as x and y-coordinates
plt.figure(figsize=(5,5))
plt.scatter(df['hydrophobic'],df['hydrophilic'])
#write the residue label for each point
for i,label in enumerate(labels):
    plt.annotate(label,(df['hydrophobic'].iloc[i],df['hydrophilic'].iloc[i]))

plt.xlabel('Hydrophobic Signature (%)',fontsize=12)
plt.ylabel('Hydrophilic Signature (%)',fontsize=12)
plt.title('Hydrophobin: each dot is a residue',fontsize=15)




#%% plot bars of the hydrophobic signal
plt.figure(figsize=(10,5))
plt.bar(df.index,df['hydrophobic'])
plt.xlabel('Residue Number',fontsize=12)
plt.ylabel('Hydrophobic Signature (%)',fontsize=12)
plt.title('Hydrophobin: Full Residues',fontsize=15)
plt.ylim(20,32)
plt.xticks(df.index,labels,rotation=90)

#%% plot bars of the hydrophilic signal
plt.figure(figsize=(10,5))
plt.bar(df.index,df['hydrophilic'])
plt.xlabel('Residue Number',fontsize=12)
plt.ylabel('Hydrophilic Signature (%)',fontsize=12)
plt.title('Hydrophobin: Full Residues',fontsize=15)
plt.ylim(0.4,1.4)
plt.xticks(df.index,labels,rotation=90)

#%% plot bars of the dewetting free energy
plt.figure(figsize=(10,5))
plt.bar(df.index,df['F_dewetting'])
plt.xlabel('Residue Number',fontsize=12)
plt.ylabel('Predicted Dewetting Free Energy (kJ/mol)',fontsize=12)
plt.title('Hydrophobin: Full Residues',fontsize=15)
plt.ylim(-0.2,3)
plt.xticks(df.index,labels,rotation=90)
# plot the horizontal line at 0
plt.axhline(y=0,color='k',linestyle='--')

# %% Create Pymol commands to color the residues based on the dewetting free energy
# pymol_object = 'hydrophobin'

# resids = [int(label[1:]) for label in labels]

# for i,resid in enumerate(resids):
#     # Note that residue indices start from 1 in your original script
#     dewetting_val = df['F_dewetting'].iloc[i]
    
#     # If it's the last residue (resid 70), exclude OC1 and OC2
#     if resid == 70:
#         print(f'alter resid {resid} and not (name OC1 or name OC2), b={dewetting_val:.3f}')
#     else:
#         print(f'alter resid {resid}, b={dewetting_val:.3f}')
        
# # For spectrum coloring, you need to specify the minimum and maximum values. 
# # I'm assuming they are the min and max of the 'F_dewetting' column. If not, you can replace them.
# low_val = df['F_dewetting'].min()
# high_val = df['F_dewetting'].max()

# # print(f'spectrum b, red_white_blue, {pymol_object}, minimum={low_val:.3f}, maximum={high_val:.3f}')
# print(f'spectrum b, red_white_blue, {pymol_object}, minimum=1, maximum=2')

# %% change beta factor with MDAnalysis
import MDAnalysis as mda

# make a dictionary to map three-letter code to one-letter code
aa_dict = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
           'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
           'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
           'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# list the residue numbers with sufficient hydration
resids = [int(label[1:]) for label in labels]
# load the structure
u = mda.Universe('omicron_RBD_withH.pdb')
for atom in u.atoms:
    print(f'{atom.resid}{atom.resname}{atom.id}')
    if atom.resid in resids:
        one_letter_code = aa_dict[atom.resname]
        atom.tempfactor = df['F_dewetting'][f'{one_letter_code}{atom.resid}']
    else:
        atom.tempfactor = -1

# save the structure
u.atoms.write('omicron_dewet_colored.pdb')


#%%