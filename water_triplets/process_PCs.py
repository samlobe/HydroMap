#%%
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import argparse
import convert_triplets
import MDAnalysis as mda

# Setting up argparse
parser = argparse.ArgumentParser(description="Generate dewetting free energy predictions and PC contributions from the triplet distributions, and output 'colored' pdb.\nExample usage: `python analyze_groups.py ../myProtein_withH.pdb`")
# Add arguments
parser.add_argument('protein', help="unprocessed protein pdb file to color, e.g. myProtein_withH.pdb\n(recommended to use a pdb of protein without solvent or ions, and with hydrogens present)")
parser.add_argument('--tripletFile', type=str, help="CSV file containing triplet data; will default to <protein_name>_triplet_data.csv")
parser.add_argument('--PCs_outputFile', type=str, help="CSV file to save the PCs to; will default to <protein_name>_PCs.csv")
args = parser.parse_args()

# Assign arguments
protein = args.protein
if not protein.endswith('.pdb'):
    protein += '.pdb'
pdb_path = protein  # Uses the protein name from command line argument

# if the protein file doesn't exist, exit
if not os.path.exists(pdb_path):
    print(f"Error: Can't find {pdb_path}")
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

# Find the data file with the triplet angles
if args.tripletFile: 
    triplet_file = args.tripletFile
else: 
    triplet_file = f'{protein_name}_triplet_data.csv'
groups_df = pd.read_csv(triplet_file,index_col=0)

PCs_df = convert_triplets.get_PCs(groups_df,'a99SBdisp')

if args.PCs_outputFile:
    PCs_df.to_csv(args.PCs_outputFile)
else:
    PCs_df.to_csv(f'{protein_name}_PCs.csv')

#%%
# convert triplet distribution to predicted dewetting free energy
# using the model from fitting to single amino acids
# dewet_df = convert_triplets.singleAA_dewet_fit(groups_df,'a99SBdisp')
# save the df
# dewet_df.to_csv(f'{protein_name}_Fdewet.csv')
# the columns are 'FDewet (kJ/mol)' and 'MDAnalysis_selection_string'

## modelling hydrophobicity in a03ws and CHARMM36m:
# dewet_df = convert_triplets.singleAA_dewet_fit(groups_df,'a03ws')
# dewet_df = convert_triplets.singleAA_dewet_fit(groups_df,'C36m')

## to apply our hydrophobicity model to other FFs, you should:
## 1) add that FFs bulk water triplet distribution to 'bulk_water_triplets.csv', and
## 2) use this function:
# dewet_df = convert_triplets.otherFF_singleAA_dewet_fit(groups_df,FF)

#%%
# convert triplet distribution to PC1, PC2, and PC3 contributions
# using Robinson / Jiao's PCs from their triplet distributions
# (they subtracted out the bulk water triplet distribution before getting PCs)
# PCs_df = convert_triplets.get_PCs(groups_df,'a99SBdisp')
# # save the df
# PCs_df.to_csv(f'{protein_name}_PCs.csv')
## the columns are 'PC1','PC2','PC3', and 'MDAnalysis_selection_string'

## if using other force fields:
# PCs_df = convert_triplets.get_PCs(groups_df,'C36m')
# PCs_df = convert_triplets.get_PCs(groups_df,'a03ws')
# PCs_df = convert_triplets.get_PCs(groups_df,FF) # if you added the FF's bulk water distribution to to 'bulk_water_triplets.csv'

#%%
# combine the dewet_df and PCs_df into one dataframe
# (this is just for convenience)

#%%
# function to color the protein based on a property (e.g. dewetting free energy or a PC)
# (actually the tempfactor for each atom is set, and later we will use ChimeraX or nglview to color the protein)
# def color_pdb(pdb_path, df_wProp_selecStr, property):
#     protein_name = pdb_path[:-4]
#     # check if the protein name has '/' in it (e.g. '../myProtein')
#     # if so read the part after the last '/' (e.g. 'myProtein')
#     if '/' in protein_name:
#        protein_name = protein_name.split('/')[-1] 

#     # load the structure
#     u = mda.Universe(pdb_path)
#     # initialize all tempfactors to -100
#     for atom in u.atoms:
#         atom.tempfactor = -100
#     # loop through the df and assign the property to the groups' atoms' tempfactors
#     for group in df_wProp_selecStr.index:
#         selection_string = df_wProp_selecStr['MDAnalysis_selection_strings'][group]
#         for atom in u.select_atoms(selection_string):
#             atom.tempfactor = np.around(df_wProp_selecStr[property][group],2)
#     # save the structure
#     u.atoms.write(f'{protein_name}_{property}_colored.pdb')
#     print(f'Outputted {protein_name}_{property}_colored.pdb')
#     # return the MDAnalysis universe object
#     return u

# u_Fdewet = color_pdb(pdb_path, dewet_df, 'Fdewet')
# u_PC1 =    color_pdb(pdb_path, PCs_df, 'PC1')
# u_PC2 =    color_pdb(pdb_path, PCs_df, 'PC2')
# u_PC3 =    color_pdb(pdb_path, PCs_df, 'PC3')

# #%% 
# # plot histograms of Fdewet, PC1, PC2, and PC3 in 2x2 grid
# fig, axes = plt.subplots(2,2,figsize=(7,7))
# axes[0,0].hist(dewet_df['Fdewet'],bins=20)
# axes[0,0].set_xlabel('Dewetting free energy (kJ/mol)',fontsize=14)
# axes[0,0].set_ylabel('Number of groups',fontsize=14)
# axes[0,1].hist(PCs_df['PC1'],bins=20)
# axes[0,1].set_xlabel('PC1',fontsize=14)
# axes[0,1].set_ylabel('Number of groups',fontsize=14)
# axes[1,0].hist(PCs_df['PC2'],bins=20)
# axes[1,0].set_xlabel('PC2',fontsize=14)
# axes[1,0].set_ylabel('Number of groups',fontsize=14)
# axes[1,1].hist(PCs_df['PC3'],bins=20)
# axes[1,1].set_xlabel('PC3',fontsize=14)
# axes[1,1].set_ylabel('Number of groups',fontsize=14)
# # make sure there aren't decimals in the y-ticks
# for ax in axes.flatten():
#     ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
# plt.tight_layout()
# plt.savefig(f'{protein_name}_histograms.png')
# print(f'Outputted {protein_name}_histograms.png')
# # plt.show() # uncomment if running this on your personal computer
# #%%
# # plot heatmaps of PC1 vs PC2, PC1 vs PC3, and PC2 vs PC3 in 1x3 grid
# fig, axes = plt.subplots(1,3,figsize=(15,5))
# axes[0].scatter(PCs_df['PC1'],PCs_df['PC2'])
# axes[0].set_xlabel('PC1',fontsize=14)
# axes[0].set_ylabel('PC2',fontsize=14)
# axes[1].scatter(PCs_df['PC1'],PCs_df['PC3'])
# axes[1].set_xlabel('PC1',fontsize=14)
# axes[1].set_ylabel('PC3',fontsize=14)
# axes[2].scatter(PCs_df['PC2'],PCs_df['PC3'])
# axes[2].set_xlabel('PC2',fontsize=14)
# axes[2].set_ylabel('PC3',fontsize=14)
# # label each point with the residue name (from groups_df.index)
# for i,res in enumerate(PCs_df.index):
#     axes[0].annotate(res,(PCs_df['PC1'].iloc[i],PCs_df['PC2'].iloc[i]))
#     axes[1].annotate(res,(PCs_df['PC1'].iloc[i],PCs_df['PC3'].iloc[i]))
#     axes[2].annotate(res,(PCs_df['PC2'].iloc[i],PCs_df['PC3'].iloc[i]))
# plt.tight_layout()
# plt.savefig(f'{protein_name}_PCs_2D.png')
# print(f'Outputted {protein_name}_PCs_2D.png')
# plt.show() # uncomment if running this on your personal computer
#%%
# now you can use ChimeraX or Pymol (or nglview) to color the protein
# With ChimeraX: open the pdb with your property-of-interest set to the tempfactor (i.e. bfactor)
# Example `color bfactor range 2.5,7 palette red-white-blue; color @@bfactor<-99 black`
# where 2.5 and 7 are the min and max values of the property (pick this based on the histograms)
# Then Go to `Tools -> Depiction -> Color Key` to add a key if you want. Label the values.
# Make a label with `2dlab text "Predicted Dewetting Free Energy"`. You can drag by selecting "Move Label" in the Right Mouse tab.

# With Pymol: open the pdb with your property-of-interest set to the tempfactor (i.e. bfactor)
# `show surface`
# `spectrum b, red_white_blue, minimum=2.5, maximum=7`
# where 2.5 and 7 are the min and max values of the property (pick this based on the histograms)
# `color black, b<-99`

# # With nglview:
# import nglview as nv
# # reset the view
# view = nv.show_mdanalysis(u_Fdewet)
# # Color the surface based on the β-factor with specific bounds
# view.add_surface("protein", color_scheme="bfactor", colorScale='rwb',colorDomain=[2.5, 7])  # color based on secondary structure
# view.display()
# # I can't get it to color the unsolvated atoms black. Maybe you can figure it out.

