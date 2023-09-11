#%%
import numpy as np
import pandas as pd

# import the bulk water triplet distribution
df_bulk = pd.read_csv('bulk_water_triplets.csv',index_col=0)
df_bulk.columns = df_bulk.columns.astype(float)
a99SBdisp_bulk = df_bulk.loc['a99SBdisp']

# import the principal components from Robinson / Jiao's PCA analysis
# which describes the solute's triplet distribution subtracted by the bulk water triplet distribution
PCs = pd.read_csv('principalComps.csv',index_col=0)
PCs.columns = PCs.columns.astype(float)

# model that predicted singleAA dewetting free energies (according to INDUS)
# using just two parameters: 100-120 degrees and 45-50 degrees
def singleAA_dewet_fit(groups_df):
    # get boolean array for groups that are solvated
    # i.e. each frame has ~10 angle measurements on average
    solvated_mask = groups_df['avg_residue_angles'] > 10

    # select tetrahedral signature (100-120 degrees)
    tetrahedral = groups_df.loc[:,'102.5':'117.5'].sum(axis=1) * 5 * 100
    # x5 for bin width (degrees); x100 to convert to percentage)
    
    # select highly coordinated signature (45-50 degrees)
    highcoord = groups_df.loc[:,'47.5'] * 5 * 100
    # x5 for bin width (degrees); x100 to convert to percentage)

    # calculate predicted dewetting free energy based on regression to singleAA INDUS data
    m1 = -0.148; m2 = 1.538; b = 3.899
    dewet_series = m1 * tetrahedral + m2 * highcoord + b
    dewet_series = dewet_series * 300 * 0.008314 # convert to kJ/mol

    # delete rows that aren't solvated
    dewet_series = dewet_series[solvated_mask]

    # delete entries of MDAnalysis selection strings for unsolvated groups
    selection_strings = groups_df['MDAnalysis_selection_strings'][solvated_mask]

    # create new dataframe (Fdewet is in kJ/mol)
    dewet_df = pd.DataFrame({'Fdewet':dewet_series,'MDAnalysis_selection_strings':selection_strings})

    return dewet_df

def get_PCs(groups_df):
    triplet_distros = groups_df.loc[:,'42.5':'177.5']
    triplet_distros.columns = triplet_distros.columns.astype(float)

    # get boolean array for groups that are solvated
    # i.e. each frame has ~10 angle measurements on average
    solvated_mask = groups_df['avg_residue_angles'] > 10

    # subtract bulk water triplet distribution
    triplet_distros_dBulk = triplet_distros.sub(a99SBdisp_bulk,axis=1)

    # calculate dot product with each PC to get the contribution of that PC to each groups triplet distribution
    PC1 = triplet_distros_dBulk.dot(PCs.loc['PC1'])*1000
    PC2 = triplet_distros_dBulk.dot(PCs.loc['PC2'])*1000
    PC3 = triplet_distros_dBulk.dot(PCs.loc['PC3'])*1000
    # scaled up by 1000 so that the differences are more intuitive to us humans

    # combine the PCs into a dataframe
    PCs_df = pd.DataFrame({'PC1':PC1,'PC2':PC2,'PC3':PC3})

    # add the MDAnalysis selection strings to the df
    PCs_df['MDAnalysis_selection_strings'] = groups_df['MDAnalysis_selection_strings']

    # delete rows that aren't solvated
    PCs_df = PCs_df[solvated_mask]
    return PCs_df

