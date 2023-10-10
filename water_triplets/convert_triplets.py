#%%
import numpy as np
import pandas as pd

# import the bulk water triplet distribution
df_bulk = pd.read_csv('bulk_water_triplets.csv',index_col=0)
df_bulk.columns = df_bulk.columns.astype(float)

# import the principal components from Robinson / Jiao's PCA analysis
# which describes the solute's triplet distribution subtracted by the bulk water triplet distribution
PCs = pd.read_csv('principalComps.csv',index_col=0)
PCs.columns = PCs.columns.astype(float)

# model that predicted singleAA dewetting free energies (according to INDUS)
# using just two parameters: 100-120 degrees and 45-50 degrees
def singleAA_dewet_fit(groups_df,FF): # for a99SBdisp
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
    if FF == 'a99SBdisp':
        m1 = -0.148; m2 = 1.538; b = 3.899
    elif FF == 'a03ws':
        m1 = -0.322; m2 = 0.295; b = 8.869
    elif FF == 'C36m':
        m1 = -0.417; m2 = 1.456; b = 7.429
    else: # throw error if FF isn't one of these 3
        raise ValueError(f"Error: {FF} isn't one of the 3 force fields we tested.\nExpecting 'a99SBdisp','a03ws',or 'C36m'.\nMaybe you meant to use the function otherFF_singleAA_dewet_fit?")

    dewet_series = m1 * tetrahedral + m2 * highcoord + b
    dewet_series = dewet_series * 300 * 0.008314 # convert to kJ/mol

    # delete rows that aren't solvated
    dewet_series = dewet_series[solvated_mask]

    # delete entries of MDAnalysis selection strings for unsolvated groups
    selection_strings = groups_df['MDAnalysis_selection_strings'][solvated_mask]

    # create new dataframe (Fdewet is in kJ/mol)
    dewet_df = pd.DataFrame({'Fdewet':dewet_series,'MDAnalysis_selection_strings':selection_strings})

    return dewet_df

# measure the principal component contributions that describe how the solute's triplet distribution differs from bulk water's triplet distribution
# according to Robinson / Jiao's PCA analysis
# FF is a string of the force field name (e.g. 'a99SBdisp','a03ws','C36m, etc.)
# which should be in the index of the bulk_water_triplets.csv file
def get_PCs(groups_df,FF):
    triplet_distros = groups_df.loc[:,'42.5':'177.5']
    triplet_distros.columns = triplet_distros.columns.astype(float)

    # get boolean array for groups that are solvated
    # i.e. each frame has ~10 angle measurements on average
    solvated_mask = groups_df['avg_residue_angles'] > 10

    # subtract bulk water triplet distribution
    try:
        bulk_distro = df_bulk.loc[FF]
    except KeyError:
        raise ValueError(f"Error: {FF}'s bulk triplet distribution cannot be found (looking in the index of bulk_water_triplets.csv).\nPlease add it to bulk_water_triplets.csv and try again.")

    triplet_distros_dBulk = triplet_distros.sub(bulk_distro,axis=1)

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

# if you want to predict hydrophobicity with a different force field, you can use this function
# this uses the model that collapsed all 3 force fields' dewetting free energies into one line
# who knows if it's accurate for your force field, but we got a Pearson's r > 0.96 and a R^2 > 0.92 for the 3 force fields we tested
def otherFF_singleAA_dewet_fit(groups_df,FF): # a FF other than a99SBdisp, a03ws, and CHARMM36m
    # get bulk water triplet distribution for the force field
    try:
        bulk_distro = df_bulk.loc[FF]
    except KeyError:
        raise ValueError(f"Error: {FF}'s bulk triplet distribution cannot be found (looking in the index of bulk_water_triplets.csv).\nPlease add it to bulk_water_triplets.csv and try again.")
    bulk_tetrahedral = bulk_distro.loc['102.5':'117.5'].sum() * 5 * 100
    bulk_highcoord = bulk_distro.loc['47.5'] * 5 * 100
    
    # get boolean array for groups that are solvated
    # i.e. each frame has ~10 angle measurements on average
    solvated_mask = groups_df['avg_residue_angles'] > 10

    # select tetrahedral signature (100-120 degrees)
    tetrahedral = groups_df.loc[:,'102.5':'117.5'].sum(axis=1) * 5 * 100
    # x5 for bin width (degrees); x100 to convert to percentage)
    #     
    # select highly coordinated signature (45-50 degrees)
    highcoord = groups_df.loc[:,'47.5'] * 5 * 100
    # x5 for bin width (degrees); x100 to convert to percentage)

    # calculate predicted dewetting free energy based on regression to singleAA INDUS data
    m1 = -0.62; m2 = 2.24; m3 = 0.11; m4 = -0.95; b = 2.54
    dewet_series = m1 * tetrahedral + m2 * highcoord + m3 * bulk_tetrahedral + m4 * bulk_highcoord + b # already in kJ/mol

    # delete rows that aren't solvated
    dewet_series = dewet_series[solvated_mask]

    # delete entries of MDAnalysis selection strings for unsolvated groups
    selection_strings = groups_df['MDAnalysis_selection_strings'][solvated_mask]

    # create new dataframe (Fdewet is in kJ/mol)
    dewet_df = pd.DataFrame({'Fdewet':dewet_series,'MDAnalysis_selection_strings':selection_strings})

    return dewet_df

