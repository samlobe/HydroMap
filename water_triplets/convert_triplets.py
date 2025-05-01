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

# measure the principal component contributions that describe how the solute's triplet distribution differs from bulk water's triplet distribution
# according to Robinson / Jiao's PCA analysis
# FF is a string of the force field name (e.g. 'a99SBdisp','a03ws','C36m, etc.)
# which should be in the index of the bulk_water_triplets.csv file
def get_PCs(groups_df,FF):
    triplet_distros = groups_df.loc[:,'42.5':'177.5']
    triplet_distros.columns = triplet_distros.columns.astype(float)

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

    return PCs_df


