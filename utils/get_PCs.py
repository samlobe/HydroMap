import pandas as pd

def get_PCs(groups_df, FF, df_bulk, PCs):
    """
    Compute the principal component contributions that describe how the solute's 
    triplet distribution differs from bulk water's, according to Robinson / Jiao's PCA.

    Parameters:
    -----------
    groups_df : pd.DataFrame
        DataFrame containing triplet distributions per group and a column 'MDAnalysis_selection_strings'.
    FF : str
        Name of the force field to use for selecting the appropriate bulk distribution.
    df_bulk : pd.DataFrame
        DataFrame of bulk water triplet distributions (index: force field names, columns: angles).
    PCs : pd.DataFrame
        DataFrame of principal components (index: ['PC1', 'PC2', 'PC3'], columns: angles).

    Returns:
    --------
    PCs_df : pd.DataFrame
        DataFrame containing PC1, PC2, PC3 values and corresponding MDAnalysis selection strings.
    """
    triplet_distros = groups_df.loc[:, '42.5':'177.5']
    triplet_distros.columns = triplet_distros.columns.astype(float)

    try:
        bulk_distro = df_bulk.loc[FF]
    except KeyError:
        raise ValueError(
            f"Error: {FF}'s bulk triplet distribution cannot be found "
            f"(looking in the index of df_bulk). Please add it and try again."
        )

    triplet_distros_dBulk = triplet_distros.sub(bulk_distro, axis=1)

    PC1 = triplet_distros_dBulk.dot(PCs.loc['PC1']) * 1000
    PC2 = triplet_distros_dBulk.dot(PCs.loc['PC2']) * 1000
    PC3 = triplet_distros_dBulk.dot(PCs.loc['PC3']) * 1000

    PCs_df = pd.DataFrame({'PC1': PC1, 'PC2': PC2, 'PC3': PC3})
    PCs_df['MDAnalysis_selection_strings'] = groups_df['MDAnalysis_selection_strings'].values

    return PCs_df
