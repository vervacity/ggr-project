# description: basic useful analyses

import os

import numpy as np
import pandas as pd


def quantile_norm(array):
    """quantile norm on a numpy array
    """
    # set up sorted mean and originals to interpolate
    array_rank_sorted = np.sort(array, axis=0)
    rank_means = np.mean(array_rank_sorted, axis=1)
    
    # now use ranks to adjust values
    for i in range(array.shape[1]):
        array[:,i] = np.interp(
            array[:,i],
            array_rank_sorted[:,i],
            rank_means)

    return array


def build_id_matching_mat(
        primary_mat_file,
        secondary_mat_file,
        secondary_out_file,
        keep_cols=["Basal_Avg", "Diff_Avg"],
        primary_id_col=None,
        id_conversion_file=None):
    """
    """
    primary_data = pd.read_csv(primary_mat_file, sep="\t", index_col=0)
    if primary_id_col is not None:
        primary_data = primary_data.set_index(primary_id_col)
    secondary_data = pd.read_csv(secondary_mat_file, sep="\t")
    if secondary_data.columns[0].startswith("Unnamed"):
        secondary_data = pd.read_csv(secondary_mat_file, sep="\t", index_col=0)
    
    if id_conversion_file is not None:
        mappings = pd.read_csv(id_conversion_file, sep="\t")
        mappings = mappings[["ensembl_gene_id", "hgnc_symbol"]]
        mappings = mappings[~pd.isna(mappings["hgnc_symbol"])]
        mappings = mappings.set_index("hgnc_symbol")
        mappings = mappings.drop_duplicates("ensembl_gene_id")
        secondary_data = secondary_data.merge(mappings, left_index=True, right_index=True, how="left")
        secondary_data = secondary_data[~pd.isna(secondary_data["ensembl_gene_id"])]
        secondary_data = secondary_data.set_index("ensembl_gene_id")

    if len(keep_cols) > 0:
        secondary_data = secondary_data[keep_cols]

    # merge to match and fill na
    merged = primary_data.merge(secondary_data, left_index=True, right_index=True, how="left")
    merged = merged.fillna(0)

    if len(keep_cols) > 0:
        merged = merged[keep_cols]

    # and save out to new file
    merged.to_csv(secondary_out_file, sep="\t", compression="gzip")
    
    return None




def plot_PCA(mat_files, out_file):
    """
    """
    pca_cmd = "plot.pca.R {} {}".format(out_file, " ".join(mat_files))
    print pca_cmd
    os.system(pca_cmd)
    
    return


