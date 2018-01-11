# Description: RNA-specific analyses


import os
import pandas as pd
import numpy as np


def make_rsem_matrix(quant_files, out_file, colname="expected_count"):
    """Given a list of RSEM files, extract the desired column and make matrix
    column for counts: expected_counts
    column for tpm: TPM
    """
    signal_mat = pd.DataFrame()
    for quant_file in quant_files:
        fields = quant_file.split("/")[-2].split(".")
        prefix = '{0}_{1}'.format(fields[0].split('-')[1], fields[4])
        
        data = pd.read_table(quant_file, sep="\t")
        signal_mat[prefix] = data[colname]
        
    signal_mat["gene_id"], _ = data["gene_id"].str.split(".", 1).str
    signal_mat = signal_mat.set_index("gene_id")
    if colname == "expected_count":
        signal_mat = signal_mat.astype(int)
    signal_mat.to_csv(out_file, compression="gzip", sep='\t')
    
    return None


def threshold_empirically_off_genes(mat_file, out_file, thresh=1.0):
    """Given a matrix file, get a distribution and determine empirical threshold
    """
    data = pd.read_table(mat_file, index_col=0)
    data_thresh = data.loc[~(data < thresh).all(axis=1)]

    data_thresh.to_csv(out_file, compression="gzip", sep='\t')
    
    return None
