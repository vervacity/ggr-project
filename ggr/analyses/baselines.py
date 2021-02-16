"""description: code for comparing analyses
"""

import os
import h5py

import numpy as np
import pandas as pd


def _get_sig_motif_list(pvals_file):
    """extract sig motifs
    """
    with h5py.File(pvals_file, "r") as hf:
        keys = sorted(hf["pvals"].keys())

    for key_idx in range(len(keys)):
        key = "pvals/{}/sig".format(keys[key_idx])

        with h5py.File(pvals_file, "r") as hf:
            key_sig = hf[key][:]
            
        if key_idx == 0:
            sig = key_sig
        else:
            sig = np.logical_or(sig, key_sig)
            
    return sig.astype(int)



def compare_sig_motifs(pvals_file_1, pvals_file_2, sig_names, out_dir, prefix):
    """given two pvals files, compare which motifs were found as sig
    and plot out as venn diagram
    """
    # for each pvals file, extract the list of sig motifs
    sig_motifs_1 = _get_sig_motif_list(pvals_file_1)
    sig_motifs_2 = _get_sig_motif_list(pvals_file_2)

    # save out together
    sig_compare = pd.DataFrame({
        sig_names[0]: sig_motifs_1,
        sig_names[1]: sig_motifs_2})
    data_file = "{}/{}.compare_sig.txt.gz".format(out_dir, prefix)
    sig_compare.to_csv(data_file, sep="\t", compression="gzip")

    # go to R to plot
    plot_file = "{}.plot.pdf".format(data_file.split(".txt")[0])
    r_cmd = "~/git/ggr-project/R/plot.compare_sig.venn.R {} {}".format(data_file, plot_file)
    print r_cmd
    os.system(r_cmd)
    
    return


def _get_best_corrs(pvals_file):
    """extract sig motifs
    """
    with h5py.File(pvals_file, "r") as hf:
        keys = sorted(hf["pvals"].keys())

    for key_idx in range(len(keys)):
        key = "pvals/{}/correlations".format(keys[key_idx])

        # get hgnc_ids and corr values
        with h5py.File(pvals_file, "r") as hf:
            hgnc_ids = hf["pvals/{}".format(keys[key_idx])].attrs["hgnc_ids"]
            key_corrs = hf[key][:]

        # make dataframe
        key_results = pd.DataFrame({key_idx: key_corrs}, index=hgnc_ids)

        # save together
        if key_idx == 0:
            results = key_results
        else:
            results = results.merge(key_results, how="outer", left_index=True, right_index=True)

    # get best results
    results = results.fillna(0)
    results = pd.DataFrame(
        {0: results.max(axis=1)})
            
    return results



def compare_tf_motif_corrs(corr_1, corr_2, sig_names, out_dir, prefix):
    """given two corr files, compare TFs for best correlations uncovered
    """
    # for each pvals file, extract the top correlations for each motif
    tf_motif_corrs_1 = _get_best_corrs(corr_1)
    tf_motif_corrs_2 = _get_best_corrs(corr_2)

    # merge
    results = tf_motif_corrs_1.merge(tf_motif_corrs_2, how="outer", left_index=True, right_index=True)
    results.columns = sig_names
    results = results.fillna(0)

    # save out and plot
    data_file = "{}/{}.compare_corr.txt.gz".format(out_dir, prefix)
    results.to_csv(data_file, sep="\t", compression="gzip")
    
    plot_file = "{}.plot.pdf".format(data_file.split(".txt")[0])
    r_cmd = "~/git/ggr-project/R/plot.compare_tf_motif_corr.R {} {}".format(data_file, plot_file)
    print r_cmd
    os.system(r_cmd)
    
    return 

