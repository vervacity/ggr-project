#!/usr/bin/env python

import os
import h5py
import logging

import numpy as np
import pandas as pd

def max_from_zero_1d(array):
    """this doesn't exist in numpy
    """
    val = max(np.max(array), np.min(array), key=abs)
    
    return val


def main():
    """analyses
    """
    logger = logging.basicConfig(level=logging.INFO)
    
    # inputs
    h5_file = "/srv/scratch/dskim89/ggr/ggr.tronn.2019-06-29.ase/ggr.analyzevariants.h5"
    asATAC_file = "/srv/scratch/dskim89/ggr/ggr.atac.2019-06-06.ase/all.results_ALL.txt"
    atac_tasks = [0,1,2,3,4,5,6,9,10,12]

    # glenn list of validated SNPs
    glenn_list = [
        "rs2664280",
        "rs2349075",
        "rs4704864",
        "rs9262142",
        "rs10749669",
        "rs4240865",
        "rs4237547",
        "rs479844"]
    
    # load reference
    asATAC_data = pd.read_csv(asATAC_file, sep="\t", index_col=0) # rsid
    asATAC_data = asATAC_data[
        ["rsid", "sig", "beta_max"]]
    
    # extract relevant info from file
    results = {}
    with h5py.File(h5_file, "r") as hf:

        # ids
        results["rsid"] = hf["variants.id.string"][:,0]
        results["variants.sig"] = np.mean(hf["variants.sig"][:], axis=-1).astype(float)
        
        # important scores
        impt_scores = hf["variants.impt_score"][:]
        delta_impt_scores = impt_scores[:,0] - impt_scores[:,1]
        results["variant.impt_score"] = np.apply_along_axis(
            max_from_zero_1d, 1, delta_impt_scores)
        results["variant.impt_score.present"] = np.mean(impt_scores != 0, axis=(1,2)).astype(float)
        results["variant.impt_score.ref"] = np.any(impt_scores > 0, axis=(1,2)).astype(float)
        
        # sign flip - the scores go in opp directions
        neg_sign_flip = ((impt_scores[:,0] < 0) == (impt_scores[:,1] > 0)).astype(int)
        pos_sign_flip = ((impt_scores[:,0] > 0) == (impt_scores[:,1] < 0)).astype(int)
        sign_flip = np.logical_xor(neg_sign_flip, pos_sign_flip).astype(int)
        results["variant.sign_flip"] = np.mean(sign_flip != 0, axis=-1).astype(float)

        # only keep impt scores where there was a sign flip
        results["variant.impt_score"] = np.apply_along_axis(
            max_from_zero_1d, 1, delta_impt_scores)

        # logits
        logits = hf["logits.norm"][:][:,:,atac_tasks]
        delta_logits = logits[:,0] - logits[:,1]
        results["delta_logit"] = np.apply_along_axis(
            max_from_zero_1d, 1, delta_logits)


    # make dataframe, group and agg
    nn_data = pd.DataFrame(results)
    nn_data = nn_data.groupby(['rsid']).mean()
    
    # set up sign
    nn_data["sign"] = (nn_data["delta_logit"] > 0).astype(int)
    nn_data["sign"][nn_data["delta_logit"] < 0] = -1
    
    # merge
    data = asATAC_data.merge(nn_data, on="rsid", how="inner") # left

    # check anything in glenn validated snps?
    data_glenn = data[data["rsid"].isin(glenn_list)]

    # look at only those with importance scores
    # ie, when impt score present, do we separate sig from non sig well?
    imptscore_file = "results.w_imptscore.txt"
    data_imptscore = data[data["variant.impt_score.present"] > 0]
    keep_cols = ["rsid", "sig", "delta_logit", "variants.sig", "variant.impt_score.present"]
    data_imptscore = data_imptscore[keep_cols]
    data_imptscore.to_csv(imptscore_file, sep="\t", index=False)

    # and also look at just sig
    sig_file = "results.sig_only.txt"
    data_sig = data[data["sig"] == 1]
    keep_cols = ["rsid", "sig", "beta_max", "variants.sig", "variant.impt_score.present"]
    data_sig = data_sig[keep_cols]
    data_sig.to_csv(sig_file, sep="\t", index=False)

    return


main()
