#!/usr/bin/env python

# fig 3b prep
# description: compile data for plotting
# data: /srv/scratch/dskim89/ggr/ggr.tronn.2019-06-11.ase

import h5py
import logging

import numpy as np
import pandas as pd

from scipy.stats import pearsonr

from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve


def _select_best_importance_score(importances, method="2"):
    """acrorss tasks, choose best (or selectively index in)
    """
    if method == "1":
        # method 1: select the position with the STRONGEST importance score overall
        abs_max_between_alleles = np.max(np.abs(importances), axis=1) # {N, task}
        max_indices = np.argmax(abs_max_between_alleles, axis=-1) # {N}
    elif method == "2":
        # method 2: select position with biggest DIFFERENCE between alleles
        allele_diff = np.abs(importances[:,0] - importances[:,1]) # {N, task}
        max_indices = np.argmax(allele_diff, axis=-1)
        
    # method 3: just choose one timepoint
        
    # then select
    importances_best = importances[:,:,max_indices].diagonal(axis1=0, axis2=2)
    importances_best = importances_best.transpose()

    if False:
        importances_best = importances[:,:,0]
    
    # subtract alt from ref
    importances_best = importances_best[:,0] - importances_best[:,1]
    
    return importances_best


def auprc(labels, probs):
    """Wrapper for sklearn AUPRC
    Args:
      labels: 1D vector of labels
      probs: 1D vector of probs
    Returns:
      auprc: area under precision-recall curve
    """
    pr_curve = precision_recall_curve(labels, probs)
    precision, recall = pr_curve[:2]

    return auc(recall, precision)

def make_recall_at_fdr(fdr):
    """Construct function to get recall at FDR
    
    Args:
      fdr: FDR value for precision cutoff
    Returns:
      recall_at_fdr: Function to generate the recall at
        fdr fraction (number of positives
        correctly called as positives out of total 
        positives, at a certain FDR value)
    """
    def recall_at_fdr(labels, probs):
        pr_curve = precision_recall_curve(labels, probs)
        precision, recall = pr_curve[:2]
        return recall[np.searchsorted(precision - fdr, 0)]

    return recall_at_fdr


def _summarize_split(pos_df, neg_df):
    """summaries
    """
    logging.info(
        "neg -> sig fraction: {}".format(
            np.sum(neg_df["sig"]) / float(neg_df.shape[0])))
    logging.info(
        "pos -> sig fraction: {}".format(
            np.sum(pos_df["sig"]) / float(pos_df.shape[0])))
    logging.info(
        "neg -> beta (abs): {}".format(
            np.percentile(np.abs(neg_df["beta_max"]),[10,25,50,75,90])))
    logging.info(
        "pos -> beta (abs): {}".format(
            np.percentile(np.abs(pos_df["beta_max"]),[10,25,50,75,90])))
    logging.info(
        "neg -> delta logit (abs): {}".format(
            np.percentile(np.abs(neg_df["delta_logit"]),[10,25,50,75,90])))
    logging.info(
        "pos -> delta logit (abs): {}".format(
            np.percentile(np.abs(pos_df["delta_logit"]),[10,25,50,75,90])))
    
    return None


def main():
    """analyses
    """
    logger = logging.basicConfig(level=logging.INFO)
    
    # inputs
    h5_file = "ggr.analyzevariants.h5"
    variant_position = 79
    asATAC_file = "/srv/scratch/dskim89/ggr/ggr.atac.2019-06-06.ase/all.results_ALL.txt"
    window = 10
    
    # extract relevant info from file
    with h5py.File(h5_file, "r") as hf:

        # backcheck - make sure variants are in the right place
        if False:
            seq = hf["sequence.active"][:,:,:,variant_position]
            diff = seq[:,0] - seq[:,1]
            diff_sum = np.sum(np.abs(diff))
            assert diff_sum == 2*diff.shape[0]
        
        # get importances at base pair position, adjust as desired
        single_bp_importances = np.sum(
            hf["sequence-weighted.active"][:,:,:,variant_position], axis=-1) # {N, 2, task}
        single_bp_importances = _select_best_importance_score(single_bp_importances)
        logging.info("single bp impts: {}".format(single_bp_importances.shape))
        
        # get importances within window
        extend_len = window / 2
        window_importances = np.sum(
            hf["sequence-weighted.active"][
                :,:,:,variant_position-extend_len:variant_position+extend_len], axis=-1)
        window_importances = np.sum(window_importances, axis=3) # {N, 2, task}
        window_importances = _select_best_importance_score(window_importances)
        logging.info("window impts: {}".format(window_importances.shape))

        # get sig variants (as calculated by gaussian null)
        sig_by_delta_logit = hf["variants.sig"][:]
        sig_by_delta_logit = np.any(sig_by_delta_logit, axis=-1)
        logging.info("sig by delta logit scoring: {}".format(sig_by_delta_logit.shape))
        
        # get logits NOTE: may need to be matched with importances
        atac_tasks = [0,1,2,3,4,5,6,9,10,12]
        logits = hf["logits.norm"][:] # {N, 2, task}
        logits = logits[:,:,atac_tasks]
        logging.info("logits: {}".format(logits.shape))
        delta_logit = logits[:,0] - logits[:,1]
        delta_logit = np.max(np.abs(delta_logit), axis=-1)

        # get sign (based on output)
        predicted_sign = np.mean(logits[:,0] - logits[:,1], axis=-1)
        predicted_sign = logits[:,0,1] - logits[:,1,1]
        #for key in sorted(hf.keys()): print key, hf[key].shape
        #predicted_sign = window_importances
        
        # split by whether there's an importance score or not
        single_impt_present = (single_bp_importances != 0).astype(int)
        window_impt_present = (window_importances != 0).astype(int)

        # get variant id
        variant_ids = hf["variants.id.string"][:,0]
        
    # combine all info into a df
    nn_scoring = pd.DataFrame({
        "rsid": variant_ids,
        "single_importance": single_bp_importances,
        "single_importance_present": single_impt_present,
        "window_importance": window_importances,
        "window_importance_present": window_impt_present,
        "delta_logit": delta_logit,
        "sign": predicted_sign,
        "sig_by_delta_logit": sig_by_delta_logit,
    })

    # check
    #print nn_scoring
    logging.info(
        "variants with importance score: {}".format(
            np.sum(nn_scoring["single_importance_present"].values)))
        
    # load reference and simplify
    asATAC_data = pd.read_csv(asATAC_file, sep="\t", index_col=0) # rsid
    asATAC_data = asATAC_data[
        ["rsid", "sig", "beta_max"]]
    logging.info("all asATAC: {}".format(asATAC_data.shape))
    logging.info("all asATAC, median beta: {}".format(
        np.median(np.abs(asATAC_data["beta_max"].values))))
    
    # merge in data
    data = asATAC_data.merge(nn_scoring, on="rsid", how="left")
    logging.info("all data: {}".format(data.shape))
    logging.info("asATAC sig: {}".format(
        np.sum(data["sig"])))

    if True:
        data["delta_logit"] = np.power(
            2, data["delta_logit"], where=data["delta_logit"]!=0)
    
    # across everything, if split by importance what are sig fractions and betas
    neg_df = data[data["single_importance"]==0]
    pos_df = data[data["single_importance"]!=0]
    logging.info(">> ALL: not impt vs single impt")
    _summarize_split(pos_df, neg_df)

    neg_df = data[data["window_importance"]==0]
    pos_df = data[data["window_importance"]!=0]
    logging.info(">> ALL: not impt vs window impt")
    _summarize_split(pos_df, neg_df)

    # now just look at the sig asATAC variants (atacQTLs)
    data_filt = data[data["sig"]!=0]
    neg_df = data_filt[data_filt["single_importance"]==0]
    pos_df = data_filt[data_filt["single_importance"]!=0]
    logging.info(">> ataQTLS: not impt vs single impt")
    _summarize_split(pos_df, neg_df)

    neg_df = data_filt[data_filt["window_importance"]==0]
    pos_df = data_filt[data_filt["window_importance"]!=0]
    logging.info(">> ataQTLS: not impt vs window impt")
    _summarize_split(pos_df, neg_df)

    # split by importance AND sig on output?
    
    

    # check AUPR (on predicting atacQTLs)
    if False:
        key = "window_importance"
        key = "single_importance"
        #key = "beta_max"
        key = "sig_by_delta_logit"
        key = "delta_logit"
        labels = data["sig"].values
        probs = np.abs(data[key])
        print np.sum(data["sig"])
        print data.shape
        print auprc(labels, probs)

        recall_at_fdr_fn = make_recall_at_fdr(0.50)
        print recall_at_fdr_fn(labels, probs)
    
    # try correlations
    logging.info("compare beta max to single impt: {}".format(
        pearsonr(data_filt["beta_max"], data_filt["single_importance"])))
    logging.info("compare beta max to single impt: {}".format(
        pearsonr(data_filt["beta_max"], np.abs(data_filt["single_importance"]))))

    logging.info("compare beta max to window impt: {}".format(
        pearsonr(data_filt["beta_max"], data_filt["window_importance"])))
    logging.info("compare beta max to window impt: {}".format(
        pearsonr(data_filt["beta_max"], np.abs(data_filt["window_importance"]))))
    
    logging.info("compare beta max to delta logit: {}".format(
        pearsonr(data_filt["beta_max"], data_filt["delta_logit"])))
    logging.info("compare beta max to delta logit: {}".format(
        pearsonr(data_filt["beta_max"], np.abs(data_filt["delta_logit"]))))

    # save out
    # for output, need:
    # sig, impt category, delta logit, beta
    #data["window_category"] = data["window_importance_present"].astype(str)
    #data.loc[data["window_importance_present"] == 0, "window_category"] = "no_impt"
    #data.loc[data["window_importance_present"] == 1, "window_category"] = "window_impt"

    #data["single_bp_category"] = data["single_importance_present"].astype(str)
    #data.loc[data["single_importance_present"] == 0, "single_bp_category"] = "no_single_impt"
    #data.loc[data["single_importance_present"] == 1, "single_bp_category"] = "single_impt"

    # make group category
    data["group"] = data["single_importance_present"].astype(str)
    data.loc[data["window_importance_present"] == 0, "group"] = "no_impt"
    data.loc[data["single_importance_present"] == 1, "group"] = "single_impt"
    data = data[data["group"]!="0"]
    window_data = data[data["window_importance_present"] == 1]
    window_data["group"] = "window_impt"
    data = pd.concat([data, window_data], axis=0)
    
    
    #data["groups"] = data["window_category"].map(str) + "_and_" + data["single_bp_category"]
    #data = data[data["groups"] != "window_input_and_no_single_impt"]
    #data.loc[data["groups"] == "no_impt_and_no_single_impt" ,"groups"] = "no_impt"
    #data.loc[data["groups"] == "window_impt_and_no_single_impt" ,"groups"] = "no_impt"
    
    data["sign"][data["sign"] < 0] = -1
    data["sign"][data["sign"] > 0] = 1
    
    data["delta_logit"] = np.abs(data["delta_logit"])
    data["beta_max_raw"] = data["beta_max"].values
    data["beta_max"] = np.abs(data["beta_max"])
    
    #data["sig"]
    data.loc[data["sig"] == 0, "sig"] = "not_sig"
    data.loc[data["sig"] == 1, "sig"] = "sig"

    #data["sig_w_group"] = data["sig"] + "_and_" + data["group"]
    
    save_data = data[[
        "sig",
        "group",
        "delta_logit",
        "sign",
        "beta_max",
        "beta_max_raw",
        "sig_by_delta_logit",
        "single_importance"]]
    save_data.to_csv("results.txt", sep="\t", index=False)
    
    
    print data_filt.columns
    
    return


main()
