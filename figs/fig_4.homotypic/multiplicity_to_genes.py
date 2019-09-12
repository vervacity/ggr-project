#!/usr/bin/env python

import os
import re
import sys
import h5py

import numpy as np
import pandas as pd

from ggr.analyses.bioinformatics import run_gprofiler

from tronn.util.utils import DataKeys
from tronn.util.formats import array_to_bed

REMOVE_SUBSTRINGS = [
    "anatomical",
    "ameboidal",
    "animal organ",
    "multicellular organism",
    "cellular developmental",
    "tube",
    "regulation of",
    "embryonic",
    "cardiovascular",
    "angiogenesis",
    "blood vessel",
    "vasculature",
    "immune",
    "defense",
    "signaling",
    "response to",
    "movement of"]

REMOVE_EXACT_STRINGS = [
    "system process",
    "system development",
    "developmental process",
    "tissue development"]


def load_data_from_multiple_h5_files(h5_files, key):
    """convenience wrapper
    """
    key_data = []
    for h5_file in h5_files:
        with h5py.File(h5_file, "r") as hf:
            key_data.append(hf[key][:])
    key_data = np.concatenate(key_data, axis=0)
    
    return key_data


def bed_to_gene_set_by_proximity(
        bed_file,
        tss_file,
        out_file,
        k_nearest=3,
        max_dist=250000):
    """assuming proximal linking captures the majority of linking
    """
    # bedtools closest
    tmp_file = "tss.overlap.tmp.txt"
    closest = "bedtools closest -d -k {} -a {} -b {} > {}".format(
        k_nearest, bed_file, tss_file, tmp_file)
    os.system(closest)

    # load results and use distance cutoff
    data = pd.read_csv(tmp_file, sep="\t", header=None)
    data = data[data[9] < max_dist]

    # clean up
    data = data.sort_values([6, 9])
    data = data.drop_duplicates(6, keep="first")
    print "num genes:", len(set(data[6].values.tolist()))
    
    # split coumn and keep traj, gene name, distance
    data["gene_id"] = data[6]
    data["dist"] = data[9]
    gene_set = data[["gene_id", "dist"]]
    
    # save out
    gene_set.to_csv(out_file, compression="gzip", sep="\t", header=True, index=False)
    
    # cleanup
    os.system("rm {}".format(tmp_file))
    
    return data


def main():
    """get gene sets
    """
    # inputs
    OUT_DIR = sys.argv[1]
    sig_pwms_file = sys.argv[2]
    tss_file = sys.argv[3]
    background_gene_set_file = sys.argv[4]
    motifs_files = sys.argv[5:]
    density_key = DataKeys.ORIG_SEQ_PWM_MAX_DENSITIES
    metadata_key = DataKeys.SEQ_METADATA
    activity_keys = ["ATAC_SIGNALS.NORM", "H3K27ac_SIGNALS.NORM"]
    
    # make a tmp_dir
    tmp_dir = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(tmp_dir))
    
    # params
    POOL_WINDOW = 20 # from scanmotifs
    MAX_COUNT = 6
    count_thresholds = range(1, MAX_COUNT+1)
    plot_timepoint_indices = {
        "ATAC_SIGNALS.NORM": [0,6,12],
        "H3K27ac_SIGNALS.NORM": [0,1,2]} # d0,3,6
    tol = 1e-5
    MIN_REGION_COUNT = 500
    
    # load sig pwms
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", index_col=0, header=0).index)

    # load all pwm names and get indices for sig pwms
    with h5py.File(motifs_files[0], "r") as hf:
        all_pwms = hf[DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_IDX].attrs["pwm_names"]
    pwm_mask = np.array([1 if pwm_name in sig_pwms else 0 for pwm_name in all_pwms])
    sig_indices = np.where(pwm_mask != 0)[0]
    
    # load densities
    densities = load_data_from_multiple_h5_files(motifs_files, density_key)
    counts = densities * POOL_WINDOW

    # load metadata
    metadata = load_data_from_multiple_h5_files(motifs_files, metadata_key)[:,0]

    # analyze per pwm and save out
    keep_pwm_names = []
    for i in range(len(sig_indices)):
        sig_idx = sig_indices[i]
        pwm_name = all_pwms[sig_idx]
        pwm_name = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name = re.sub(".UNK.0.A", "", pwm_name)
        keep_pwm_names.append(pwm_name)
        print sig_idx, pwm_name

        if "TP53" not in pwm_name:
            continue
        
        # extract features
        pwm_counts = counts[:,0,sig_idx]

        # first mask values that aren't (mod 0.05),
        # they're edge values
        keep_mask = (np.mod(pwm_counts, 1) < tol).astype(int)
        #print np.sum(pwm_counts != 0)
        pwm_counts *= keep_mask
        #print np.sum(pwm_counts != 0)
        #quit()
        
        # now go through each level
        summary_df = None
        for count_threshold in count_thresholds:
            thresholded_counts = pwm_counts == count_threshold # vs greater than?
            thresholded_indices = np.where(thresholded_counts)[0]

            if thresholded_indices.shape[0] > MIN_REGION_COUNT:
                thresholded_metadata = metadata[thresholded_indices]
                print ">> num_regions:", thresholded_metadata.shape[0]
                # build the BED file
                bed_file = "{}/{}.counts.greater_equal-{}.bed.gz".format(
                    tmp_dir, pwm_name, count_threshold)
                array_to_bed(
                    thresholded_metadata,
                    bed_file, interval_key="active", merge=True)
                
                # then match to proximal gene set
                gene_set_file = "{}.gene_set.txt.gz".format(bed_file.split(".bed")[0])
                bed_to_gene_set_by_proximity(
                    bed_file,
                    tss_file,
                    gene_set_file,
                    k_nearest=2, # 2?
                    max_dist=25000)

                # and run gprofiler
                gprofiler_file = "{}.go_gprofiler.txt".format(gene_set_file.split(".txt.gz")[0])
                if not os.path.isfile(gprofiler_file):
                    run_gprofiler(
                        gene_set_file,
                        background_gene_set_file,
                        tmp_dir,
                        header=True)

                # read in gprofiler file and clean
                thresholded_summary = pd.read_csv(gprofiler_file, sep="\t")
                thresholded_summary = thresholded_summary[
                    (thresholded_summary["domain"] == "BP") |
                    (thresholded_summary["domain"] == "CC") |
                    (thresholded_summary["domain"] == "rea")]
                #thresholded_summary = thresholded_summary[thresholded_summary["domain"] == "rea"]
                #thresholded_summary = thresholded_summary[thresholded_summary["domain"] == "BP"]
                thresholded_summary = thresholded_summary[
                    ["term.id", "p.value", "term.name"]]
                thresholded_summary["count_threshold"] = count_threshold
                print "Num terms:", thresholded_summary.shape[0]

                # add to summary
                if summary_df is None:
                    summary_df = thresholded_summary.copy()
                else:
                    # filter first?
                    if thresholded_summary.shape[0] > 3:
                        summary_df = summary_df[
                            summary_df["term.id"].isin(thresholded_summary["term.id"].values)]
                    summary_df = pd.concat([summary_df, thresholded_summary], axis=0)
                    summary_df = summary_df.sort_values("term.name")

        # extra filtering?
        # keep only top 10 by p.value from max count?
        if False:
            keep_num = 100
            for count_threshold in reversed(count_thresholds):
                print count_threshold
                filter_data = summary_df[summary_df["count_threshold"] == count_threshold]
                if filter_data.shape[0] == 0:
                    continue
                filter_data = filter_data.sort_values("p.value", ascending=True)
                keep_terms = filter_data.iloc[0:keep_num]["term.id"].values
                summary_df = summary_df[summary_df["term.id"].isin(keep_terms)]
                break

        summary_df["log10pval"] = -np.log10(summary_df["p.value"].values)
        summary_df = summary_df.drop("p.value", axis=1)
        summary_file = "{}/{}.summary.txt.gz".format(tmp_dir, pwm_name)
        summary_df.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")
                
        # remove substrings
        for substring in REMOVE_SUBSTRINGS:
            keep = [False if substring in term else True
                    for term in summary_df["term.name"]]
            keep_indices = np.where(keep)[0]
            summary_df = summary_df.iloc[keep_indices]

        # remove exact strings
        keep = [False if term in REMOVE_EXACT_STRINGS else True
                for term in summary_df["term.name"]]
        keep_indices = np.where(keep)[0]
        summary_df = summary_df.iloc[keep_indices]
        summary_df = summary_df.sort_values(["term.name", "count_threshold"])
        summary_file = "{}/{}.summary.filt.txt.gz".format(tmp_dir, pwm_name)
        summary_df.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")

        
        # plot?
        # TODO think about how to filter down GO terms
        
        quit()

    
    return


main()
