#!/usr/bin/env python

import os
import re
import sys
import h5py

import numpy as np
import pandas as pd

from tronn.util.utils import DataKeys


def nonzero_colmeans(array):
    """don't consider zeros when calculating mean
    """
    value_sums = np.sum(array, axis=0)
    presence_sums = np.sum(array!=0, axis=0)
    distr = np.true_divide(
        value_sums, presence_sums, where=presence_sums!=0)
    
    return distr


def main():
    """analyze spacing
    """
    # ==============================
    # SETUP
    # ==============================
    
    # inputs
    SCRIPT_DIR = sys.argv[1]
    OUT_DIR = sys.argv[2]
    motifs_file = sys.argv[3]
    sig_pwms_file = sys.argv[4]

    # set up dirs
    os.system("mkdir -p {}".format(OUT_DIR))
    TMP_DIR = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(TMP_DIR))
    MORE_PLOTS_DIR = "{}/plots".format(OUT_DIR)
    os.system("mkdir -p {}".format(MORE_PLOTS_DIR))

    # keys
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    max_idx_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_IDX
    pwm_scores_key = DataKeys.WEIGHTED_SEQ_PWM_SCORES_THRESH
    raw_pwm_scores_key = DataKeys.ORIG_SEQ_PWM_SCORES_THRESH
    signal_keys = ["ATAC_SIGNALS", "H3K27ac_SIGNALS"]
    
    # params
    left_clip = 420 + 12
    num_positions = 200
    final_extend_len = 160

    # ==============================
    # LOAD DATA
    # ==============================
    
    # read in the sig pwms
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", header=0).iloc[:,0])

    # read in the list of all pwm names
    with h5py.File(motifs_file, "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    
    # read in max vals
    with h5py.File(motifs_file, "r") as hf:
        max_vals = hf[max_val_key][:] # {N, pwm, 1}

    # ==============================
    # ANALYZE
    # ==============================
        
    # analyze each pwm
    for pwm_idx in range(len(sig_pwms)):
        pwm_name = sig_pwms[pwm_idx]
        pwm_name_clean = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name_clean = re.sub(".UNK.0.A", "", pwm_name_clean)
        pwm_global_idx = np.where(
            [1 if pwm_name in global_name else 0
             for global_name in all_pwms])[0][0]
        print pwm_name_clean, pwm_global_idx

        # skip
        #continue
        #if "GRHL" not in pwm_name_clean:
        #    continue
        
        # check to see which have this pwm and adjust indices
        example_indices = np.where(max_vals[:,pwm_global_idx,0] != 0)[0]
        print example_indices.shape
        
        # collect necessary data, sliced
        pwm_signals = {}
        with h5py.File(motifs_file, "r") as hf:
            pwm_max_vals = hf[max_val_key][example_indices,pwm_global_idx,0]
            pwm_max_indices = hf[max_idx_key][example_indices,pwm_global_idx,0]
            pwm_scores = hf[pwm_scores_key][example_indices,:,:,pwm_global_idx]
            raw_pwm_scores = hf[raw_pwm_scores_key][example_indices,0,:,pwm_global_idx]
            for key in signal_keys:
                pwm_signals[key] = hf[key][example_indices,:]

        # adjust scores
        pwm_scores = np.max(pwm_scores, axis=1) # {N, seqlen}
        score_len = pwm_scores.shape[1]
        
        # set up results arrays
        pwm_aligned_array = np.zeros((example_indices.shape[0], 2*num_positions))
        pwm_aligned_raw_array = np.zeros((example_indices.shape[0], 2*num_positions))

        # and adjust with offset
        for example_idx in range(pwm_aligned_array.shape[0]):
            # get offset and put in with offset
            offset = num_positions - (int(pwm_max_indices[example_idx]) - left_clip)
            pwm_aligned_array[example_idx, offset:(offset+score_len)] = pwm_scores[example_idx]
            pwm_aligned_raw_array[example_idx, offset:(offset+score_len)] = raw_pwm_scores[example_idx]
            
        # and then fix to center
        means = np.mean(pwm_aligned_array, axis=0)
        center_idx = np.argmax(means)
        clip_start = center_idx - final_extend_len
        clip_stop = center_idx + final_extend_len
        
        # clip to completely center the anchor pwm
        pwm_aligned_array = pwm_aligned_array[:,clip_start:clip_stop]
        pwm_aligned_raw_array = pwm_aligned_raw_array[:,clip_start:clip_stop]
        
        # positions
        positions = np.arange(pwm_aligned_array.shape[1]) - int(pwm_aligned_array.shape[1]/2.)
        mid_idx = int(pwm_aligned_array.shape[1]/2.)
        
        # save out aligned array
        aligned_spacing_file = "{}/genome.{}.spacings.txt.gz".format(TMP_DIR, pwm_name_clean)
        pd.DataFrame(data=pwm_aligned_array, columns=positions).to_csv(
            aligned_spacing_file, sep="\t", header=True, index=False, compression="gzip")
        
        # finally get the distributions (save out and plot individuals)
        # NOTE: comparing two slightly different things - for the weighted scores, we want
        # to know if the weights correspond to where they are located to each other - do the scores
        # also integrate in positional info. For the raw pwm scores, we need to consider frequency
        # in the genome, so we do a mean across all examples (include zeros).
        # in most ideal case, the NN highlights positional info that the pwm scores could not get
        spacing_distr_file = "{}/genome.spacing.{}.distribution.txt.gz".format(TMP_DIR, pwm_name_clean)
        active_distr = nonzero_colmeans(pwm_aligned_array)
        active_distr[mid_idx] = 0
        all_distr = nonzero_colmeans(pwm_aligned_raw_array)
        all_distr[mid_idx] = 0
        spacing_distributions = pd.DataFrame({
            "position": positions,
            "active": active_distr / np.sum(active_distr[110:210]),
            "all": all_distr / np.sum(all_distr[110:210])})
        spacing_distributions.to_csv(spacing_distr_file, sep="\t", header=True, index=False, compression="gzip")

        # plot
        # TODO adjust this so doesn't rely on the array shape dim val
        plot_prefix = "{}/genome.spacing.{}".format(OUT_DIR, pwm_name_clean)
        plot_cmd = "{}/plot.spacing.freq.indiv.R {} {} {}".format(
            SCRIPT_DIR, spacing_distr_file, pwm_aligned_array.shape[0], plot_prefix)
        print plot_cmd
        os.system(plot_cmd)

        # plot of motif density by pwm hits (what could you discover by counting in genome)
        spacing_distr_file = "{}/genome.spacing.{}.distribution.genomic_freq.txt.gz".format(
            TMP_DIR, pwm_name_clean)
        all_distr = np.mean(pwm_aligned_raw_array!=0, axis=0)
        all_distr[mid_idx] = 0
        spacing_distributions = pd.DataFrame({
            "position": positions,
            "all": all_distr})
        spacing_distributions.to_csv(spacing_distr_file, sep="\t", header=True, index=False, compression="gzip")
        plot_prefix = "{}/genome.spacing.{}.genomic_frequency".format(OUT_DIR, pwm_name_clean)
        plot_cmd = "{}/plot.spacing.freq.indiv.R {} {} {}".format(
            SCRIPT_DIR, spacing_distr_file, pwm_aligned_array.shape[0], plot_prefix)
        print plot_cmd
        os.system(plot_cmd)

        # TODO - don't bother?
        # plot separated into plots with exact motif counts
        min_region_num = 500
        if False:
            for count in range(2, 6):
                # get examples with exactly n motifs
                exact_count_indices = np.where(np.sum(pwm_aligned_array > 0, axis=1) == count)[0]
                if exact_count_indices.shape[0] < min_region_num:
                    continue
                exact_count_array = pwm_aligned_array[exact_count_indices]
                exact_count_raw_array = pwm_aligned_raw_array[exact_count_indices]

                # get distributions, save, and plot
                spacing_distr_file = "{}/genome.spacing.{}.distribution.exactly_{}.txt.gz".format(
                    TMP_DIR, pwm_name_clean, count)
                active_distr = nonzero_colmeans(exact_count_array)
                active_distr[mid_idx] = 0
                all_distr = nonzero_colmeans(exact_count_raw_array)
                all_distr[mid_idx] = 0
                spacing_distributions = pd.DataFrame({
                    "position": positions,
                    "active": active_distr / np.sum(active_distr),
                    "all": all_distr / np.sum(all_distr)})
                spacing_distributions.to_csv(
                    spacing_distr_file, sep="\t", header=True, index=False, compression="gzip")

                # plot
                plot_prefix = "{}/genome.spacing.{}.exactly_{}".format(MORE_PLOTS_DIR, pwm_name_clean, count)
                plot_cmd = "{}/plot.spacing.freq.indiv.R {} {} {}".format(
                    SCRIPT_DIR, spacing_distr_file, exact_count_array.shape[0], plot_prefix)
                print plot_cmd
                os.system(plot_cmd)
            
        # now go through signals to get N, spacing, task
        # use HITS vs the actual score (unclear how that would multiply with signal)
        pwm_aligned_hits = (pwm_aligned_array > 0).astype(int) # {N, 321}
        for key in signal_keys:
            if "ATAC" in key:
                plot_indices = [0,6,12]
            elif "H3K27ac" in key:
                plot_indices = [0,1,2]
            else:
                plot_indices = []
            
            # do by task
            for task_idx in range(pwm_signals[key].shape[1]):
                if task_idx not in plot_indices:
                    continue
                
                task_signals = pwm_signals[key][:,task_idx] # {N}
                task_signals = np.expand_dims(task_signals, axis=-1)
                task_results = np.multiply(task_signals, pwm_aligned_hits) # {N, 321}
                
                # flatten, remove zeros
                task_df = pd.DataFrame(data=task_results, columns=positions.tolist())
                task_melt_df = task_df.melt()
                task_melt_df = task_melt_df[task_melt_df["value"] != 0]

                # save out
                out_file = "{}/genome.spacing.{}.{}.taskidx-{}.txt.gz".format(
                    TMP_DIR, pwm_name_clean, key, task_idx)
                task_melt_df.to_csv(out_file, sep="\t", header=True, index=False, compression="gzip")

                # plot
                plot_prefix = "{}/genome.spacing.{}.{}.taskidx-{}".format(
                    MORE_PLOTS_DIR, pwm_name_clean, key, task_idx)
                plot_cmd = "Rscript {}/plot.spacing.signals.indiv.R {} {}".format(
                    SCRIPT_DIR, out_file, plot_prefix)
                print plot_cmd
                os.system(plot_cmd)
                
    # ==============================
    # SUMMARIZE
    # ==============================
        
    # now merge all information
    out_file = "{}/genome.spacing.SUMMARY.txt.gz".format(TMP_DIR)
    all_pwms_aligned_array = np.zeros((len(sig_pwms), 2*final_extend_len))
    for pwm_idx in range(len(sig_pwms)):
        pwm_name = sig_pwms[pwm_idx]
        pwm_name_clean = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name_clean = re.sub(".UNK.0.A", "", pwm_name_clean)
        pwm_global_idx = np.where(
            [1 if pwm_name in global_name else 0
             for global_name in all_pwms])[0][0]

        # open spacing distr file and pull info
        spacing_distr_file = "{}/genome.spacing.{}.distribution.txt.gz".format(TMP_DIR, pwm_name_clean)
        data = pd.read_csv(spacing_distr_file, sep="\t")
        all_pwms_aligned_array[pwm_idx] = data["active"].values
        
    # save out
    all_spacings_df = pd.DataFrame(data=all_pwms_aligned_array, columns=data["position"], index=sig_pwms)
    all_spacings_df.to_csv(out_file, sep="\t", header=True, index=True, compression="gzip")

    # and plot
    plot_file = "{}/genome.spacing.SUMMARY.pdf".format(OUT_DIR)
    plot_cmd = "{}/plot.spacing.freq.summary.R {} {}".format(
        SCRIPT_DIR, out_file, plot_file)
    print plot_cmd
    os.system(plot_cmd)

    
    return


main()
