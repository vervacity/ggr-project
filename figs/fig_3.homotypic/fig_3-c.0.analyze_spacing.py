#!/usr/bin/env python

import os
import re
import sys
import h5py

import numpy as np
import pandas as pd

from tronn.util.utils import DataKeys


def main():
    """analyze spacing
    """
    # dirs
    tmp_dir = "tmp"
    os.system("mkdir -p {}".format(tmp_dir))
    more_plots_dir = "plots"
    os.system("mkdir -p {}".format(more_plots_dir))
    SCRIPT_DIR = "/users/dskim89/git/ggr-project/figs/fig_3.homotypic"
    
    # inputs
    motifs_file = sys.argv[1]
    sig_pwms_file = sys.argv[2]

    # keys
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL #"sequence-weighted.active.pwm-scores.thresh.max.val"
    max_idx_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_IDX #"sequence-weighted.active.pwm-scores.thresh.max.idx"
    pwm_scores_key = DataKeys.WEIGHTED_SEQ_PWM_SCORES_THRESH #"sequence-weighted.active.pwm-scores.thresh"
    raw_pwm_scores_key = DataKeys.ORIG_SEQ_PWM_SCORES_THRESH #"sequence.active.pwm-scores.thresh"
    signal_keys = [
        "ATAC_SIGNALS",
        "H3K27ac_SIGNALS"]
    
    # params
    left_clip = 420 + 12
    num_positions = 200
    final_extend_len = 160
    out_prefix = "fig_5-c"

    # read in the sig pwms
    sig_pwms = list(pd.read_csv(sig_pwms_file, header=None).iloc[:,0])

    # read in the list of all pwm names
    with h5py.File(motifs_file, "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    
    # read in max vals
    with h5py.File(motifs_file, "r") as hf:
        max_vals = hf[max_val_key][:] # {N, pwm, 1}
    
    # analyze each pwm
    for pwm_idx in range(len(sig_pwms)):

        # name and global index
        pwm_name = sig_pwms[pwm_idx]
        pwm_name_clean = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name_clean = re.sub(".UNK.0.A", "", pwm_name_clean)
        pwm_global_idx = np.where(
            [1 if pwm_name in global_name else 0
             for global_name in all_pwms])[0][0]
        print pwm_name_clean, pwm_global_idx

        #if "TP53" not in pwm_name_clean:
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
        
        # finally get the distributions (save out and plot individuals)
        spacing_distr_file = "{}/{}.{}.spacings.txt.gz".format(tmp_dir, out_prefix, pwm_name_clean)
        spacing_distributions = pd.DataFrame({
            "position": positions,
            "active": np.mean(pwm_aligned_array, axis=0),
            "all": np.mean(pwm_aligned_raw_array, axis=0)})
        # TODO some kind of normalization here
        spacing_distributions.to_csv(spacing_distr_file, sep="\t", header=True, index=False, compression="gzip")

        # plot
        plot_prefix = "{}.{}".format(out_prefix, pwm_name_clean)
        plot_cmd = "Rscript {}/fig_3-c.1.plot.spacing.indiv.R {} {}".format(
            SCRIPT_DIR, spacing_distr_file, plot_prefix)
        print plot_cmd
        os.system(plot_cmd)
        
        # now go through signals to get N, spacing, task
        # use HITS vs the actual score (unclear how that would multiply with signal)
        pwm_aligned_hits = (pwm_aligned_array > 0).astype(int) # {N, 321}
        for key in signal_keys:
            # do by task
            for task_idx in range(pwm_signals[key].shape[1]):
                task_signals = pwm_signals[key][:,task_idx] # {N}
                task_signals = np.expand_dims(task_signals, axis=-1)
                task_results = np.multiply(task_signals, pwm_aligned_hits) # {N, 321}
                
                # flatten, remove zeros
                task_df = pd.DataFrame(data=task_results, columns=positions.tolist())
                task_melt_df = task_df.melt()
                task_melt_df = task_melt_df[task_melt_df["value"] != 0]

                # save out
                out_file = "{}/{}.{}.{}.taskidx-{}.txt.gz".format(
                    tmp_dir, out_prefix, pwm_name_clean, key, task_idx)
                task_melt_df.to_csv(out_file, sep="\t", header=True, index=False, compression="gzip")

                # plot
                plot_prefix = "{}/{}.{}.{}.taskidx-{}".format(
                    more_plots_dir, out_prefix, pwm_name_clean, key, task_idx)
                plot_cmd = "Rscript {}/fig_3-c.2.plot.signals.indiv.R {} {}".format(
                    SCRIPT_DIR, out_file, plot_prefix)
                print plot_cmd
                os.system(plot_cmd)

        # TODO try the fourier here?
        # save out aligned array, clip as needed?
        
        
    # now merge all information
    out_file = "{}.ALL.spacings.txt.gz".format(out_prefix)
    all_pwms_aligned_array = np.zeros((len(sig_pwms), 2*final_extend_len))
    for pwm_idx in range(len(sig_pwms)):

        # name and global index
        pwm_name = sig_pwms[pwm_idx]
        pwm_name_clean = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name_clean = re.sub(".UNK.0.A", "", pwm_name_clean)
        pwm_global_idx = np.where(
            [1 if pwm_name in global_name else 0
             for global_name in all_pwms])[0][0]
        print pwm_name_clean, pwm_global_idx

        # open spacing distr file
        spacing_distr_file = "{}/{}.{}.spacings.txt.gz".format(tmp_dir, out_prefix, pwm_name_clean)
        data = pd.read_csv(spacing_distr_file, sep="\t")
        print data

        # pull info
        all_pwms_aligned_array[pwm_idx] = data["active"].values
        
    # save out
    all_spacings_df = pd.DataFrame(data=all_pwms_aligned_array, columns=positions, index=sig_pwms)
    all_spacings_df.to_csv(out_file, sep="\t", header=True, index=True, compression="gzip")
    
    
    return


main()
