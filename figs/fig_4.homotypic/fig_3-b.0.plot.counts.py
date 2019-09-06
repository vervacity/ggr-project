#!/usr/bin/env python

import os
import sys
import h5py

import numpy as np
import pandas as pd

from tronn.util.utils import DataKeys


def check_nonlinear(y):
    """given 1d array, take first two values to fit line
    """
    # get a linear fit
    slope = y[1] - y[2] # x is 2 - 1 = 1
    intercept = y[1] - slope # x is 1
    
    # then get max val and count
    max_idx = np.argmax(y)
    actual = y[max_idx]
    count_for_max = max_idx + 1
    
    # expected linear val
    expected = slope * count_for_max + intercept

    # if actual > expected, nonlinear
    nonlinear = 0
    if actual > expected:
        #print "nonlinear"
        nonlinear = 1
    
    return nonlinear



def main():
    """aggregate into counts
    """
    # inputs
    motifs_file = sys.argv[1]
    sig_pwms_file = sys.argv[2]
    density_key = DataKeys.ORIG_SEQ_PWM_MAX_DENSITIES
    activity_keys = ["ATAC_SIGNALS.NORM", "H3K27ac_SIGNALS.NORM"]

    manual_keep = ["SMAD3"]

    # make a tmp_dir
    tmp_dir = "tmp"
    os.system("mkdir -p {}".format(tmp_dir))
    
    # params
    POOL_WINDOW = 20 # from scanmotifs
    MAX_COUNT = 6
    plot_timepoint_indices = {
        "ATAC_SIGNALS.NORM": [0,6,12],
        "H3K27ac_SIGNALS.NORM": [0,1,2]} # d0,3,6
    activation_thresh = 0.6 # at what point do you consider it activated
    tol = 1e-5
    
    # load sig pwms
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", index_col=0, header=0).index)

    # load all pwm names and get indices for sig pwms
    with h5py.File(motifs_file, "r") as hf:
        all_pwms = hf[DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_IDX].attrs["pwm_names"]
    pwm_mask = np.array([1 if pwm_name in sig_pwms else 0 for pwm_name in all_pwms])
    sig_indices = np.where(pwm_mask != 0)[0]
    
    # load densities
    with h5py.File(motifs_file, "r") as hf:
        densities = hf[density_key][:]
    counts = densities * POOL_WINDOW

    # load activities
    activities = {}
    with h5py.File(motifs_file, "r") as hf:
        for activity_key in activity_keys:
            activities[activity_key] = hf[
                activity_key][:,plot_timepoint_indices[activity_key]]
            
    # set up results matrices
    heatmap_results = {}
    for activity_key in activity_keys:
        heatmap_results[activity_key] = {}
        heatmap_results[activity_key]["count"] = np.zeros(
            (len(plot_timepoint_indices[activity_key]),
             len(sig_indices),
             MAX_COUNT)) # {timepoint, pwm, count}
        heatmap_results[activity_key]["nonlinear"] = np.zeros(
            (len(sig_indices))) # {pwm}
        
    pwm_results = {}
    for activity_key in activity_keys:
        pwm_results[activity_key] = {} # {N, timepoint+count}
    
    # analyze per pwm and save out
    keep_pwm_names = []
    for i in range(len(sig_indices)):
        sig_idx = sig_indices[i]
        pwm_name = all_pwms[sig_idx]
        keep_pwm_names.append(pwm_name)
        print sig_idx, pwm_name
        
        # extract features
        pwm_counts = counts[:,0,sig_idx]

        # masks:
        # first mask values that aren't (mod 0.05)
        # they're edge values
        # then mask places where there is no density (ie, no motifs)
        keep_mask = np.mod(pwm_counts, 1) < tol
        keep_mask *= pwm_counts > 0
        
        # now collect for each activity key
        for activity_key in activity_keys:
            pwm_activity = activities[activity_key]

            # mask for activity (any), and get remaining indices
            activity_mask = np.any(pwm_activity > 0, axis=-1) * keep_mask
            activity_indices = np.where(activity_mask)[0]

            # convert to pandas
            activity_df = pd.DataFrame(
                data=pwm_activity,
                columns=plot_timepoint_indices[activity_key])
            activity_df["counts"] = pwm_counts.astype(int)

            # filter
            activity_df = activity_df.iloc[activity_indices]

            # save out and plot these results
            pwm_example_counts_file = "{}/{}.{}.counts.per_example.txt.gz".format(tmp_dir, pwm_name, activity_key)
            activity_df.to_csv(pwm_example_counts_file, sep="\t", compression="gzip")
            plot_cmd = "Rscript ~/git/ggr-project/figs/fig_3.homotypic/fig_3-b.0.plot.example_counts.R {}".format(
                pwm_example_counts_file)
            #print plot_cmd
            #os.system(plot_cmd)

            # filter max count
            activity_df = activity_df[activity_df["counts"] <= MAX_COUNT]
            
            # groupby
            activity_per_count_df = activity_df.groupby("counts").median()
            for count in range(1, MAX_COUNT+1):
                if count not in activity_per_count_df.index:
                    activity_per_count_df.loc[count] = 0
            activity_per_count = activity_per_count_df.transpose().values
            heatmap_results[activity_key]["count"][:,i,:] = activity_per_count

    # now only keep those that have strong activation
    for key in heatmap_results.keys():
        threshold = np.percentile(
            heatmap_results[key]["count"].flatten(), 90)
        keep = np.any(heatmap_results[key]["count"] >= threshold, axis=(0,2)).astype(int)

        # check manual
        if True:
            manual = np.zeros_like(keep)
            for i in range(len(keep_pwm_names)):
                for substring in manual_keep:
                    if substring in keep_pwm_names[i]:
                        manual[i] = 1
            keep = np.logical_or(manual, keep).astype(int)
                    
        heatmap_results[key]["nonlinear"][:] = keep

        # save out list of keep pwms
        pwm_names_file = "{}.pwms.keep.txt".format(key)
        pd.DataFrame(
            np.array(keep_pwm_names)[keep == 1]).to_csv(pwm_names_file, header=False, index=False)
        
    # with all results collected, save out and plot
    h5_results_file = "motifs.counts_v_activity.h5"
    if not os.path.isfile(h5_results_file):
        for key in heatmap_results.keys():
            print key
            with h5py.File(h5_results_file, "a") as out:
                for sub_key in heatmap_results[key].keys():
                    save_key = "{}/{}".format(key, sub_key)
                    out.create_dataset(save_key, data=heatmap_results[key][sub_key])
                    out[save_key].attrs["pwm_names"] = keep_pwm_names

                
    # plot all
    plot_cmd = "Rscript ~/git/ggr-project/figs/fig_3.homotypic/plot_by_density.R {} counts.ALL FALSE {}".format(
        h5_results_file, " ".join(activity_keys))
    print plot_cmd
    os.system(plot_cmd)

    # plot some
    plot_cmd = "Rscript ~/git/ggr-project/figs/fig_3.homotypic/plot_by_density.R {} counts.filt TRUE {}".format(
        h5_results_file, " ".join(activity_keys))
    print plot_cmd
    os.system(plot_cmd)

    
    return


main()
