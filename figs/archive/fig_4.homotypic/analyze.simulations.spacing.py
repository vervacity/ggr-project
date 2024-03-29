#!/usr/bin/env python

import os
import re
import sys
import h5py
import glob

import numpy as np
import pandas as pd

def main():
    """analyze density sims
    """
    # ==============================
    # SETUP
    # ==============================
    
    # inputs
    SCRIPT_DIR = sys.argv[1]
    SIM_DIR = sys.argv[2]
    OUT_DIR = sys.argv[3]
    prefix = "simulations.spacing"
    
    # set up dirs
    os.system("mkdir -p {}".format(OUT_DIR))
    TMP_DIR = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(TMP_DIR))
    MORE_PLOTS_DIR = "{}/plots".format(OUT_DIR)
    os.system("mkdir -p {}".format(MORE_PLOTS_DIR))
    
    # keys
    predictions_key = "logits.norm"
    sample_idx_key = "simul.pwm.sample_idx"
    pwm_dist_key = "simul.pwm.dist"
    string_key = "grammar.string"

    # params
    plot_indices = [0,6,12,13,14,15]

    # ==============================
    # ANALYZE
    # ==============================
    
    # analyze each pwm
    results_files = sorted(
        glob.glob("{}/*/ggr.predictions.h5".format(SIM_DIR)))
    all_results = None
    for results_file in results_files:
        
        # pwm name
        orig_pwm_name = results_file.split("/")[-2]
        pwm_name = re.sub("HCLUST-\\d+_", "", orig_pwm_name)
        pwm_name = re.sub(".UNK.0.A", "", pwm_name)
        print pwm_name

        # TODO analyze directionality? <- but then will need to find these in genome...
        # for now, only analyze FWD pwm
        keep_grammar_string = ",".join(
            ["{}+".format(orig_pwm_name),
             "{}+".format(orig_pwm_name)])
        
        # read in data
        with h5py.File(results_file, "r") as hf:
            sample_ids = hf[sample_idx_key][:]
            pwm_dist = hf[pwm_dist_key][:]
            predictions = hf[predictions_key][:]
            grammar_string = hf[string_key][:,0]

        # get sample ids
        unique_sample_ids = sorted(list(set(sample_ids.tolist())))
        
        # divide by task and analyze per task
        num_tasks = predictions.shape[1]
        for task_idx in range(num_tasks):

            # results_file
            task_prefix = "{}.{}.taskidx-{}".format(
                prefix, pwm_name, task_idx)
            results_file = "{}/{}.aggregated.txt.gz".format(TMP_DIR, task_prefix)

            # run analysis if results file doesn't exist
            if not os.path.isfile(results_file):

                # get task data
                task_data = pd.DataFrame({
                    "sample_id": sample_ids,
                    "pwm_dist": pwm_dist,
                    "prediction": predictions[:,task_idx],
                    "grammar_string": grammar_string})
                task_data_norm = None

                # for each sample, subtract background
                for sample_idx in range(len(unique_sample_ids)):
                    sample_id = unique_sample_ids[sample_idx]
                    sample_data = task_data[task_data["sample_id"] == sample_id]

                    # confirm that background exists
                    if sample_data[sample_data["grammar_string"] == "BACKGROUND"].shape[0] == 0:
                        continue

                    # confirm that there is at least 1 other sample
                    if sample_data[sample_data["grammar_string"] != "BACKGROUND"].shape[0] == 0:
                        continue

                    # delete background
                    background_val = sample_data["prediction"][sample_data["grammar_string"] == "BACKGROUND"].values[0]
                    sample_data = sample_data[sample_data["grammar_string"] != "BACKGROUND"]
                    sample_data["prediction"] = sample_data["prediction"] - background_val

                    # if the sequence got worse (negative after delete), don't save out
                    if np.mean(sample_data["prediction"].values) <= 0:
                        continue

                    # for now, only keep FWD pwm setup
                    if True:
                        sample_data = sample_data[sample_data["grammar_string"] == keep_grammar_string]
                        sample_data = sample_data.drop(["grammar_string"], axis=1)

                    # save out
                    if sample_idx == 0:
                        task_data_norm = sample_data
                    else:
                        task_data_norm = pd.concat([task_data_norm, sample_data], axis=0)

                # and save out to plot in R
                if task_data_norm is None:
                    continue
                task_data_norm.to_csv(results_file, header=True, index=False, sep="\t", compression="gzip")

            # plot
            if task_idx in plot_indices:
                plot_file = "{}/{}.pdf".format(MORE_PLOTS_DIR, task_prefix)
                plot_cmd = "{}/plot.spacing.simulations.indiv.R {} {}".format(
                    SCRIPT_DIR, results_file, plot_file)
                print plot_cmd
                os.system(plot_cmd)

    # ==============================
    # SUMMARIZE
    # ==============================
    # build the matrix by considering: fraction of samples that activated, also
    # did the average across the samples activate enough
    h5_results_file = "{}/simulations.spacing_v_activity.h5".format(OUT_DIR)
    keys = ["ATAC", "H3K27ac"]
    task_indices = [0, 6, 12]
    averages = []
    sample_fracts = []
    for task_idx in task_indices:
        results_files = sorted(
            glob.glob("{}/*taskidx-{}*".format(TMP_DIR, task_idx)))
        task_averages = []
        task_sample_fracts = []
        pwm_names = []
        for results_file in results_files:
            pwm_name = os.path.basename(results_file).split(".")[2]
            data = pd.read_csv(results_file, sep="\t")
            # get average results
            data_tmp = data[["prediction", "pwm_dist"]]
            data_tmp = data_tmp.groupby(["pwm_dist"]).median().transpose()
            data_tmp.index = [pwm_name]
            task_averages.append(data_tmp)

            # also track how many samples
            sample_num = set(data["sample_id"].values.tolist())
            task_sample_fracts.append(len(sample_num))
            pwm_names.append(pwm_name)

        # save out averages
        # note that groupby has sort=True by default, so should be ordered
        # so can concatenate in list without needing to merge on name
        task_averages = pd.concat(task_averages, axis=0)
        averages.append(task_averages.values)

        # save out counts
        task_sample_fracts = pd.DataFrame(
            {"taskidx-{}".format(task_idx): task_sample_fracts}, index=pwm_names)
        sample_fracts.append(task_sample_fracts.values)

    # consolidate and save out
    averages = np.stack(averages, axis=0)
    sample_fracts = np.squeeze(np.stack(sample_fracts, axis=0), axis=-1)
    
    with h5py.File(h5_results_file, "a") as out:
        key = "ATAC/count"
        out.create_dataset(key, data=averages)
        out[key].attrs["pwm_names"] = pwm_names
        key = "ATAC/sample_count"
        out.create_dataset(key, data=sample_fracts)
        out[key].attrs["pwm_names"] = pwm_names
    
    
    # H3K27ac: tasks 13,14,15
    task_indices = [13, 14, 15]
    averages = []
    sample_fracts = []
    for task_idx in task_indices:
        results_files = sorted(
            glob.glob("{}/*taskidx-{}*".format(TMP_DIR, task_idx)))
        task_averages = []
        task_sample_fracts = []
        pwm_names = []
        for results_file in results_files:
            pwm_name = os.path.basename(results_file).split(".")[2]
            data = pd.read_csv(results_file, sep="\t")
            # get average results
            data_tmp = data[["prediction", "pwm_dist"]]
            data_tmp = data_tmp.groupby(["pwm_dist"]).median().transpose()
            data_tmp.index = [pwm_name]
            task_averages.append(data_tmp)

            # also track how many samples
            sample_num = set(data["sample_id"].values.tolist())
            task_sample_fracts.append(len(sample_num))
            pwm_names.append(pwm_name)

        # save out averages
        # TODO potentially need to do a merge, not an append
        task_averages = pd.concat(task_averages, axis=0)
        averages.append(task_averages.values)

        # save out counts
        task_sample_fracts = pd.DataFrame(
            {"taskidx-{}".format(task_idx): task_sample_fracts}, index=pwm_names)
        sample_fracts.append(task_sample_fracts.values)

    # consolidate and save out
    averages = np.stack(averages, axis=0)
    sample_fracts = np.squeeze(np.stack(sample_fracts, axis=0), axis=-1)

    
    with h5py.File(h5_results_file, "a") as out:
        key = "H3K27ac/count"
        out.create_dataset(key, data=averages)
        out[key].attrs["pwm_names"] = pwm_names
        key = "H3K27ac/sample_count"
        out.create_dataset(key, data=sample_fracts)
        out[key].attrs["pwm_names"] = pwm_names
        
    
    # plot
    plot_cmd = "{}/plot.spacing.simulations.summary.R {} {}/simulations.spacing FALSE {}".format(
        SCRIPT_DIR, h5_results_file, OUT_DIR, "ATAC H3K27ac")
    print plot_cmd
    os.system(plot_cmd)

    return

main()
