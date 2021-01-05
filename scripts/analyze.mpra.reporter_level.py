
import os
import h5py

import numpy as np
import pandas as pd

from tronn.plot.visualization import scale_scores
from tronn.plot.visualization import plot_weights_group


def _check_avg_traj(example_data, rule_id, amp_thresh=0.0):
    """filter for those that actually match pattern
    """
    days = ["d0", "d3", "d6"]
    reps = ["b1", "b2", "b3", "b4"]

    val_headers = []
    for day in days:
        for rep in reps:
            val_headers.append("{}_{}".format(day, rep))

    # analyze first by rep (focus on traj)
    example_results = {}
    for rep in reps:
        rep_headers = [header for header in example_data[val_headers].columns
                       if rep in header]
        rep_data = example_data[rep_headers]
        
        # ignore missing barcodes in rep
        rep_data = rep_data[rep_data.sum(axis=1) != 0]

        # average
        rep_traj = rep_data.mean(axis=0)
        example_results[rep] = rep_traj.values

        
    # summary
    example_results = pd.DataFrame(example_results)

    # ignore nan
    if example_results.isnull().values.any():
        return False

    # check if rule is early or mid/late
    early_rules = ["TRAJ_LABELS-0", "TRAJ_LABELS-8-10-11", "TRAJ_LABELS-9"]
    is_early_rule = False
    for early_rule in early_rules:
        if early_rule in rule_id:
            is_early_rule = True

    # then check if max time val matches
    max_timepoint_idx = np.argmax(
        np.mean(example_results.values, axis=1))
    if is_early_rule:
        result = (max_timepoint_idx == 0)
    else:
        result = (max_timepoint_idx != 0)

    # and check that max val is higher than amp_thresh
    #result = (np.max(example_results.values) > amp_thresh) and result
        
    #if result:
        #print np.max(example_results.values)
        #print np.mean(example_results.values, axis=1)
    
    return result



def main():
    """reporter level analysis of elements in the MPRA
    """
    
    # inputs
    MPRA_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results"
    mpra_data_file = "{}/results/ggr.signal.w_metadata.mat.txt.gz".format(MPRA_DIR)

    # rules of interest
    rules = [
        "ggr.TRAJ_LABELS-0_x_TRAJ_LABELS-8-10-11.50_x_26", # ATF1/CREB1 x LEF1
        "ggr.TRAJ_LABELS-2.grammar-99.annot-316", # TP63 x ATF1
        "ggr.TRAJ_LABELS-2.grammar-59.annot-272", # CEBPA x TFAP2
        "ggr.TRAJ_LABELS-2.grammar-60.annot-274", # CEBPD x TFAP2
        "ggr.TRAJ_LABELS-1_x_TRAJ_LABELS-14_x_TRAJ_LABELS-2.42_x_21_x_58", # KLF x ZNF
        "ggr.TRAJ_LABELS-0_x_TRAJ_LABELS-9.16_x_125" # HOXA1 x ERG
    ]
    
    # load data
    mpra_data = pd.read_csv(mpra_data_file, sep="\t")

    # set up rule id and only keep GGR
    mpra_data["rule"] = mpra_data["example_id"].str.replace("\.\d+$", "", regex=True)
    mpra_data = mpra_data[mpra_data["rule"].str.contains("ggr")]

    # reduce down to only those of interest
    mpra_data = mpra_data[mpra_data["rule"].isin(rules)]

    # value headers
    days = ["d0", "d3", "d6"]
    reps = ["b1", "b2", "b3", "b4"]
    val_headers = []
    for day in days:
        for rep in reps:
            val_headers.append("{}_{}".format(day, rep))
    
    # for each rule, filter for best representative element
    # mostly here, filter for trajectory and then manually look later
    # NOTE: this is about 19 x 6 ~ 114 plots to look through
    for rule_id in rules:
        print rule_id
        rule_data = mpra_data[mpra_data["rule"] == rule_id]
        
        # remove zeros
        rule_data = rule_data[
            rule_data[val_headers].sum(axis=1) != 0]

        # analyze at individual example level
        example_ids = list(set(rule_data["example_id"].values))
        good_example_total = 0
        good_example_ids = []
        for example_id in example_ids:
            
            example_data = rule_data[rule_data["example_id"] == example_id]

            # only care about endog version for checks
            example_data = example_data[
                example_data["combos"] == "0,0"]
            
            # filter if not good at a summary level
            is_correct_avg_traj = _check_avg_traj(example_data, rule_id, amp_thresh=1.0)
            if not is_correct_avg_traj:
                continue
            
            good_example_ids.append(example_id)
            good_example_total += 1
            
        print good_example_total

        # plot all out on one sheet
        good_rule_data = rule_data[rule_data["example_id"].isin(good_example_ids)]
        good_rule_data = good_rule_data[
            good_rule_data["combos"] == "0,0"]

        # extract and save to text file
        rule_data_file = "{}.endog_data.txt.gz".format(rule_id)
        good_rule_data.to_csv(rule_data_file, sep="\t", compression="gzip")

        # plot in R
        # plot each rep separately
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.rule_indiv_example.R {}".format(rule_data_file)
        print plot_cmd
        #os.system(plot_cmd)

        #print good_rule_data.columns
        #quit()

    # manual checks:
    # MPRA activity pattern - most important
    # Epigenetic pattern (ATAC, H3K27ac) - this is mostly guaranteed by filtering, most important to check on genome browser
    # NN prediction
    # Most consistent across barcodes
    # High amplitude
    # Genome browser
    # any available ChIP-seq


    # selections - no more than 15
    selections = {
        "ggr.TRAJ_LABELS-0_x_TRAJ_LABELS-8-10-11.50_x_26": [163, 1023, 1249], # ATF1/CREB1 x LEF1
        "ggr.TRAJ_LABELS-0_x_TRAJ_LABELS-9.16_x_125": [33, 945, 1102], # HOXA1 x ERG
        "ggr.TRAJ_LABELS-1_x_TRAJ_LABELS-14_x_TRAJ_LABELS-2.42_x_21_x_58": [9, 13, 132, 187], # KLF x ZNF
        "ggr.TRAJ_LABELS-2.grammar-59.annot-272": [], # CEBPA x TFAP2
        "ggr.TRAJ_LABELS-2.grammar-60.annot-274": [278, 343, 508, 669], # CEBPD x TFAP2
        "ggr.TRAJ_LABELS-2.grammar-99.annot-316": [92, 314] # TP63 x ATF1
    }

    # set up keep ids, and plot out importance scores
    interactions_dir = "/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/dmim.shuffle.complete/synergy"
    keep_ids = []
    for rule_id in selections.keys():
        for example_fileidx in selections[rule_id]:
            keep_id = "{}.{}".format(rule_id, example_fileidx)
            keep_ids.append(keep_id)

            # also print out importance scores - across time and mutated
            synergy_file = "{}/{}/ggr.synergy.h5".format(
                interactions_dir, rule_id)
            with h5py.File(synergy_file, "r") as hf:
                #for key in sorted(hf.keys()): print key, hf[key].shape

                # across time
                orig_importances = hf["sequence-weighted"][example_fileidx][:,420:580]
                match_importances = hf["sequence-weighted.active"][example_fileidx]
                importances = scale_scores(orig_importances, match_importances)
                plot_file = "{}.plot.pdf".format(keep_id)
                print plot_file
                plot_weights_group(importances, plot_file)

                # TODO mutated

            # TODO also confirm correct region
                
    # filter for selected examples
    selected_data = mpra_data[mpra_data["example_id"].isin(keep_ids)]
    
    # remove zeros
    selected_data = selected_data[
        selected_data[val_headers].sum(axis=1) != 0]
    
    # only care about endog version for checks
    selected_data = selected_data[
        selected_data["combos"] == "0,0"]
    print len(list(set(selected_data["example_id"].values)))
    
    #print selected_data
    
    
    

    return


main()
