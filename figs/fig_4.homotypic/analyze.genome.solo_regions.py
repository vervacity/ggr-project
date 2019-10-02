#!/usr/bin/env python

import os
import re
import sys
import h5py

import numpy as np
import pandas as pd

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import regions_to_genes

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


def load_data_from_multiple_h5_files(h5_files, key, example_indices=None):
    """convenience wrapper
    example indices is a list of arrays
    """
    key_data = []
    for h5_idx in range(h5_files):
        h5_file = h5_files[h5_idx]
        h5_indices = example_indices[h5_idx]
        with h5py.File(h5_file, "r") as hf:
            if example_indices is not None:
                key_data.append(hf[key][example_indices])
            else:
                key_data.append(hf[key][:])
    key_data = np.concatenate(key_data, axis=0)
    
    return key_data


def main():
    """analyze spacing
    """
    # inputs
    SCRIPT_DIR = "/users/dskim89/git/ggr-project/figs/fig_4.homotypic"
    OUT_DIR = sys.argv[1]
    sig_pwms_file = sys.argv[2]
    links_file = sys.argv[3]
    tss_file = sys.argv[4]
    background_gene_set_file = sys.argv[5]
    motifs_files = sys.argv[6:]

    # dirs
    tmp_dir = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(tmp_dir))
    more_plots_dir = "{}/plots".format(OUT_DIR)
    os.system("mkdir -p {}".format(more_plots_dir))
    
    # keys
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    #max_idx_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_IDX
    pwm_scores_key = DataKeys.WEIGHTED_SEQ_PWM_SCORES_THRESH
    #raw_pwm_scores_key = DataKeys.ORIG_SEQ_PWM_SCORES_THRESH
    importances_key = DataKeys.WEIGHTED_SEQ_ACTIVE
    metadata_key = DataKeys.SEQ_METADATA
    signal_keys = [
        "ATAC_SIGNALS",
        "H3K27ac_SIGNALS"]
    
    # params
    left_clip = 420 + 12
    num_positions = 200
    final_extend_len = 160
    MIN_REGION_COUNT = 100
    impt_clip_start = 12
    impt_clip_end = 148
    pwm_window = 20
    fract_thresh = 0.9

    
    # read in the sig pwms
    print "WARNING - check pwms file"
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", header=0).iloc[:,0])

    # read in the list of all pwm names
    with h5py.File(motifs_files[0], "r") as hf:
        all_pwms = hf[pwm_scores_key].attrs["pwm_names"]

    # get max vals (use to filter for which sequences have pwm)
    max_vals_sets = []
    for motifs_file in motifs_files:
        with h5py.File(motifs_file, "r") as hf:
            max_vals_sets.append(hf[max_val_key][:]) # {N, pwm, 1}
    
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

        # skip
        #continue
        if "GATA3" not in pwm_name_clean:
            continue
        
        # check to see which have this pwm and adjust indices
        example_indices = []
        for max_vals in max_vals_sets:
            example_indices.append(np.where(max_vals[:,pwm_global_idx,0] != 0)[0])

        # debug
        if False:
            total = 0
            for example_index_set in example_indices:
                total += example_index_set.shape[0]
            print total
        
        # collect necessary data, sliced
        pwm_scores = []
        metadata = []
        importances = [] # TODO
        for motifs_file_idx in range(len(motifs_files)):
            print motifs_files[motifs_file_idx]
            with h5py.File(motifs_files[motifs_file_idx], "r") as hf:
                pwm_scores.append(
                    hf[pwm_scores_key][example_indices[motifs_file_idx], :, :, pwm_global_idx])
                metadata.append(
                    hf[metadata_key][example_indices[motifs_file_idx], 0])
                importances_reduced = np.sum(
                    hf[importances_key][example_indices[motifs_file_idx],:,:,:], axis=(1,3))
                importances_reduced = importances_reduced[:,impt_clip_start:impt_clip_end]
                importances.append(importances_reduced)
        pwm_scores = np.concatenate(pwm_scores, axis=0)
        metadata = np.concatenate(metadata, axis=0)
        importances = np.concatenate(importances, axis=0)
        
        # adjust scores
        pwm_scores = np.max(pwm_scores, axis=1) # {N, seqlen}
        score_len = pwm_scores.shape[1]
        
        # now basically want to calculate a fraction of scores that fall in pwm window
        pwm_impt_coverage = np.zeros((pwm_scores.shape[0]))
        pwm_hit_count = np.zeros((pwm_scores.shape[0]))
        for example_idx in range(pwm_scores.shape[0]):
            # figure out how many positions were marked
            impt_sum_total = np.sum(importances[example_idx] > 0)
            
            # per index, get window around and collect importance scores
            pwm_pos_indices = np.where(pwm_scores[example_idx] > 0)[0]
            pwm_hit_count[example_idx] = pwm_pos_indices.shape[0]
            
            # collect importance scores within window
            impt_sum_pwm_overlap = 0
            for pwm_pos_idx in pwm_pos_indices:
                start = max(pwm_pos_idx - pwm_window/2, 0)
                end = min(pwm_pos_idx + pwm_window/2, score_len)
                impt_sum_pwm_overlap += np.sum(importances[example_idx,start:end] > 0)

            # and check
            fract_covered = impt_sum_pwm_overlap / float(impt_sum_total)
            pwm_impt_coverage[example_idx] = fract_covered

        # filter
        keep_indices = np.where(pwm_impt_coverage >= fract_thresh)[0]
        metadata = metadata[keep_indices]
        pwm_hit_count = pwm_hit_count[keep_indices]
        
        # now use fract to choose a cutoff and only keep those above cutoff
        import ipdb
        ipdb.set_trace()

        # look at different count levels
        count_thresholds = range(1, 6)
        for count_threshold in count_thresholds:

            # get matching regions
            thresholded_indices = np.where(pwm_hit_count >= count_threshold)[0]
            thresholded_metadata = metadata[thresholded_indices]
            if thresholded_metadata.shape[0] < MIN_REGION_COUNT:
                continue

            # make bed
            bed_file = "{}/{}.count-{}.bed.gz".format(OUT_DIR, pwm_name_clean, count_threshold)
            array_to_bed(
                thresholded_metadata,
                bed_file, interval_key="active", merge=True)

            # get gene set enrichment
            gene_set_file = "{}.gene_set.txt.gz".format(bed_file.split(".bed")[0])
            regions_to_genes(
                bed_file,
                links_file,
                tss_file,
                gene_set_file,
                filter_by_score=0.5)

            # gprofiler
            gprofiler_dir = "{}/enrichments.{}".format(OUT_DIR, pwm_name_clean)
            run_shell_cmd("mkdir -p {}".format(gprofiler_dir))
            if run_enrichments:
                run_gprofiler(
                    gene_set_file, background_gene_file,
                    gprofiler_dir, ordered=True, header=True)
            

        quit()
        summary_df = None
        for position_range in position_ranges:
            print ">> range:", position_range
            prefix = "{}/{}.range_{}-{}".format(
                tmp_dir, pwm_name_clean, position_range[0], position_range[1]-1)
            position_range = np.array(position_range)
            
            # check positive side
            pos_position_range = position_range + mid_idx
            #print pos_position_range
            positive_present = np.sum(
                pwm_aligned_array[:,pos_position_range[0]:pos_position_range[1]] > 0, axis=1) > 0

            # check negative side
            neg_position_range = -np.flip(position_range) + mid_idx
            #print neg_position_range
            negative_present = np.sum(
                pwm_aligned_array[:,neg_position_range[0]:neg_position_range[1]] > 0, axis=1) > 0

            # combine
            distance_present = np.logical_or(positive_present, negative_present)
            distance_indices = np.where(distance_present)[0]
            print "num regions:", distance_indices.shape[0]

            pwm_counts = np.sum(pwm_aligned_array[distance_indices] > 0, axis=1)
            pwm_count_indices = np.where(pwm_counts == 2)[0]
            print pwm_count_indices.shape
            
            
            if distance_indices.shape[0] > MIN_REGION_COUNT:
                #thresholded_metadata = metadata[distance_indices]
                thresholded_metadata = metadata[distance_indices]#[pwm_count_indices]
                # build the BED file
                bed_file = "{}.bed.gz".format(prefix)
                array_to_bed(
                    thresholded_metadata,
                    bed_file, interval_key="active", merge=True)
                
                # track ATAC/H3K27ac?
                
                # then match to proximal gene set
                # TODO - think about adjusting max dist OR use distance based links
                gene_set_file = "{}.gene_set.txt.gz".format(bed_file.split(".bed")[0])
                regions_to_genes_through_links(
                    bed_file,
                    links_file,
                    gene_set_file)

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
                thresholded_summary["range"] = "-".join(
                    [str(val) for val in position_range.tolist()])
                print "term count:", thresholded_summary.shape[0]

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

        summary_df["log10pval"] = -np.log10(summary_df["p.value"].values)
        summary_df = summary_df.drop("p.value", axis=1)
        summary_file = "{}/{}.summary.txt.gz".format(tmp_dir, pwm_name_clean)
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
        summary_df = summary_df.sort_values(["term.name", "range"])
        summary_file = "{}/{}.summary.filt.txt.gz".format(tmp_dir, pwm_name_clean)
        summary_df.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")
        
        quit()

    return


main()
