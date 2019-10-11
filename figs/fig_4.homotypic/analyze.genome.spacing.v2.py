#!/usr/bin/env python

import os
import re
import sys
import h5py

import numpy as np
import pandas as pd

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import regions_to_genes

from tronn.interpretation.syntax import analyze_syntax
from tronn.interpretation.syntax import recombine_syntax_results

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


def _get_gene_set_enrichment(metadata, links_file, tss_file, background_gene_set_file, prefix, tmp_dir):
    """
    """
    # build the BED file
    bed_file = "{}.bed.gz".format(prefix)
    array_to_bed(
        metadata,
        bed_file, interval_key="active", merge=True)
    
    # then match to proximal gene set
    gene_set_file = "{}.gene_set.txt.gz".format(bed_file.split(".bed")[0])
    regions_to_genes(
        bed_file,
        links_file,
        tss_file,
        gene_set_file,
        filter_by_score=0.5,
        filter_genes=[]) # just dynamic?

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
    thresholded_summary = thresholded_summary[
        ["term.id", "p.value", "term.name"]]
    print "term count:", thresholded_summary.shape[0]

    return thresholded_summary


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

    # signals
    signal_keys = ["ATAC_SIGNALS", "H3K27ac_SIGNALS"]
    
    # read in the sig pwms
    print "WARNING - check pwms file"
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", header=0).iloc[:,0])
    
    # read in the list of all pwm names
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    with h5py.File(motifs_files[0], "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    num_pwms = len(all_pwms)
    print num_pwms
        
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
        pwm_indices = [pwm_global_idx, pwm_global_idx]

        # reverse is in second half of indices
        rc_idx = pwm_global_idx + num_pwms

        if "RELA" not in pwm_name_clean:
            continue

        # results dict
        results = {}

        # go through orientations
        possible_orientations = ["FWD", "REV"]
        for orientation_a in possible_orientations:
            for orientation_b in possible_orientations:
                # set up new pwm indices
                pwm_indices = []
                
                # figure out whether fwd or rev
                if orientation_a == "FWD":
                    pwm_indices.append(pwm_global_idx)
                elif orientation_a == "REV":
                    pwm_indices.append(rc_idx)
                
                # figure out whether fwd or rev
                if orientation_b == "FWD":
                    pwm_indices.append(pwm_global_idx)
                elif orientation_b == "REV":
                    pwm_indices.append(rc_idx)
                    
                print pwm_indices
                assert len(pwm_indices) == 2
                results_key = "{}_{}".format(orientation_a, orientation_b)
                print results_key
                
                # build aligned array around first pwm
                aligned_results = analyze_syntax(
                    motifs_files, pwm_indices) # TODO add in solo filter here?
                
                if aligned_results is None:
                    continue

                if orientation_b == "FWD":
                    results[results_key] = aligned_results[pwm_global_idx]
                elif orientation_b == "REV":
                    results[results_key] = aligned_results[rc_idx]

        # 1) AGNOSTIC:  both(anchor)_x_both - merge all
        # 2) AGNOS_IN: both(anchor)_x_in - fwd_REV(+) + rev_REV(+) + fwd_FWD(-) + rev_FWD(-)
        # 3) AGNOS_OUT: both(anchor)_x_out - fwd_FWD(+) + rev_FWD(+) + fwd_REV(-) + rev_REV(-)
        # 4) SAME_DIR: same_x_same - fwd_FWD + rev_REV
        # 5) IN: in_x_in - fwd_REV(+) + rev_FWD(-)
        # 6) OUT: out_x_out - rev_FWD(+) + fwd_REV(-)
        # NOTE: ignoring rev_BOTH and fwd_BOTH - same as 2,3

        adjustments = {
            "BOTH_BOTH": (
                ["FWD_FWD", "FWD_REV", "REV_FWD", "REV_REV"],
                [".", ".", ".", "."]),
            "BOTH_IN": (
                ["FWD_FWD", "FWD_REV", "REV_FWD", "REV_REV"],
                ["-", "+", "-", "+"]),
            "BOTH_OUT": (
                ["FWD_FWD", "FWD_REV", "REV_FWD", "REV_REV"],
                ["+", "-", "+", "-"]),
            "SAME": (
                ["FWD_FWD", "REV_REV"],
                [".", "."]),
            "IN": (
                ["FWD_REV", "REV_FWD"],
                ["+", "-"]),
            "OUT": (
                ["FWD_REV", "REV_FWD"],
                ["-", "+"])}
        
        results_adjusted = {}
        enrichments_summary = None
        for adj_key in sorted(adjustments.keys()):
            print adj_key
            results_adjusted[adj_key] = recombine_syntax_results(
                results,
                adjustments[adj_key][0],
                adjustments[adj_key][1],
                signal_keys)
            
            # plot results
            plot_prefix = "{}/{}.{}".format(
                OUT_DIR, pwm_name_clean, adj_key)

            # plot pwm scores
            score_spacing_distr_file = "{}.genome.active_pwm_scores.avg.txt.gz".format(
                plot_prefix)
            aligned_scores = results_adjusted[adj_key]["scores"]
            num_examples = aligned_scores.shape[0]
            aligned_scores = aligned_scores.mean(axis=0).transpose()
            aligned_scores.loc[0] = 0
            aligned_scores = aligned_scores.divide(aligned_scores.sum()).reset_index()
            aligned_scores.columns = ["position", "active"]
            aligned_scores.to_csv(
                score_spacing_distr_file,
                sep="\t", header=True, index=False, compression="gzip")
            plot_cmd = "{}/plot.spacing.freq.indiv.R {} {} {}".format(
                SCRIPT_DIR, score_spacing_distr_file, num_examples, plot_prefix)
            print plot_cmd
            os.system(plot_cmd)

            # plot signals
            for signal_key in signal_keys:
                signal_prefix = "{}.{}".format(plot_prefix, signal_key)
                tasks = results_adjusted[adj_key][signal_key].keys()
                for task in tasks:
                    task_prefix = "{}.{}".format(signal_prefix, task)
                    # save out
                    out_file = "{}.genome.txt.gz".format(
                        task_prefix)
                    results_adjusted[adj_key][signal_key][task].to_csv(
                        out_file,
                        sep="\t", header=True, index=False, compression="gzip")
                    # plot
                    plot_cmd = "Rscript {}/plot.spacing.signals.indiv.R {} {}".format(
                        SCRIPT_DIR, out_file, task_prefix)
                    print plot_cmd
                    
                    os.system(plot_cmd)
                    
            # and take region sets and get gene set enrichments
            syntax_enrichments = _get_gene_set_enrichment(
                results_adjusted[adj_key]["scores"].index.values,
                links_file,
                tss_file,
                background_gene_set_file,
                plot_prefix,
                OUT_DIR)
            syntax_enrichments["syntax"] = adj_key
            
            # add to summary
            if enrichments_summary is None:
                enrichments_summary = syntax_enrichments.copy()
            else:
                enrichments_summary = pd.concat([enrichments_summary, syntax_enrichments], axis=0)
                enrichments_summary = enrichments_summary.sort_values("term.name")
                
        # clean up summary
        enrichments_summary["log10pval"] = -np.log10(enrichments_summary["p.value"].values)
        enrichments_summary = enrichments_summary.drop("p.value", axis=1)
        summary_file = "{}/{}.summary.txt.gz".format(tmp_dir, pwm_name_clean)
        enrichments_summary.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")
                
        # remove substrings
        for substring in REMOVE_SUBSTRINGS:
            keep = [False if substring in term else True
                    for term in enrichments_summary["term.name"]]
            keep_indices = np.where(keep)[0]
            enrichments_summary = enrichments_summary.iloc[keep_indices]

        # remove exact strings
        keep = [False if term in REMOVE_EXACT_STRINGS else True
                for term in enrichments_summary["term.name"]]
        keep_indices = np.where(keep)[0]
        enrichments_summary = enrichments_summary.iloc[keep_indices]
        enrichments_summary = enrichments_summary.sort_values(["term.name", "log10pval"])
        summary_file = "{}/{}.summary.filt.txt.gz".format(OUT_DIR, pwm_name_clean)
        enrichments_summary.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")

        quit()

    return


main()
