#!/usr/bin/env python

import os
import re
import sys
import h5py

import numpy as np
import pandas as pd


def main():
    """get gene sets
    """
    # inputs
    OUT_DIR = sys.argv[1]
    sig_pwms_file = sys.argv[2]
    multiplicity_dir = sys.argv[3]
    spacing_dir = sys.argv[4]
    grammar_file = sys.argv[5]
    
    # dirs
    os.system("mkdir -p {}".format(OUT_DIR))
    tmp_dir = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(tmp_dir))
    
    # load sig pwms
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", index_col=0, header=0).index)

    # TODO load manual file
    grammars = pd.read_csv(grammar_file, sep="\t", header=0, index_col=0)
    
    # analyze per pwm and save out
    grammar_enrichment_summary = None
    for sig_idx in range(len(sig_pwms)):
        pwm_name = sig_pwms[sig_idx]
        pwm_name = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name = re.sub(".UNK.0.A", "", pwm_name)
        print sig_idx, pwm_name
        
        # summary files
        multiplicity_summary = "{}/{}.summary.filt.txt.gz".format(multiplicity_dir, pwm_name)
        spacing_summary = "{}/tmp/{}.summary.filt.txt.gz".format(spacing_dir, pwm_name)

        multiplicity_present = os.path.isfile(multiplicity_summary)
        spacing_present = os.path.isfile(spacing_summary)
        if not (multiplicity_present or spacing_present):
            continue
        
        # merge, melt? do this in R?
        plot_file = "{}/enrichments.{}.summary.pdf".format(OUT_DIR, pwm_name)
        plot_cmd = "Rscript ~/git/ggr-project/figs/fig_4.homotypic/plot.enrichments.R {} {} {}".format(multiplicity_summary, spacing_summary, plot_file)
        #print plot_cmd
        #os.system(plot_cmd)
        
        # if present in manual file, read in and select correct multiplicity and spacing
        # read in GO terms, filter, and merge
        # also need to get ATAC signal - should have gotten this earlier?
        if pwm_name in grammars.index.values:

            # check multiplicity, extract info
            multiplicity = grammars.loc[pwm_name]["multiplicity"]
            if not np.isnan(multiplicity):
                multiplicity_file = "{}/tmp/{}.counts.greater_equal-{}.gene_set.go_gprofiler.txt".format(
                    multiplicity_dir, pwm_name, int(multiplicity))
                multiplicity_data = pd.read_csv(multiplicity_file, sep="\t", header=0)
                keep = ["p.value", "term.id", "term.name", "domain"]
                multiplicity_data = multiplicity_data[keep]
                multiplicity_data["pwm_name"] = pwm_name
                multiplicity_data["multiplicity"] = int(multiplicity)
                multiplicity_data["spacing"] = np.nan
                
                if grammar_enrichment_summary is None:
                    grammar_enrichment_summary = multiplicity_data.copy()
                else:
                    grammar_enrichment_summary = pd.concat(
                        [grammar_enrichment_summary, multiplicity_data], axis=0)
                
            # check spacing, extract info
            spacing = str(grammars.loc[pwm_name]["spacing"])
            print spacing
            #if not np.isnan(spacing):
            if spacing != "nan":
                spacing_file = "{}/tmp/{}.range_{}.gene_set.go_gprofiler.txt".format(
                    spacing_dir, pwm_name, spacing)
                spacing_data = pd.read_csv(spacing_file, sep="\t", header=0)
                keep = ["p.value", "term.id", "term.name", "domain"]
                spacing_data = spacing_data[keep]
                spacing_data["pwm_name"] = pwm_name
                spacing_data["multiplicity"] = np.nan
                spacing_data["spacing"] = spacing
            
                if grammar_enrichment_summary is None:
                    grammar_enrichment_summary = spacing_data.copy()
                else:
                    grammar_enrichment_summary = pd.concat(
                        [grammar_enrichment_summary, spacing_data], axis=0)

    # filter for just GO/REAC
    grammar_enrichment_summary = grammar_enrichment_summary[
        (grammar_enrichment_summary["domain"] == "BP") |
        (grammar_enrichment_summary["domain"] == "CC") |
        (grammar_enrichment_summary["domain"] == "rea")]
    
    # save out grammar enrichment results
    grammar_summary_file = "{}/enrichments.GRAMMAR.txt.gz".format(OUT_DIR)
    grammar_enrichment_summary.to_csv(grammar_summary_file, sep="\t", compression="gzip")

    # and plot these
    
    return

main()
