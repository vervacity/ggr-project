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
    # TODO manual file for specific summaries
    
    # dirs
    tmp_dir = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(tmp_dir))
    
    # load sig pwms
    sig_pwms = list(pd.read_csv(sig_pwms_file, sep="\t", index_col=0, header=0).index)

    # TODO load manual file
    
    # analyze per pwm and save out
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
        plot_file = "{}/{}.summary.pdf".format(OUT_DIR, pwm_name)
        plot_cmd = "Rscript ~/git/ggr-project/figs/fig_4.homotypic/plot.enrichments.R {} {} {}".format(multiplicity_summary, spacing_summary, plot_file)
        print plot_cmd
        os.system(plot_cmd)
        
        # if present in manual file, read in and select correct multiplicity and spacing
        quit()

    return

main()
