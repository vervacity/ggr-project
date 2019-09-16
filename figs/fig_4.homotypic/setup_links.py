#!/usr/bin/env python

import os
import sys

import numpy as np
import pandas as pd

from ggr.analyses.linking import setup_proximity_links


def main():
    """set up links for downstream analyses
    """
    # args
    ref_regions_bed_file = sys.argv[1]
    tss_file = sys.argv[2]
    ref_regions_signal_file = sys.argv[3]
    rna_signal_file = sys.argv[4]
    OUT_DIR = sys.argv[5]

    # params
    k_nearest = 2
    max_dist = 25000 # Gasperini et al Cell 2019
    corr_coeff_thresh = -100 #0
    pval_thresh = 1 #0.10
    
    # set up links
    links_file = "{}/links.proximity.k-{}.d-{}.bed.gz".format(
        OUT_DIR, k_nearest, max_dist)
    setup_proximity_links(
        ref_regions_bed_file,
        tss_file,
        links_file,
        region_signal_file=ref_regions_signal_file,
        rna_signal_file=rna_signal_file,
        k_nearest=k_nearest,
        max_dist=max_dist,
        corr_coeff_thresh=corr_coeff_thresh,
        pval_thresh=pval_thresh)
    
    return

main()
