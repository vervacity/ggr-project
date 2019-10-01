#!/usr/bin/env python

import os
import sys
import json

from ggr.analyses.linking import build_correlation_matrix
from ggr.analyses.linking import build_confusion_matrix


def main():
    """plot traj heatmaps
    """
    # pull in outputs json
    GGR_DIR = sys.argv[1]
    json_file = "{}/outputs.json".format(GGR_DIR)
    with open(json_file, "r") as fp:
        outputs = json.load(fp)
        
    # params
    out_dir = "."
    prefix = "fig_1.linking"

    # plot correlation matrix
    atac_clusters_file = "{}/{}".format(
        GGR_DIR,
        outputs["results"]["atac"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"].lstrip("./"))
    atac_mat_file = "{}/{}".format(
        GGR_DIR,
        outputs["data"]["atac.counts.pooled.rlog.dynamic.mat"].lstrip("./"))
    rna_clusters_file = "{}/{}".format(
        GGR_DIR,
        outputs["results"]["rna"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"].lstrip("./"))
    rna_mat_file = "{}/{}".format(
        GGR_DIR,
        outputs["data"][
            "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.mat"].lstrip("./"))

    # run analysis
    correlation_matrix_file = "{}/{}.correlation_mat.txt.gz".format(
        out_dir, prefix)
    build_correlation_matrix(
        atac_clusters_file,
        atac_mat_file,
        rna_clusters_file,
        rna_mat_file,
        correlation_matrix_file)

    # plot traj enrichment confusion matrix
    link_dir = "proximity"
    traj_dir = "{}/results/linking/{}/traj.linking".format(GGR_DIR, link_dir)
    build_confusion_matrix(
        atac_clusters_file,
        rna_clusters_file,
        traj_dir,
        "{}.traj_enrich.{}".format(prefix, link_dir))
    
    return

main()
