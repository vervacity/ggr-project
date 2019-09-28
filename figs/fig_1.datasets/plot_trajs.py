#!/usr/bin/env python

import sys
import json

from ggr.analyses.timeseries import plot_clusters


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
    prefix = "fig_1"
    histones = ["H3K27ac", "H3K4me1", "H3K27me3"]
    histone_colors = ["red", "orange", "green"]
    
    # plot dynamic ATAC with matching profile maps of H3K27ac/H3K4me1/H3K27me3
    plot_clusters(
        "{}/{}".format(GGR_DIR, outputs["results"]["atac"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"].lstrip("./")),
        "{}/{}".format(GGR_DIR, outputs["results"]["epigenome"]["dynamic"][
            "atac.epigenome_ordered.subsample.list"].lstrip("./")),
        "{}/{}".format(GGR_DIR, outputs["data"][
            "atac.counts.pooled.rlog.dynamic.traj.mat"].lstrip("./")),
        out_dir, prefix)

    # histones
    rowsep_file = "{}/plots/ggr.epigenome.row_seps.txt".format(
        outputs["results"]["epigenome"]["dir"].lstrip("./"))
    for histone_idx in range(len(histones)):
        histone = histones[histone_idx]
        histone_mat_file = "{}/plots/ggr.epigenome.{}_overlap.point.mat.gz".format(
            outputs["results"]["epigenome"]["dir"].lstrip("./"),
            histone)
        histone_plot_file = "{}/{}.{}.profile_heatmap.pdf".format(
            out_dir, prefix, histone)
        plot_cmd = (
            "plot.profile_heatmaps.R {} {} {} {} 1,100 101,200 201,300").format(
                histone_mat_file, rowsep_file, histone_plot_file,
                histone_colors[histone_idx])
        print plot_cmd
        os.system(plot_cmd)
    
    quit()
    # plot stable ATAC with matching profile maps of H3K27ac/H3K4me1/H3K27me3


    # plot dynamic RNA
    
    
    return

main()
