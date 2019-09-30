#!/usr/bin/env python

import os
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
    histone_colors = ["Reds", "Oranges", "Greens"]
    
    # plot dynamic ATAC with matching profile maps of H3K27ac/H3K4me1/H3K27me3
    if True:
        plot_clusters(
            "{}/{}".format(GGR_DIR, outputs["results"]["atac"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"].lstrip("./")),
            "{}/{}".format(GGR_DIR, outputs["results"]["epigenome"]["dynamic"][
                "atac.epigenome_ordered.subsample.list"].lstrip("./")),
            "{}/{}".format(GGR_DIR, outputs["data"][
                "atac.counts.pooled.rlog.dynamic.traj.mat"].lstrip("./")),
            out_dir, "{}.epigenome.ATAC".format(prefix),
            plot_individual=False)
    # row seps file produced by plotting
    rowsep_file = "{}/{}.epigenome.ATAC.row_seps.txt".format(out_dir, prefix)

    if True:
        # histones
        for histone_idx in range(len(histones)):
            histone = histones[histone_idx]
            histone_mat_file = "{}/{}/plots/ggr.epigenome.dynamic.{}_overlap.point.mat.gz".format(
                GGR_DIR,
                outputs["results"]["epigenome"]["dynamic"]["dir"].lstrip("./"),
                histone)
            histone_plot_file = "{}/{}.{}.dynamic.profile_heatmap.pdf".format(
                out_dir, prefix, histone)
            plot_cmd = (
                "plot.profile_heatmaps.R {} {} {} {} {} 1,100 101,200 201,300").format(
                    histone_mat_file, histone, rowsep_file, histone_plot_file,
                    histone_colors[histone_idx])
            print plot_cmd
            os.system(plot_cmd)

    # plot stable ATAC with matching profile maps of H3K27ac/H3K4me1/H3K27me3

    # plot dynamic RNA
    if True:
        plot_clusters(
            "{}/{}".format(GGR_DIR, outputs["results"]["rna"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"].lstrip("./")),
            "{}/{}".format(GGR_DIR, outputs["results"]["rna"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.subsample.list"].lstrip("./")),
            "{}/{}".format(GGR_DIR, outputs["data"][
                "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"].lstrip("./")),
            out_dir, "{}.RNA".format(prefix),
            plot_individual=False)
    
    return

main()
