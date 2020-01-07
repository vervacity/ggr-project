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

    # -----------------------
    # dynamic ATAC
    # -----------------------
    
    # plot dynamic ATAC with matching profile maps of H3K27ac/H3K4me1/H3K27me3
    if False:
        plot_clusters(
            "{}/{}".format(GGR_DIR, outputs["results"]["atac"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"].lstrip("./")),
            "{}/{}".format(GGR_DIR, outputs["results"]["epigenome"]["dynamic"][
                "atac.epigenome_ordered.subsample.list"].lstrip("./")),
            "{}/{}".format(GGR_DIR, outputs["data"][
                "atac.counts.pooled.rlog.dynamic.traj.mat"].lstrip("./")),
            out_dir, "{}.epigenome.ATAC".format(prefix),
            plot_individual=True)
    # row seps file produced by plotting
    rowsep_file = "{}/{}.epigenome.ATAC.row_seps.txt".format(out_dir, prefix)
    
    if False:
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
            
    # -----------------------
    # stable ATAC
    # -----------------------

    stable_sets = [
        "dynamic_histones",
        "stable_histones"]
    stable_key_prefixes = [
        "atac.epigenome_ordered.dynamic_histones.clusters.subsample.mat",
        "atac.epigenome_ordered.stable_histones.clusters.subsample.mat"]
    new_prefixes = [
        "epigenome.ATAC.stable.dynamic_histones",
        "epigenome.ATAC.stable.stable_histones"]
    fake_clusters_key = "atac.stable.fake_cluster.list"

    # go through dynamic histones and stable histones
    for i in range(len(stable_key_prefixes)):
        stable_set = stable_sets[i]
        stable_key = stable_key_prefixes[i]
        stable_prefix = new_prefixes[i]
        
        # plot stable ATAC with matching profile maps of H3K27ac/H3K4me1/H3K27me3
        if False:
            plot_clusters(
                "{}/{}".format(GGR_DIR, outputs["results"]["epigenome"]["stable"][
                    fake_clusters_key].lstrip("./")),
                "{}/{}.plot.mat".format(GGR_DIR, outputs["results"]["epigenome"]["stable"][
                    stable_key].split(".mat")[0].lstrip("./")),
                "{}/{}".format(GGR_DIR, outputs["data"]["atac.counts.pooled.rlog.stable.mat"].lstrip("./")),
                out_dir, "{}.{}".format(prefix, stable_prefix),
                plot_individual=False)
        
        # row seps file produced by plotting
        rowsep_file = "{}/{}.{}.row_seps.txt".format(out_dir, prefix, stable_prefix)

        if False:
            # histones
            for histone_idx in range(len(histones)):
                histone = histones[histone_idx]
                histone_mat_file = "{}/{}/plots/ggr.epigenome.stable.{}_overlap.{}.point.mat.gz".format(
                    GGR_DIR,
                    outputs["results"]["epigenome"]["stable"]["dir"].lstrip("./"),
                    histone,
                    stable_set)
                histone_plot_file = "{}/{}.stable.{}.{}.profile_heatmap.pdf".format(
                    out_dir, prefix, histone, stable_set)
                plot_cmd = (
                    "plot.profile_heatmaps.R {} {} {} {} {} 1,100 101,200 201,300").format(
                        histone_mat_file, histone, rowsep_file, histone_plot_file,
                        histone_colors[histone_idx])
                print plot_cmd
                os.system(plot_cmd)
    
    # -----------------------
    # dynamic RNA
    # -----------------------

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
            plot_individual=True)
    
    return

main()
