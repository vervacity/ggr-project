# workflows integrating epigenomic datasets

import os
import glob
import logging

from ggr.util.utils import run_shell_cmd

from ggr.util.bioinformatics import make_deeptools_heatmap

from ggr.util.bed_utils import id_to_bed

from ggr.analyses.epigenome import get_histone_overlaps

from ggr.analyses.filtering import sort_by_clusters
from ggr.analyses.filtering import get_ordered_subsample

from ggr.analyses.timeseries import plot_clusters

def run_dynamic_epigenome_workflow(
        args,
        prefix,
        #out_dir,
        #prefix,
        #activating_histones,
        #activating_histone_files,
        #all_histones,
        histone_assignment="nearest"):
    """Order clusters by hclust and histones more simple
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get dynamic epigenome")

    # assertions
    
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "dynamic"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    # -----------------------------------------
    # ANALYSIS 0 - get overlap with histones
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    logger.info("ANALYSIS: overlap dynamic ATAC with histones")
    atac_dynamic_bed_key = "atac.dynamic.bed"
    histone_overlap_dir = "{}/overlap_histone".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(histone_overlap_dir))

    histone_overlap_mat_key = "atac.dynamic.overlap_histones.mat"
    out_results[histone_overlap_mat_key] = "{}/{}.overlaps.mat.txt.gz".format(
        histone_overlap_dir, prefix)

    # set up activating histone file sets
    histones = args.inputs["chipseq"][args.cluster]["histones"]["ordered_names"]
    histone_files = [
        (histone,
         args.outputs["results"]["histones"][histone]["timeseries"]["ggr.histone.{}.enumerated.bed".format(histone)],
         args.inputs["params"]["histones"][histone]["overlap_extend_len"])
        for histone in histones]
    activating_marks = histones[0:2]
    activating_mark_files = histone_files[0:2]
    
    if not os.path.isfile(out_results[histone_overlap_mat_key]):
        # NOTE: this currently requires bedtools 2.23.0
        # fix this for upward compatibility
        get_histone_overlaps(
            out_data[atac_dynamic_bed_key],
            activating_mark_files,
            histone_overlap_dir,
            out_results[histone_overlap_mat_key],
            args.inputs["annot"][args.cluster]["chromsizes"],
            histone_assignment=histone_assignment)

    # -----------------------------------------
    # ANALYSIS 1 - reorder ATAC with all clusters
    # inputs: dynamic ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: reorder ATAC regions using ATAC/histone clusters")
    atac_clusters_key = "clusters.reproducible.hard.reordered.list"
    atac_cluster_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][atac_clusters_key]
    atac_cluster_col = "cluster"

    histone_clusters_key = histone_overlap_mat_key
    histone_cluster_file = out_results[histone_clusters_key]
    histone_cluster_col = "histone_cluster"

    cluster_files = [
        (atac_cluster_file, atac_cluster_col),
        (histone_cluster_file, histone_cluster_col)]

    # reorder full list of ids
    atac_ordered_key = "atac.epigenome_ordered.list"
    out_results[atac_ordered_key] = "{}/{}.ordered.txt.gz".format(
        results_dir, prefix)
    atac_ordered_clusters_key = "atac.epigenome_ordered.mat"
    out_results[atac_ordered_clusters_key] = "{}/{}.ordered.clusters.txt".format(
        results_dir, prefix)

    if not os.path.isfile(out_results[atac_ordered_key]):
        sort_by_clusters(
            cluster_files,
            out_results[atac_ordered_clusters_key],
            out_results[atac_ordered_key])

    # and generate a subsample for plotting too
    atac_ordered_subsample_key = "{}.subsample.list".format(
        atac_ordered_key.split(".list")[0])
    out_results[atac_ordered_subsample_key] = "{}.subsampled.txt".format(
        out_results[atac_ordered_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_key]):
        get_ordered_subsample(
            out_results[atac_ordered_clusters_key],
            out_results[atac_ordered_subsample_key])

    # produce a BED too
    atac_ordered_subsample_bed_key = "{}.bed".format(
        atac_ordered_subsample_key.split(".list")[0])
    out_results[atac_ordered_subsample_bed_key] = "{}.bed".format(
        out_results[atac_ordered_subsample_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_bed_key]):
        id_to_bed(
            out_results[atac_ordered_subsample_key],
            "{}.gz".format(out_results[atac_ordered_subsample_bed_key]))

        # and unzip
        run_shell_cmd("gunzip {}.gz".format(out_results[atac_ordered_subsample_bed_key]))
        
    # -----------------------------------------
    # ANALYSIS 2 - plot this out
    # inputs: dynamic ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: plot out the heatmaps with the ordering in the subsample")
    mat_key = "atac.counts.pooled.rlog.dynamic.mat"
    plot_dir = "{}/plots".format(results_dir)
    
    # plot ATAC
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))
        plot_clusters(
            atac_cluster_file,
            out_results[atac_ordered_subsample_key],
            out_data[mat_key],
            plot_dir,
            prefix)
        # TODO need to plot color bars too?
        
        
        
        # plot the histone signal profiles with deeptools
        histone_colors = args.inputs["chipseq"][args.cluster]["histones"]["ordered_deeptools_colors"]
        for histone_idx in range(len(histones)):
            histone = histones[histone_idx]
            histone_color = histone_colors[histone_idx]
            histone_bigwigs = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_bigwig_glob"])))

            out_prefix = "{}/{}.{}_overlap".format(
                plot_dir,
                prefix,
                histone)
            out_file = "{}.heatmap.profile.pdf".format(out_prefix)
        
            if not os.path.isfile(out_file):
                make_deeptools_heatmap(
                    out_results[atac_ordered_subsample_bed_key],
                    histone_bigwigs,
                    out_prefix,
                    sort=False,
                    referencepoint="center",
                    color=histone_color)

    # don't forget to store outputs back
    
            
    return args


def run_stable_epigenome_workflow(args, prefix):
    """
    """


    return args


def runall(args, prefix):
    """Integrate epigenomic datasets
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("MASTER_WORKFLOW: run epigenomic analyses")

    # set up data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "epigenome"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]

    # -----------------------------------------
    # ANALYSIS 0 - integrate histones into dynamic ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    logger.info("ANALYSIS: integrate dynamic ATAC w histones")
    
    args = run_dynamic_epigenome_workflow(args, "{}.dynamic".format(prefix))

    # -----------------------------------------
    # ANALYSIS 1 - integrate histones into stable ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    
    
    
    # -----------------------------------------
    # ANALYSIS 2 - run bioinformatics tooks
    # inputs: BED files
    # outputs: TF enrichments, GO enrichments
    # -----------------------------------------
    



    

    return args
