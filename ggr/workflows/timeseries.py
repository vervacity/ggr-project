# description: workflows for timeseries clustering

import os
import glob
import logging

from ggr.util.utils import run_shell_cmd

from ggr.analyses.filtering import filter_for_ids
from ggr.analyses.filtering import get_ordered_subsample

from ggr.analyses.counting import split_count_matrix_by_replicate

from ggr.analyses.timeseries import run_dpgp
from ggr.analyses.timeseries import filter_null_and_small_clusters
from ggr.analyses.timeseries import get_reproducible_clusters
from ggr.analyses.timeseries import split_mat_to_clusters
from ggr.analyses.timeseries import plot_clusters
from ggr.analyses.timeseries import reorder_clusters
from ggr.analyses.timeseries import split_clusters


def run_cluster_plotting_subworkflow(
        out_data,
        out_results,
        out_dir,
        prefix,
        mat_key,
        cluster_list_key):
    """subsample and plot
    """
    # assertions
    assert out_data.get(mat_key) is not None
    assert out_results.get(cluster_list_key) is not None

    run_shell_cmd("mkdir -p {}".format(out_dir))
    
    # subsample
    cluster_list_subsample_key = "{}.subsample.list".format(
        cluster_list_key.split(".list")[0])
    out_results[cluster_list_subsample_key] = "{}.subsampled.txt".format(
        out_results[cluster_list_key].split(".txt")[0])
    if not os.path.isfile(out_results[cluster_list_subsample_key]):
        get_ordered_subsample(
            out_results[cluster_list_key],
            out_results[cluster_list_subsample_key])

    # plotting
    plot_dir = "{}/plots".format(out_dir)
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))
        plot_clusters(
            out_results[cluster_list_key],
            out_results[cluster_list_subsample_key],
            out_data[mat_key],
            "{}/plots".format(out_dir),
            prefix)

    # splitting
    individual_cluster_files = glob.glob("{}/*cluster_*".format(out_dir))
    if len(individual_cluster_files) == 0:
        split_clusters(out_results[cluster_list_key])

    return out_data, out_results


def run_cluster_reordering_subworkflow(
        out_data,
        out_results,
        out_dir,
        prefix,
        mat_key,
        cluster_list_key):
    """Reorder clusters based on hclust - try to order by time
    """
    # assertions
    assert out_data.get(mat_key) is not None
    assert out_results.get(cluster_list_key) is not None

    run_shell_cmd("mkdir -p {}".format(out_dir))
    
    # run reordering
    cluster_list_reordered_key = "{}.reordered.list".format(
        cluster_list_key.split(".list")[0])
    out_results[cluster_list_reordered_key] = "{}/{}.reordered.clustering.txt".format(
        out_dir, os.path.basename(out_results[cluster_list_key]).split(".clustering")[0])
    if not os.path.isfile(out_results[cluster_list_reordered_key]):
        reorder_clusters(
            out_results[cluster_list_key],
            out_data[mat_key],
            out_results[cluster_list_reordered_key])

    # run plotting subworkflow
    out_data, out_results = run_cluster_plotting_subworkflow(
        out_data,
        out_results,
        out_dir,
        prefix,
        mat_key,
        cluster_list_reordered_key)
    
    return out_data, out_results


def run_reproducibility_dpgp_workflow(
        args, 
        prefix,
        datatype_key="rna",
        pooled_mat_key="counts.pooled.mat",
        rep1_mat_key="counts.rep1.mat",
        rep2_mat_key="counts.rep2.mat"):
    """run DP-GP algorithm and then do reproducibility checks. 
    Requires replicates.
    """
    # logging and folder setup
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get reproducible DP-GP trajectories")
    
    # assertions
    assert args.outputs["data"].get(pooled_mat_key) is not None
    assert args.outputs["data"].get(rep1_mat_key) is not None
    assert args.outputs["data"].get(rep2_mat_key) is not None
    assert args.outputs["data"].get("dir") is not None
    assert args.outputs["results"][datatype_key]["timeseries"].get("dir") is not None

    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]
    
    results_dirname = "dp_gp"
    results_dir = "{}/{}".format(
        args.outputs["results"][datatype_key]["timeseries"]["dir"], 
        results_dirname)
    args.outputs["results"][datatype_key]["timeseries"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][datatype_key]["timeseries"][results_dirname]

    # ------------------------------------------------
    # ANALYSIS 0 - run DP GP
    # input: normalized count matrix
    # output: clustering files
    # ------------------------------------------------
    logger.info("ANALYSIS: run DP GP trajectory analysis")
    dpgp_pooled_dir = "{}/pooled".format(results_dir)
    clusters_raw_handle = "clusters.raw.list"
    out_results[clusters_raw_handle] = "{0}/{1}.pooled_optimal_clustering.txt".format(
        dpgp_pooled_dir, prefix)
    if not os.path.isfile(out_results[clusters_raw_handle]):
        run_dpgp(
            out_data[pooled_mat_key],
            "{}.pooled".format(prefix), 
            dpgp_pooled_dir, 
            results_dir)

    # subsample and plot
    out_data, out_results = run_cluster_plotting_subworkflow(
        out_data,
        out_results,
        dpgp_pooled_dir,
        "{}.pooled".format(prefix),
        pooled_mat_key,
        clusters_raw_handle)
    
    # reorder
    reorder_dir = "{}/reordered".format(dpgp_pooled_dir)
    out_data, out_results = run_cluster_reordering_subworkflow(
        out_data,
        out_results,
        reorder_dir,
        "{}.raw.reordered".format(prefix),
        pooled_mat_key,
        clusters_raw_handle)
    clusters_raw_reordered_handle = "{}.reordered.list".format(
        clusters_raw_handle.split(".list")[0])
    
    # ------------------------------------------------
    # ANALYSIS 1 - filter null and small clusters
    # input: cluster list and mat file
    # output: filtered clusters
    # ------------------------------------------------
    logger.info("ANALYSIS: filter out null and small clusters")
    nullfilt_dir = "{}/nullfilt".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(nullfilt_dir))
    clusters_nullfilt_handle = "clusters.null_filt.list"
    out_results[clusters_nullfilt_handle] = "{0}/{1}.nullfilt.clustering.txt".format(
        nullfilt_dir, prefix)
    if not os.path.isfile(out_results[clusters_nullfilt_handle]):
        filter_null_and_small_clusters(
            out_results[clusters_raw_reordered_handle],
            out_data[pooled_mat_key],
            out_results[clusters_nullfilt_handle],
            ci=0.95,
            size_cutoff=10)

    # and reorder
    reorder_dir = "{}/reordered".format(nullfilt_dir)
    out_data, out_results = run_cluster_reordering_subworkflow(
        out_data,
        out_results,
        reorder_dir,
        "{}.nullfilt.reordered".format(prefix),
        pooled_mat_key,
        clusters_nullfilt_handle)

    # ------------------------------------------------
    # ANALYSIS 2 - run reproducibility
    # input: cluster list, mat files (rep1, rep2, pooled)
    # output: filtered hard and soft clusters
    # ------------------------------------------------
    logger.info("ANALYSIS: check reproducibility of trajectories")
    reproducible_dir = "{}/reproducible".format(results_dir)
    clusters_reproducible_soft_handle = "clusters.reproducible.soft.list"
    out_results[clusters_reproducible_soft_handle] = (
        "{0}/soft/{1}.reproducible.soft.clustering.txt").format(
            reproducible_dir, prefix)
    clusters_reproducible_hard_handle = "clusters.reproducible.hard.list"
    out_results[clusters_reproducible_hard_handle] = (
        "{0}/hard/{1}.reproducible.hard.clustering.txt").format(
            reproducible_dir, prefix)
    if not os.path.isfile(out_results[clusters_reproducible_soft_handle]):
        get_reproducible_clusters(
            out_results[clusters_nullfilt_handle],
            out_data[pooled_mat_key],
            out_data[rep1_mat_key],
            out_data[rep2_mat_key],
            out_results[clusters_reproducible_soft_handle],
            out_results[clusters_reproducible_hard_handle],
            reproducible_dir,
            "{}.reproducible".format(prefix),
            ci=0.95,
            corr_cutoff=0.05,
            epsilon=0.10)

    # and reorder
    reorder_dir = "{}/hard/reordered".format(reproducible_dir)
    out_data, out_results = run_cluster_reordering_subworkflow(
        out_data,
        out_results,
        reorder_dir,
        "{}.reproducible.hard.reordered".format(prefix),
        pooled_mat_key,
        clusters_reproducible_hard_handle)

    # store data and results
    args.outputs["data"] = out_data
    args.outputs["results"][datatype_key]["timeseries"][results_dirname] = out_results
    
    return args


def run_timeseries_workflow(args, prefix, datatype_key="rna", mat_key="counts.mat"):
    """Run GGR timeseries workflow
    Use for both ATAC and RNA to have same 
    process applied to both

    Note that this workflow assumes a format of
    example_id x timepoint, where the bio reps are interleaved

    Args:
      args: a dictionary of filenames, etc
      params: all hyperparams used
    """
    # logging and folder setup
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: timeseries analysis")
    
    # assertions
    assert args.outputs["data"].get(mat_key) is not None
    assert args.outputs["data"].get("dir") is not None
    assert args.outputs["results"][datatype_key].get("dir") is not None
    assert args.inputs["params"].get("pairwise_deseq2_fdr") is not None
    logger.debug("using count matrix {}".format(args.outputs["data"][mat_key]))
    
    # inputs and outputs
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    logging.debug("data going to {}".format(data_dir))
    out_data = args.outputs["data"]
    
    results_dirname = "timeseries"
    results_dir = "{}/{}".format(args.outputs["results"][datatype_key]["dir"], results_dirname)
    args.outputs["results"][datatype_key][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    logging.debug("results going to {}".format(results_dir))
    out_results = args.outputs["results"][datatype_key][results_dirname]

    # ------------------------------------------------
    # ANALYSIS 0 - run pairwise DESeq to get dynamic IDs
    # input: count matrix
    # output: dynamic id set 
    # ------------------------------------------------
    logger.info("ANALYSIS: run pairwise deseq2 with multiple hypothesis correction")
    run_shell_cmd("mkdir -p {}/deseq2".format(results_dir))
    out_results["dynamic_ids.list"] = "{}/deseq2/{}.dynamic.ids.txt.gz".format(
        results_dir, prefix)
    if not os.path.isfile(out_results["dynamic_ids.list"]):
        run_timeseries_deseq2 = "timeseries.pairwise_deseq2.R {0} {1} {2} {3} full".format(
            out_data[mat_key],
            "{0}/deseq2/{1}".format(results_dir, prefix),
            args.inputs["params"]['pairwise_deseq2_fdr'],
            out_results['dynamic_ids.list'])
        run_shell_cmd(run_timeseries_deseq2)
        
    # ------------------------------------------------
    # ANALYSIS 1 - split into replicates and pooled
    # input: count matrix
    # output: pooled, rep1, rep2 count matrices
    # ------------------------------------------------
    logger.info("ANALYSIS: Separating count mat into reps/pooled...")
    matrix_types = ["rep1", "rep2", "pooled"]
    counts_prefix = out_data[mat_key].split(".mat")[0]
    matrix_handles = []
    for matrix_type in matrix_types:
        handle = "{}.{}.mat".format(mat_key.split(".mat")[0], matrix_type)
        matrix_handles.append(handle)
        out_data[handle] = "{}.{}.mat.txt.gz".format(counts_prefix, matrix_type)

    pooled_handle = "{}.{}.mat".format(mat_key.split(".mat")[0], matrix_types[-1])
    if not os.path.isfile(out_data[pooled_handle]):
        split_count_matrix_by_replicate(
            out_data[mat_key],
            *[out_data[handle] for handle in matrix_handles])
        
    # ------------------------------------------------
    # ANALYSIS 2 - normalize by rlog
    # input: count matrices
    # output: rlog normalized pooled, rep1, rep2 count matrices
    # ------------------------------------------------
    logger.info("ANALYSIS: normalizing count files...")
    pooled_rlog_mat_file = "{}.rlog.mat.txt.gz".format(out_data[pooled_handle].split(".mat")[0])
    if not os.path.isfile(pooled_rlog_mat_file):
        run_rlogs = "normalize.rlog_master_and_transfer.R {0} {1} {2} {3}".format(
            out_data[mat_key], *[out_data[handle] for handle in matrix_handles])
        run_shell_cmd(run_rlogs)
        
    # ------------------------------------------------
    # ANALYSIS 3 - split matrices to only keep the dynamic ids
    # input: normalized count matrices
    # output: normalized count matrices with just dynamic ids
    # ------------------------------------------------ 
    logger.info("ANALYSIS: splitting dynamic and stable...")
    dynamic_mat_handles = []
    for matrix_handle in matrix_handles:

        # set up handle
        rlog_handle = "{}.rlog.mat".format(matrix_handle)
        out_data[rlog_handle] = "{}.rlog.mat.txt.gz".format(
            out_data[matrix_handle].split('.mat')[0])

        # filter for dynamic ids (NOTE: set up BED file earlier)
        dynamic_handle = "{}.dynamic.mat".format(rlog_handle)
        dynamic_mat_handles.append(dynamic_handle)
        out_data[dynamic_handle] = "{}.dynamic.mat.txt.gz".format(
            out_data[rlog_handle].split(".mat")[0])
        if not os.path.isfile(out_data[dynamic_handle]):
            filter_for_ids(
                out_data[rlog_handle],
                out_results["dynamic_ids.list"],
                out_data[dynamic_handle])

    # store outputs in prep for DP GP clustering
    args.outputs["data"] = out_data
    args.outputs["results"][datatype_key][results_dirname] = out_results

    # ------------------------------------------------
    # ANALYSIS 4 - run reproducible dpgp workflow
    # input: normalized count matrices
    # output: dynamic clusters
    # ------------------------------------------------ 
    key = "{}.counts.pc.pooled.rlog.dynamic.mat".format(datatype_key)
    args = run_reproducibility_dpgp_workflow(
        args, 
        prefix,
        datatype_key=datatype_key,
        pooled_mat_key=dynamic_mat_handles[2],
        rep1_mat_key=dynamic_mat_handles[0],
        rep2_mat_key=dynamic_mat_handles[1])
    
    if False:
        # ------------------------------------------------
        # ANALYSIS 4b - merge "similar" clusters? do this for both hard/soft clusters?
        # input: cluster list, mat files (rep1, rep2, pooled)
        # output: filtered hard and soft clusters
        # notes: impose hclust on the clusters, merge from
        #  leaves until they are no longer similar (ie distribution
        #  test)
        # ------------------------------------------------
        # bootstrap test 
        merged_dir = "{}/merged".format(results_dir)
        clusters_merged_handle = "clusters.merged.list"
        out_results[clusters_merged_handle] = "{0}/{1}.merged.clustering.txt".format(
            merged_dir, prefix)
        if not os.path.isfile(out_results[clusters_merged_handle]):
            pass

    # don't use below at the moment, since last return point was args
    # before return, make sure to save outputs to args
    #args.outputs["data"] = out_data
    #args.outputs["results"][datatype_key][results_dirname] = out_results
    
    return args
