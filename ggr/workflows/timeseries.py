# Description: workflow starting from a count matrix to trajectories

import os
import logging

from ggr.util.utils import run_shell_cmd

from ggr.analyses.filtering import filter_for_ids

from ggr.analyses.counting import split_count_matrix_by_replicate
from ggr.analyses.timeseries import get_consistent_dpgp_trajectories


def run_timeseries_workflow(args, prefix, mat_key="counts.mat"):
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
    assert args.inputs["params"].get("pairwise_deseq2_fdr") is not None
    logger.debug("using count matrix {}".format(args.outputs["data"][mat_key]))
    
    # inputs and outputs
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    logging.debug("data going to {}".format(data_dir))
    out_data = args.outputs["data"]
    
    results_dirname = "timeseries"
    results_dir = "{}/{}".format(args.outputs["results"]["rna"]["dir"], results_dirname)
    args.outputs["results"]["rna"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    logging.debug("results going to {}".format(results_dir))
    out_results = args.outputs["results"]["rna"][results_dirname]

    # ------------------------------------------------
    # ANALYSIS 0 - run pairwise DESeq to get dynamic IDs
    # input: count matrix
    # output: dynamic id set 
    # ------------------------------------------------
    logger.info("ANALYSIS: run pairwise deseq2 with multiple hypothesis correction")
    run_shell_cmd("mkdir -p {}/deseq2".format(results_dir))
    out_results["dynamic_ids.list"] = "{}/deseq2/{}.dynamic.ids.txt.gz".format(results_dir, prefix)
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
        split_count_matrix_by_replicate(out_data[mat_key], *[out_data[handle] for handle in matrix_handles])
        
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
        dynamic_handle = "{}.dynamic".format(rlog_handle)
        dynamic_mat_handles.append(dynamic_handle)
        out_data[dynamic_handle] = "{}.dynamic.mat.txt.gz".format(
            out_data[rlog_handle].split(".mat")[0])
        if not os.path.isfile(out_data[dynamic_handle]):
            filter_for_ids(
                out_data[rlog_handle],
                out_results["dynamic_ids.list"],
                out_data[dynamic_handle])

    # run trajectories using consistent DP_GP clustering
    # gives soft and hard clusters
    out_results["dpgp.clusters"] = get_consistent_dpgp_trajectories(
        [out_data[handle] for handle in dynamic_mat_handles],
        "{}/dp_gp".format(results_dir), # todo: register this folder?
        prefix,
        raw_cluster_min_size=args.inputs["params"]["rna_cluster_min_size"], # TODO convert this into a fract
        raw_cluster_reject_null_ci_interval=args.inputs["params"]["raw_cluster_reject_null_ci_interval"],
        rep_to_cluster_ci_interval=args.inputs["params"]["rep_to_cluster_ci_interval"],
        rep_to_cluster_corr_cutoff=args.inputs["params"]["rep_to_cluster_corr_cutoff"],
        epsilon=args.inputs["params"]["rep_to_cluster_epsilon"])

    quit()

    

    # hclust to get dendrogram, merge branches if they don't pass a statistical difference test? KL?
    
    
    
    # reorder in a deterministic (ish) way....


    # plotting
    

    # before return, make sure to save outputs to args



    
    return args
