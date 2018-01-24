# Description: workflows for expression analyses

import os
import glob
import logging

from ggr.util.utils import run_shell_cmd
from ggr.analyses.filtering import filter_for_ids

from ggr.analyses.rna import make_rsem_matrix
from ggr.analyses.rna import threshold_empirically_off_genes
from ggr.analyses.rna import filter_clusters_for_tfs
from ggr.analyses.rna import run_rdavid

from ggr.workflows.timeseries import run_timeseries_workflow


def run_rna_qc_workflow(args, params, out_dir, tmp_dir, prefix):
    """run QC
    """

    # PCA

    # correlation


    # heatmap of highly variable genes?

    

    return


def run_expressed_threshold_workflow(args, prefix, mat_key="rna.counts.pc.mat"):
    """Workflow for determining an empirical "ON" 
    (ie expressed in at least one sample) threshold for the sample.
    Requires an input count matrix, ideally with replicates 
    (interleaved columns)

    Args:
      args: input dictionary of files
      mat_key: key to the matrix that you want to calculate the empirical thresh
      out_dir: where to save the final output matrix
      prefix: prefix to append to tmp and final files
    """
    # logging and folder setup
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: filter out empirically OFF genes")

    # assertions
    assert args.outputs["data"].get(mat_key) is not None
    assert args.outputs["data"].get("dir") is not None
    assert args.inputs["params"].get("rna_empirical_rlog_thresh") is not None

    # inputs and outputs
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    logging.debug("data going to {}".format(data_dir))
    out_data = args.outputs["data"]

    results_dirname = "expression_filtering"
    results_dir = "{}/{}".format(args.outputs["results"]["rna"]["dir"], results_dirname)
    args.outputs["results"]["rna"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    logging.debug("results going to {}".format(results_dir))
    out_results = args.outputs["results"]["rna"][results_dirname]
    
    # rlog
    rlog_key = "{}.rlog.mat".format(mat_key.split(".mat")[0])
    out_results[rlog_key] = "{}/{}.rlog.mat.txt.gz".format(
        results_dir, os.path.basename(out_data[mat_key]).split(".mat")[0])
    if not os.path.isfile(out_results[rlog_key]):
        rlog_norm = "normalize.rlog.R {0} {1}".format(
            out_data[mat_key],
            out_results[rlog_key])
        run_shell_cmd(rlog_norm)

    # and then look at the empirical distribution
    density_plot_key = "{}.distr.plot".format(rlog_key.split(".mat")[0])
    out_results[density_plot_key] = "{}.density.pdf".format(
        out_results[rlog_key].split(".mat")[0])
    if not os.path.isfile(out_results[density_plot_key]):
        plot_distribution = "viz.plot_quant_distribution.R {} {} {}".format(
            out_results[rlog_key],
            args.inputs["params"]["rna_empirical_rlog_thresh"],
            out_results[density_plot_key])
        run_shell_cmd(plot_distribution)
        
    # and now threshold out genes that fall below this
    rlog_expressed_key = "{}.expressed.mat".format(rlog_key.split(".mat")[0])
    out_results[rlog_expressed_key] = "{}.expressed.mat.txt.gz".format(
        out_results[rlog_key].split(".mat")[0])
    if not os.path.isfile(out_results[rlog_expressed_key]):
        threshold_empirically_off_genes(
            out_results[rlog_key],
            out_results[rlog_expressed_key],
            thresh=args.inputs["params"]["rna_empirical_rlog_thresh"])
        
    # make into a gene id list
    geneids_expressed_key = "gene_ids.expressed.list" # TODO fix this key
    out_results[geneids_expressed_key] = "{}.txt.gz".format(
        out_results[rlog_expressed_key].split(".mat")[0])
    if not os.path.isfile(out_results[geneids_expressed_key]):
        get_ids = (
            "zcat {} | "
            "awk -F '\t' '{{ print $1 }}' | "
            "grep -v -e '^$' | "
            "gzip -c > {}").format(
                out_results[rlog_expressed_key],
                out_results[geneids_expressed_key])
        run_shell_cmd(get_ids)

    # and use to filter out the initial input count matrix
    expressed_key = "{}.expressed.mat".format(mat_key.split(".mat")[0])
    out_data[expressed_key] = "{}.expressed.mat.txt.gz".format(
        out_data[mat_key].split(".mat")[0])
    if not os.path.isfile(out_data[expressed_key]):
        filter_for_ids(
            out_data[mat_key],
            out_results[geneids_expressed_key],
            out_data[expressed_key])

    return args


def runall(args, prefix):
    """all workflows for expression data
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("MASTER_WORKFLOW: run all rna analyses")

    # pull out specific inputs
    inputs = args.inputs["rna"][args.cluster]
    
    # assertions
    assert inputs.get("data_dir") is not None
    assert inputs.get("rsem_glob") is not None

    # set up data dir
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]

    # set up results dir
    results_dirname = "rna"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]
    
    # ----------------------------------------------------
    # ANALYSIS 0 - create data matrices of expected counts
    #  from RSEM quant files. Remove media influenced
    #  timepoints.
    # input: RSEM files
    # output: matrix of counts and TPMs
    # ----------------------------------------------------
    logger.info("ANALYSIS: get matrices from RSEM quant files")
    matrix_types = ["rna.counts", "rna.tpm"]
    colnames = ["expected_count", "TPM"]
    rna_quant_files = sorted(
        glob.glob("{0}/{1}".format(
            inputs["data_dir"],
            inputs["rsem_glob"])))
    assert len(rna_quant_files) > 0
    # remove files from media influenced timepoints
    for timepoint_string in args.inputs["params"]["media_timepoints"]:
        rna_quant_files = [
            filename for filename in rna_quant_files
            if timepoint_string not in filename]
    rna_matrices = []
    for i in xrange(len(matrix_types)):
        matrix_type = matrix_types[i]
        colname = colnames[i]
        handle = "{}.mat".format(matrix_type)
        rna_matrices.append(handle)
        out_data[handle] = "{}/{}.{}.mat.txt.gz".format(
            data_dir, prefix, matrix_type.split(".")[1])
        if not os.path.isfile(out_data[handle]):
            make_rsem_matrix(
                rna_quant_files,
                out_data[handle],
                colname=colname)
            
    # ----------------------------------------------------
    # ANALYSIS 1 - filter for just protein coding genes
    # input: count matrix
    # output: count matrix of just protein coding
    # ----------------------------------------------------
    logger.info("ANALYSIS: filter matrices for protein coding genes")
    rna_pc_matrices = []
    for matrix_type in matrix_types:
        handle = "{}.mat".format(matrix_type)
        pc_handle = "{}.pc.mat".format(matrix_type)
        rna_pc_matrices.append(pc_handle)
        out_data[pc_handle] = "{}.pc.mat.txt.gz".format(
            out_data[handle].split(".mat")[0])
            #data_dir, prefix, matrix_type)
        if not os.path.isfile(out_data[pc_handle]):
            filter_for_ids(
                out_data[handle],
                args.outputs["annotations"]["geneids.pc.list"],
                out_data[pc_handle])

    # ----------------------------------------------------
    # ANALYSIS 2 - filter for empirically ON genes
    # input: count matrix of protein coding
    # output: count matrix of ON genes 
    # ----------------------------------------------------
    args = run_expressed_threshold_workflow(
        args,
        prefix,
        mat_key="rna.counts.pc.mat")
    
    # ----------------------------------------------------
    # ANALYSIS 3 - run timeseries analysis on these genes
    # input: count matrix of expressed protein coding genes
    # output: gene trajectories
    # ----------------------------------------------------
    args = run_timeseries_workflow(
        args,
        prefix,
        datatype_key="rna",
        mat_key="rna.counts.pc.expressed.mat")

    # TODO: maybe need to have a link to the cluster path and dir?
    cluster_key = "clusters.reproducible.hard.reordered.list"
    out_data = args.outputs["data"]
    out_results = args.outputs["results"][results_dirname]
    
    # ------------------------------------------------
    # ANALYSIS 4 - mark key TFs/biomarkers per cluster
    # input: cluster file
    # output: key TFs per cluster file
    # ------------------------------------------------ 
    cluster_results = out_results["timeseries"]["dp_gp"]
    cluster_tfs_key = "{}.tfs_only.list".format(cluster_key.split(".list")[0])
    cluster_results[cluster_tfs_key] = (
        "{}.tfs_only.hgnc.clustering.txt").format(
        cluster_results[cluster_key].split(".clustering")[0])
    if not os.path.isfile(cluster_results[cluster_tfs_key]):
        filter_list = args.inputs["params"]["known_tfs_hgnc"]
        filter_list += args.inputs["params"]["known_biomarkers_hgnc"]
        filter_clusters_for_tfs(
            cluster_results[cluster_key],
            cluster_results[cluster_tfs_key],
            args.outputs["annotations"]["geneids.mappings.mat"],
            filter_list)
    out_results["timeseries"]["dp_gp"] = cluster_results
        
    # ------------------------------------------------
    # ANALYSIS 5 - Gene Ontology enrichments
    # input: gene clusters
    # output: annotations
    # ------------------------------------------------
    single_cluster_files = sorted(
        glob.glob("{}/reproducible/hard/reordered/*cluster_*txt.gz".format(
            out_results["timeseries"]["dp_gp"]["dir"])))
    
    # rDAVID
    cluster_dir = "{}/reproducible/hard/reordered".format(
        out_results["timeseries"]["dp_gp"]["dir"])
    go_dir = "{}/go.rdavid".format(cluster_dir)
    if not os.path.isdir(go_dir):
        run_shell_cmd("mkdir -p {}".format(go_dir))
        for single_cluster_file in single_cluster_files:
            run_rdavid(
                single_cluster_file,
                args.outputs["annotations"]["geneids.pc.list"],
                go_dir)
    
    logger.info("MASTER_WORKFLOW: DONE")

    # for figures:
    # supplements: GO enrichment plots
    # supplements QC
    # main figure - promoter analyses: show ATAC, histones (this might go into supplements)
    # transcription factor sub heatmap and line plots for key TFs (basically ordered TF cascade for known TFs)
    
    # eventually 3d linking to enhancers
    
    
    return args




