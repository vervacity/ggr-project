# Description: workflows for expression analyses

import os
import glob
import logging

from ggr.util.utils import run_shell_cmd
from ggr.analyses.filtering import filter_for_ids

from ggr.analyses.rna import make_rsem_matrix
from ggr.analyses.rna import threshold_empirically_off_genes
from ggr.analyses.rna import filter_clusters_for_tfs
from ggr.analyses.rna import filter_clusters_for_tss
from ggr.analyses.rna import get_nearest_regions
from ggr.analyses.rna import get_neighborhood
from ggr.analyses.rna import get_highly_expressed_genes
from ggr.analyses.rna import convert_gene_list_to_tss
from ggr.analyses.rna import add_clusters_to_tss_file
from ggr.analyses.rna import run_rdavid
from ggr.analyses.rna import run_gsea_on_series

from ggr.analyses.motifs import add_expressed_genes_to_metadata

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

    # TODO
    # annotation - pull in Andrew preprocessed scRNA-seq
    # basal vs differentiated state and create gene sigs (.gmt) for gsea downstream
    gsea_annot_dir = "{}/gsea".format(args.outputs["annotations"]["dir"])
    if not os.path.isdir(gsea_annot_dir):
        run_shell_cmd("mkdir -p {}".format(gsea_annot_dir))
        
        
    quit()
    
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

    # now that we have this, add this info to motifs
    pwm_metadata_key = "pwms.metadata.nonredundant"
    pwm_metadata_plus_expr_key = "{}.expressed".format(pwm_metadata_key)
    args.outputs["annotations"][pwm_metadata_plus_expr_key] = "{}.expressed.txt".format(
        args.outputs["annotations"][pwm_metadata_key].split(".txt")[0])
    if not os.path.isfile(args.outputs["annotations"][pwm_metadata_plus_expr_key]):
        add_expressed_genes_to_metadata(
            args.outputs["annotations"][pwm_metadata_key],
            args.outputs["annotations"][pwm_metadata_plus_expr_key],
            out_data["rna.counts.pc.expressed.mat"],
            args.outputs["annotations"]["geneids.mappings.mat"])

    # filter TSS file for expressed genes
    expressed_tss_key = "tss.pc.expressed"
    expressed_tss_file = "{}/{}.tss.expressed.bed.gz".format(
        args.outputs["annotations"]["dir"], prefix)
    out_results["expression_filtering"][expressed_tss_key] = expressed_tss_file
    if not os.path.isfile(expressed_tss_file):
        convert_gene_list_to_tss(
            args.outputs["results"]["rna"]["expression_filtering"]["gene_ids.expressed.list"],
            args.outputs["annotations"]["tss.pc.bed"],
            expressed_tss_file)

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

    # TODO - run GSEA on differential results (over d0 baseline)
    # produce 1 final results table that is gene_set x timepoint
    gsea_dir = "{}/gsea".format(
        args.outputs["results"][results_dirname]["timeseries"]["dir"])
    args.outputs["results"][results_dirname]["timeseries"]["gsea"] = {
        "dir": gsea_dir}
    #if not os.path.isdir(gsea_dir):
    if True:
        run_shell_cmd("mkdir -p {}".format(gsea_dir))
        deseq_files = sorted(glob.glob(
            "{}/deseq2/*over_d00_resultsAll.txt.gz".format(
                args.outputs["results"][results_dirname]["timeseries"]["dir"])))
        deseq_files = [filename for filename in deseq_files if "d05" not in filename]
        gsea_results_file = "{}/{}.gsea_results.txt.gz".format(gsea_dir, prefix)
        run_gsea_on_series(
            deseq_files,
            args.inputs["annot"][args.cluster]["gsea_gene_sets"],
            gsea_results_file,
            id_conversion_file=args.outputs["annotations"]["geneids.mappings.mat"],
            tmp_dir=gsea_dir)

    quit()
        
    

    
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

    # get BED files of TSS positions for each cluster
    cluster_tss_dir = "{}/reproducible/hard/reordered/tss".format(
        out_results["timeseries"]["dp_gp"]["dir"])
    cluster_tss_prefix = "{}/{}.hard.reordered".format(cluster_tss_dir, prefix)
    if not os.path.isdir(cluster_tss_dir):
        run_shell_cmd("mkdir {}".format(cluster_tss_dir))
        filter_clusters_for_tss(
            cluster_results[cluster_key],
            cluster_tss_prefix,
            args.outputs["annotations"]["tss.pc.bed"])

    # for these BED files get single nearest upstream accessible region (ie promoter)
    cluster_tss_nearest_dir = "{}/reproducible/hard/reordered/tss_nearest".format(
        out_results["timeseries"]["dp_gp"]["dir"])
    if not os.path.isdir(cluster_tss_nearest_dir):
        run_shell_cmd("mkdir {}".format(cluster_tss_nearest_dir))
        tss_files = glob.glob("{}/*.bed.gz".format(cluster_tss_dir))
        for tss_file in tss_files:
            out_file = "{}/{}.nearest.bed.gz".format(
                cluster_tss_nearest_dir,
                os.path.basename(tss_file).split(".bed")[0])
            get_nearest_regions(
                tss_file,
                out_data["atac.master.bed"],
                out_file,
                "-D a")

    # for these BED files get neighborhood of accessible regions
    cluster_tss_neighborhood_dir = "{}/reproducible/hard/reordered/tss_neighborhood".format(
        out_results["timeseries"]["dp_gp"]["dir"])
    if not os.path.isdir(cluster_tss_neighborhood_dir):
        run_shell_cmd("mkdir {}".format(cluster_tss_neighborhood_dir))
        tss_files = glob.glob("{}/*.bed.gz".format(cluster_tss_dir))
        for tss_file in tss_files:
            out_file = "{}/{}.neighborhood.bed.gz".format(
                cluster_tss_neighborhood_dir,
                os.path.basename(tss_file).split(".bed")[0])
            get_neighborhood(
                tss_file,
                out_data["atac.master.bed"],
                out_file,
                500000, # 1 Mb locus
                args.inputs["annot"][args.cluster]["chromsizes"])
        # also make global neighborhood of dynamic regions
        neighborhood_all = (
            "zcat {0}/*neighborhood.bed.gz | "
            "sort -k1,1 -k2,2n | "
            "bedtools merge -i stdin | "
            " gzip -c > {0}/{1}.hard.reordered.neighborhood.all.bed.gz").format(
                cluster_tss_neighborhood_dir,
                prefix)
        run_shell_cmd(neighborhood_all)

    # also get top decile of dynamic genes (decile based on MAX expression)
    # and pull those TSSs for use in downstream experimental results
    highly_expressed_dynamic_dir = "{}/reproducible/hard/reordered/highly_expressed".format(
        out_results["timeseries"]["dp_gp"]["dir"])
    if not os.path.isdir(highly_expressed_dynamic_dir):
        run_shell_cmd("mkdir {}".format(highly_expressed_dynamic_dir))
        # first pull a list of most highly expressed
        highly_expressed_genes_file = "{}/highly_expressed.genes.mat.txt.gz".format(
            highly_expressed_dynamic_dir)
        highly_expressed_tss_file = "{}.tss.bed.gz".format(highly_expressed_genes_file.split(".mat")[0])
        highly_expressed_tss_acc_file = "{}.acc.bed.gz".format(highly_expressed_tss_file.split(".bed")[0])
        get_highly_expressed_genes(
            "{}/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.mat.txt.gz".format(
                out_results["timeseries"]["dir"]),
            out_results["timeseries"]["dp_gp"]["clusters.reproducible.hard.reordered.list"],
            highly_expressed_genes_file,
            percentile=90)
        # then get TSS file
        convert_gene_list_to_tss(
            highly_expressed_genes_file,
            args.outputs["annotations"]["tss.pc.bed"],
            highly_expressed_tss_file)
        # and then get nearest ATAC peak
        get_nearest_regions(
            highly_expressed_tss_file,
            out_data["atac.master.bed"],
            highly_expressed_tss_acc_file,
            "-D a")

    # build a TSS file with trajectory numbers marked
    marked_tss_key = "tss.w_clusters"
    marked_tss_file = "{}/reproducible/hard/reordered/tss/{}.traj.tss.bed.gz".format(
        out_results["timeseries"]["dp_gp"]["dir"], prefix)
    out_results["timeseries"]["dp_gp"][marked_tss_key] = marked_tss_file
    if not os.path.isfile(marked_tss_file):
        add_clusters_to_tss_file(
            cluster_results[cluster_key],
            out_results["expression_filtering"][expressed_tss_key],
            marked_tss_file)
        
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




