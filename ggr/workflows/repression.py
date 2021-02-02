# workflow analyzing repression in skin differentiation


import os
import gzip
import glob
import logging

import numpy as np
import pandas as pd

from ggr.util.utils import run_shell_cmd

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import regions_to_genes_w_correlation_filtering



def runall(args, prefix):
    """Integrate epigenomic datasets
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("MASTER_WORKFLOW: run repression analyses")

    # set up data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "repression"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]

    # ---------------------------
    # SETUP: relevant files
    # ---------------------------
    
    # ATAC traj files
    atac_clusters_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    atac_mat_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.traj.mat"]
    atac_mat = pd.read_csv(atac_mat_file, sep="\t", index_col=0)
    atac_mat = atac_mat.drop("d05", axis=1)
    
    # RNA traj files
    rna_clusters_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    rna_mat_file = args.outputs["data"][
        "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
    rna_mat = pd.read_csv(rna_mat_file, sep="\t", index_col=0)
    
    # get linking file
    links_dir = "{}/linking/proximity".format(args.outputs["results"]["dir"])
    interactions_file = "{}/ggr.linking.ALL.overlap.interactions.txt.gz".format(links_dir)

    # ---------------------------
    # ANALYSIS: dynamic region sets - ATAC
    # ---------------------------

    # work dir
    work_dir = "{}/atac.dynamic".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(work_dir))

    # regions
    region_bed_dir = "{}/reproducible/hard/reordered/bed".format(
        args.outputs["results"]["atac"]["timeseries"]["dp_gp"]["dir"])
    region_sets = sorted(glob.glob("{}/*.bed.gz".format(region_bed_dir)))


    # for each region set, run analysis
    for region_set in region_sets:

        # debug
        continue
        
        # setup
        region_dir = "{}/{}".format(
            work_dir, os.path.basename(region_set).split(".bed")[0])
        run_shell_cmd("mkdir -p {}".format(region_dir))

        # make sure it has ATAC oriented IDs
        tmp_file = "{}/{}".format(region_dir, os.path.basename(region_set))
        regions = pd.read_csv(region_set, sep="\t", header=None)
        regions[3] = regions[0] + ":" + regions[1].astype(str) + "-" + regions[2].astype(str)
        regions.to_csv(tmp_file, sep="\t", compression="gzip", header=False, index=False)
        
        # get genes
        gene_file = "{}/genes.txt.gz".format(region_dir)
        regions_to_genes_w_correlation_filtering(
            tmp_file,
            interactions_file,
            args.outputs["annotations"]["tss.pc.bed"],
            gene_file,
            atac_mat,
            rna_mat,
            tmp_dir=region_dir,
            filter_by_score=0.5, # proximity corr thresh
            corr_thresh=0,
            corr_direction="negative")
        
        # and run enrichments
        run_gprofiler(
            gene_file,
            args.outputs["data"]["rna.counts.pc.expressed.mat"],
            region_dir, ordered=True, header=True)

        # homer analysis?
        
        # clean up
        #run_shell_cmd("rm {}".format(tmp_file))

    
    # ---------------------------
    # ANALYSIS: H3K27me3, ATAC-centric
    # ---------------------------
    
    # work dir
    work_dir = "{}/H3K27me3".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(work_dir))

    # regions
    region_bed_dir = "{}/clusters/by_state/bed".format(
        args.outputs["results"]["epigenome"]["stable"]["dir"])
    region_sets = [
        "{}/ggr.epigenome.stable.atac-cluster_0.state-cluster_997.bed.gz".format(
            region_bed_dir)]

    # use stable atac
    atac_mat_file = args.outputs["data"]["atac.counts.pooled.rlog.stable.mat"]
    atac_mat = pd.read_csv(atac_mat_file, sep="\t", index_col=0)
    atac_mat = atac_mat.drop("d05", axis=1)
    
    # for each region set, run analysis
    for region_set in region_sets:

        # debug
        continue

        # setup
        region_dir = "{}/{}".format(
            work_dir, os.path.basename(region_set).split(".bed")[0])
        run_shell_cmd("mkdir -p {}".format(region_dir))

        # make sure it has ATAC oriented IDs
        tmp_file = "{}/{}".format(region_dir, os.path.basename(region_set))
        regions = pd.read_csv(region_set, sep="\t", header=None)
        regions[3] = regions[0] + ":" + regions[1].astype(str) + "-" + regions[2].astype(str)
        regions.to_csv(tmp_file, sep="\t", compression="gzip", header=False, index=False)
        
        # get genes
        gene_file = "{}/genes.txt.gz".format(region_dir)
        regions_to_genes_w_correlation_filtering(
            tmp_file,
            interactions_file,
            args.outputs["annotations"]["tss.pc.bed"],
            gene_file,
            atac_mat,
            rna_mat,
            tmp_dir=region_dir,
            filter_by_score=0.5, # proximity corr thresh
            corr_thresh=0,
            corr_direction="positive")
        
        # and run enrichments
        run_gprofiler(
            gene_file,
            args.outputs["data"]["rna.counts.pc.expressed.mat"],
            region_dir, ordered=True, header=True)

        # homer analysis?
        
        # clean up
        #run_shell_cmd("rm {}".format(tmp_file))

    
    # ---------------------------
    # ANALYSIS: H3K27me3, dynamic
    # ---------------------------

    # adjust signal mats
    h3k27me3_mat_file = args.outputs["data"]["H3K27me3.counts.pooled.rlog.dynamic.mat"]
    h3k27me3_mat = pd.read_csv(h3k27me3_mat_file, sep="\t", index_col=0)
    rna_for_h3k27me3_mat = rna_mat[["d00", "d30", "d60"]]
    
    # take the dynamic clusters (not 9 or 4) and analyze those
    h3k27me3_clusters_bed = args.outputs[
        "results"]["histones"]["H3K27me3"]["timeseries"]["ggr.histone.H3K27me3.enumerated.bed"]
    h3k27me3_clusters = pd.read_csv(h3k27me3_clusters_bed, sep="\t", header=None)
    clusters = list(set(h3k27me3_clusters[3]))
    h3k27me3_clusters[4] = h3k27me3_clusters[0] + ":" + h3k27me3_clusters[1].astype(str) + "-" + h3k27me3_clusters[2].astype(str)

    for cluster_id in clusters:
        if cluster_id == 4:
            continue

        # setup
        region_dir = "{}/cluster_{}".format(work_dir, cluster_id)
        run_shell_cmd("mkdir -p {}".format(region_dir))
        
        # make a bed id file
        tmp_file = "{}/cluster_{}.bed.gz".format(region_dir, cluster_id)
        regions = h3k27me3_clusters[h3k27me3_clusters[3] == cluster_id]
        if regions.shape[0] < 100:
            continue
        regions = regions.drop(3, axis=1)
        regions.to_csv(tmp_file, sep="\t", compression="gzip", header=False, index=False)
    
        # get genes
        gene_file = "{}/genes.txt.gz".format(region_dir)
        regions_to_genes_w_correlation_filtering(
            tmp_file,
            interactions_file,
            args.outputs["annotations"]["tss.pc.bed"],
            gene_file,
            h3k27me3_mat,
            rna_for_h3k27me3_mat,
            tmp_dir=region_dir,
            filter_by_score=0.5, # proximity corr thresh
            corr_thresh=0,
            corr_direction="negative")
        
        # and run enrichments
        run_gprofiler(
            gene_file,
            args.outputs["data"]["rna.counts.pc.expressed.mat"],
            region_dir, ordered=True, header=True)

    
    

    quit()
    



    return args
