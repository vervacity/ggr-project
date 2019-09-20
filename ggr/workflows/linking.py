"""ggr workflow for atac analyses
"""

import os
import glob
import signal
import logging

from ggr.util.utils import run_shell_cmd
from ggr.util.utils import parallel_copy

from ggr.util.bed_utils import merge_regions
from ggr.util.bed_utils import id_to_bed

from ggr.analyses.linking import setup_proximity_links

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import build_correlation_matrix
from ggr.analyses.linking import traj_bed_to_gene_set_by_proximity
from ggr.analyses.linking import build_confusion_matrix
from ggr.analyses.linking import get_replicate_consistent_links



def run_replicate_consistent_links_workflow(
        args,
        prefix,
        results_dir,
        link_key="distance",
        out_dir="abc.distance"):
    """take in ABC links and get replicate consistent link set
    """
    # set up
    link_dir = args.inputs["conformation"][args.cluster][link_key]["data_dir"]
    link_days = ["d0", "d3", "d6"]
    replicated_link_dir = "{}/{}".format(results_dir, out_dir)
    os.system("mkdir -p {}".format(replicated_link_dir))

    # go through each day
    for link_day in link_days:
        pooled_file = glob.glob(
            "{}/*{}*pooled*.txt.gz".format(link_dir, link_day))
        print pooled_file
        assert len(pooled_file) == 1
        pooled_file = pooled_file[0]
        rep_link_files = sorted(
            glob.glob("{}/*{}*b*.gz".format(link_dir, link_day)))
        link_prefix = "{}/{}.{}".format(
            replicated_link_dir,
            prefix,
            link_day)
        get_replicate_consistent_links(
            pooled_file,
            rep_link_files,
            link_prefix)

    return args


def runall(args, prefix):
    """all workflows for atac-seq data
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: run linking analyses")

    # set up data
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]
    
    results_dirname = "linking"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]
    
    # -------------------------------------------
    # ANALYSIS - naive proximity linking
    # input: ATAC regions and TSS file
    # output: linking by naive proximity
    # -------------------------------------------
    logger.info("ANALYSIS: naive proximity linking")
    k_nearest = 3
    max_dist = 25000 # Gasperini et al Cell 2019
    corr_coeff_thresh = 0 #0
    pval_thresh = 0.1 #0.10
    
    # set up links <- do it day specific?
    links_file = "{}/links.proximity.k-{}.d-{}.bed.gz".format(
        OUT_DIR, k_nearest, max_dist)
    setup_proximity_links(
        ref_regions_bed_file,
        tss_file,
        links_file,
        region_signal_file=ref_regions_signal_file,
        rna_signal_file=rna_signal_file,
        k_nearest=k_nearest,
        max_dist=max_dist,
        corr_coeff_thresh=corr_coeff_thresh,
        pval_thresh=pval_thresh)


    quit()
    
    # -------------------------------------------
    # ANALYSIS - different styles of ABC linking
    # input: ABC links
    # output: replicate consistent link sets
    # -------------------------------------------
    logger.info("ANALYSIS: build replicate consistent ABC links - distance-based")
    if False:
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="distance",
            out_dir="abc.distance")

    logger.info("ANALYSIS: build replicate consistent ABC links - cell type avg HiC")
    run_replicate_consistent_links_workflow(
        args,
        prefix,
        results_dir,
        link_key="hic.celltype_avg",
        out_dir="abc.hic.celltype_avg")


    quit()


    
    # -------------------------------------------
    # ANALYSIS - get correlation matrix between ATAC and RNA
    # input: atac clusters/data, rna clusters/data
    # output: correlation matrix
    # -------------------------------------------
    logger.info("ANALYSIS: correlation between ATAC/RNA clusters")

    # set up dir
    atac_linking_dir = "{}/atac".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(atac_linking_dir))

    # get needed files
    atac_clusters_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    atac_mat_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.mat"]

    rna_clusters_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    rna_mat_file = args.outputs["data"][
        "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.mat"]

    # run analysis
    correlation_matrix_file = "{}/{}.correlation_mat.txt.gz".format(
        atac_linking_dir, prefix)
    if not os.path.isfile(correlation_matrix_file):
        build_correlation_matrix(
            atac_clusters_file,
            atac_mat_file,
            rna_clusters_file,
            rna_mat_file,
            correlation_matrix_file)
    
    # -------------------------------------------
    # ANALYSIS - overlap ATAC trajectories with RNA
    # input: atac trajectories, tss file
    # output: gene sets, enrichments, confusion matrix
    # -------------------------------------------
    logger.info("ANALYSIS: try naive proximal linking")

    # use the TSS file that only has expressed genes
    tss_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"]["tss.w_clusters"]
    
    # generate the gene sets for each ATAC trajectory
    atac_traj_dir = "{}/reproducible/hard/reordered/bed".format(
        args.outputs["results"]["atac"]["timeseries"]["dp_gp"]["dir"])
    traj_bed_files = sorted(glob.glob("{}/*cluster*bed.gz".format(atac_traj_dir)))
    gene_set_files = []
    for traj_bed_file in traj_bed_files:
        gene_set_file = "{}/{}.linked_genes.txt.gz".format(
            atac_linking_dir,
            os.path.basename(traj_bed_file).split(".bed")[0])
        if not os.path.isfile(gene_set_file):
            bed_to_gene_set_by_proximity(
                traj_bed_file,
                tss_file,
                gene_set_file,
                k_nearest=2,
                max_dist=100000)
        gene_set_files.append(gene_set_file)
        
    # run enrichments (gprofiler)
    background_gene_set_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    atac_linked_genes_enrich_dir = "{}/enrichments".format(atac_linking_dir)
    if not os.path.isdir(atac_linked_genes_enrich_dir):
        run_shell_cmd("mkdir -p {}".format(atac_linked_genes_enrich_dir))
        for gene_set_file in gene_set_files:
            run_gprofiler(
                gene_set_file,
                background_gene_set_file,
                atac_linked_genes_enrich_dir)
    
    # collect into a confusion matrix
    cluster_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"]["clusters.reproducible.hard.reordered.list"]
    confusion_matrix_file = "{}/{}.confusion_mat.txt.gz".format(atac_linking_dir, prefix)
    if not os.path.isfile(confusion_matrix_file):
        build_confusion_matrix(
            traj_bed_files,
            gene_set_files,
            cluster_file,
            confusion_matrix_file)

    # TODO split this by TFs vs non TFs <- make annotation TSS files

    # -------------------------------------------
    # ANALYSIS - get replicate consistent distance based links (ABC method)
    # input: ABC distance-based links
    # output: replicate consistent links
    # -------------------------------------------
    logger.info("ANALYSIS: build replicate consistent ABC links - distance-based")
    if False:
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="distance",
            out_dir="abc.distance")

    logger.info("ANALYSIS: build replicate consistent ABC links - cell type avg HiC")
    run_replicate_consistent_links_workflow(
        args,
        prefix,
        results_dir,
        link_key="hic.celltype_avg",
        out_dir="abc.hic.celltype_avg")

    
    quit()
    
    # set up
    link_dir = args.inputs["conformation"][args.cluster]["distance"]["data_dir"]
    link_days = ["d0", "d3", "d6"]
    replicated_link_dir = "{}/abc.distance".format(results_dir)
    os.system("mkdir -p {}".format(replicated_link_dir))
    
    for link_day in link_days:
        pooled_file = glob.glob(
            "{}/*{}*pooled*.txt.gz".format(link_dir, link_day))
        print pooled_file
        assert len(pooled_file) == 1
        pooled_file = pooled_file[0]
        rep_link_files = sorted(
            glob.glob("{}/*{}*b*.gz".format(link_dir, link_day)))
        link_prefix = "{}/{}.{}".format(
            replicated_link_dir,
            prefix,
            link_day)
        get_replicate_consistent_links(
            pooled_file,
            rep_link_files,
            link_prefix)


        
    # TODO copy over to data dir?
    quit()
    
    # -------------------------------------------
    # ANALYSIS - get replicate consistent HiC based links (ABC method)
    # input: ABC HiC links
    # output: replicate consistent links
    # -------------------------------------------
    logger.info("ANALYSIS: build replicate consistent HiChIP link set")

    # set up
    link_dir = args.inputs["conformation"][args.cluster]["hic.keratinocyte"]["data_dir"]
    #link_days = ["d0", "d3", "d6"]
    link_days = ["d3", "d6"]
    replicated_link_dir = "{}/abc.hic.keratinocyte".format(results_dir)
    os.system("mkdir -p {}".format(replicated_link_dir))
    
    # get replicate consistent links
    for link_day in link_days:
        pooled_file = glob.glob(
            "{}/*{}*pooled*.txt.gz".format(link_dir, link_day))
        assert len(pooled_file) == 1
        pooled_file = pooled_file[0]
        rep_link_files = sorted(
            glob.glob("{}/*{}*b*.gz".format(link_dir, link_day)))
        link_prefix = "{}/{}.{}".format(
            replicated_link_dir,
            prefix,
            link_day)
        get_replicate_consistent_links(
            pooled_file,
            rep_link_files,
            link_prefix)
        
    quit()
    

    
    # and then rerun select enrichments?
    

    return args
