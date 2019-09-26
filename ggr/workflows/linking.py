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

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import build_correlation_matrix
from ggr.analyses.linking import link_by_distance
from ggr.analyses.linking import link_by_distance_and_corr
from ggr.analyses.linking import get_replicate_consistent_links
from ggr.analyses.linking import get_timepoint_consistent_links
from ggr.analyses.linking import region_clusters_to_genes

# TODO deprecate these
from ggr.analyses.linking import traj_bed_to_gene_set_by_proximity
from ggr.analyses.linking import build_confusion_matrix


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
            link_prefix,
            idr_thresh=args.inputs["params"]["linking"][out_dir]["idr_thresh"])

    # and reconcile across timepoints
    all_prefix = "{}/{}.ALL".format(replicated_link_dir, prefix)
    all_interactions_file = "{}.interactions.txt.gz".format(all_prefix)
    if not os.path.isfile(all_interactions_file):
        links_files = sorted(
            glob.glob("{}/*interactions.txt.gz".format(
                replicated_link_dir)))
        get_timepoint_consistent_links(
            links_files, all_prefix)

    return args


def run_traj_linking_workflow(args, prefix, links_dir):
    """link ATAC trajectories to RNA trajectories with linking file
    """
    logging.info("ANALYSIS: correlation between ATAC/RNA clusters")

    # setup
    traj_dir = "{}/traj.linking".format(links_dir)
    os.system("mkdir -p {}".format(traj_dir))

    # get needed files
    atac_clusters_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    atac_mat_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.mat"]

    rna_clusters_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    # TODO consider using ALL expressed genes?
    rna_mat_file = args.outputs["data"][
        "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
    rna_expressed_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    
    # get linking file
    interactions_file = "{}/{}.ALL.overlap.interactions.txt.gz".format(links_dir, prefix)

    # tss file
    tss_file = args.outputs["annotations"]["tss.pc.bed"]
    
    # first, per cluster, get gene set enrichments
    region_clusters_to_genes(
        atac_clusters_file,
        interactions_file,
        tss_file,
        traj_dir,
        run_enrichments=True,
        background_gene_file=rna_expressed_file,
        filter_by_signal_corr=False,
        region_signal_file=atac_mat_file,
        rna_signal_file=rna_mat_file)

    # TODO
    # then, use the gene set files and overlap with rna clusters
    # to get counts per cluster
    
    return


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
    # ANALYSIS - get correlation matrix between ATAC and RNA
    # input: atac clusters/data, rna clusters/data
    # output: correlation matrix
    # -------------------------------------------
    logger.info("ANALYSIS: correlation between ATAC/RNA clusters")

    # set up dir
    corr_dir = "{}/traj.corr".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(corr_dir))

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
        corr_dir, prefix)
    if not os.path.isfile(correlation_matrix_file):
        build_correlation_matrix(
            atac_clusters_file,
            atac_mat_file,
            rna_clusters_file,
            rna_mat_file,
            correlation_matrix_file)
    
    # -------------------------------------------
    # ANALYSIS - naive proximity linking
    # input: ATAC regions and TSS file
    # output: linking by naive proximity
    # -------------------------------------------
    logger.info("ANALYSIS: naive proximity linking")
    days = ["d00", "d30", "d60"]
    
    # naive proximity
    proximity_links_dir = "{}/proximity".format(results_dir)
    if not os.path.isdir(proximity_links_dir):
        run_shell_cmd("mkdir -p {}".format(proximity_links_dir))
        method_type = "proximity"
        k_nearest = args.inputs["params"]["linking"][method_type]["k_nearest"]
        max_dist = args.inputs["params"]["linking"][method_type]["max_dist"]
        links_files = []
        for day in days:
            day_links_file = "{}/links.proximity.{}.k-{}.max_d-{}.txt.gz".format(
                proximity_links_dir, day, k_nearest, max_dist)
            if not os.path.isfile(day_links_file):
                # get tss file
                day_tss_file = args.outputs["data"]["tss.expressed_{}".format(day)]

                # get region file that has master region IDs (intersect with master file)
                day_orig_bed_file = glob.glob("{}/*{}*.gz".format(
                    args.outputs["results"]["atac"]["timepoint_region_dir"],
                    day))[0]
                day_regions_file = "{}/regions.{}.bed.gz".format(proximity_links_dir, day)
                intersect_cmd = (
                    "bedtools intersect -u -a {} -b {} | "
                    "gzip -c > {}").format(
                        args.outputs["data"]["atac.master.bed"],
                        day_orig_bed_file,
                        day_regions_file)
                run_shell_cmd(intersect_cmd)

                # set up naive proximity links
                link_by_distance(
                    day_regions_file,
                    day_tss_file,
                    day_links_file,
                    k_nearest=k_nearest,
                    max_dist=max_dist,
                    annotate_gene_ids=False)
            links_files.append(day_links_file)

        # and reconcile across timepoints
        all_prefix = "{}/{}.ALL".format(proximity_links_dir, prefix)
        all_interactions_file = "{}.interactions.txt.gz".format(all_prefix)
        if not os.path.isfile(all_interactions_file):
            get_timepoint_consistent_links(links_files, all_prefix)

    # proximity with corr
    proximity_corr_links_dir = "{}/proximity.corr".format(results_dir)
    if not os.path.isdir(proximity_corr_links_dir):
        run_shell_cmd("mkdir -p {}".format(proximity_corr_links_dir))
        method_type = "proximity.corr"
        k_nearest = args.inputs["params"]["linking"][method_type]["k_nearest"]
        max_dist = args.inputs["params"]["linking"][method_type]["max_dist"]
        corr_coeff_thresh = args.inputs["params"]["linking"][method_type]["corr_coeff_thresh"]
        pval_thresh = args.inputs["params"]["linking"][method_type]["pval_thresh"]
        links_files = []
        for day in days:
            day_links_file = "{}/links.proximity.{}.k-{}.max_d-{}.corr-{}.txt.gz".format(
                proximity_corr_links_dir, day, k_nearest, max_dist, corr_coeff_thresh)
            if not os.path.isfile(day_links_file):
                # get inputs
                day_tss_file = args.outputs["data"]["tss.expressed_{}".format(day)]
                day_regions_file = "{}/regions.{}.bed.gz".format(proximity_links_dir, day)

                # set up proximity links w correlation filter
                link_by_distance_and_corr(
                    day_regions_file,
                    day_tss_file,
                    day_links_file,
                    args.outputs["data"]["atac.counts.pooled.rlog.mat"],
                    args.outputs["data"]["rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat"],
                    k_nearest=k_nearest,
                    max_dist=max_dist,
                    corr_coeff_thresh=corr_coeff_thresh,
                    pval_thresh=pval_thresh,
                    is_ggr=True)
            links_files.append(day_links_file)

        # and reconcile across timepoints
        all_prefix = "{}/{}.ALL".format(proximity_corr_links_dir, prefix)
        all_interactions_file = "{}.interactions.txt.gz".format(all_prefix)
        if not os.path.isfile(all_interactions_file):
            get_timepoint_consistent_links(links_files, all_prefix)

    # proximity with corr: ONLY dynamic
    proximity_corr_dynamic_links_dir = "{}/proximity.corr.dynamic".format(results_dir)
    if not os.path.isdir(proximity_corr_dynamic_links_dir):
        run_shell_cmd("mkdir -p {}".format(proximity_corr_dynamic_links_dir))
        method_type = "proximity.corr.dynamic"
        k_nearest = args.inputs["params"]["linking"][method_type]["k_nearest"]
        max_dist = args.inputs["params"]["linking"][method_type]["max_dist"]
        corr_coeff_thresh = args.inputs["params"]["linking"][method_type]["corr_coeff_thresh"]
        pval_thresh = args.inputs["params"]["linking"][method_type]["pval_thresh"]
        links_files = []
        for day in days:
            day_links_file = "{}/links.proximity.{}.k-{}.max_d-{}.corr-{}.txt.gz".format(
                proximity_corr_dynamic_links_dir, day, k_nearest, max_dist, corr_coeff_thresh)
            if not os.path.isfile(day_links_file):
                # get inputs
                day_tss_file = args.outputs["data"]["tss.expressed_{}".format(day)]
                day_regions_file = "{}/regions.{}.bed.gz".format(proximity_links_dir, day)

                # set up proximity links w correlation filter
                link_by_distance_and_corr(
                    day_regions_file,
                    day_tss_file,
                    day_links_file,
                    args.outputs["data"]["atac.counts.pooled.rlog.dynamic.mat"],
                    args.outputs["data"]["rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.mat"],
                    k_nearest=k_nearest,
                    max_dist=max_dist,
                    corr_coeff_thresh=corr_coeff_thresh,
                    pval_thresh=pval_thresh,
                    is_ggr=True)
            links_files.append(day_links_file)

        # and reconcile across timepoints
        all_prefix = "{}/{}.ALL".format(proximity_corr_dynamic_links_dir, prefix)
        all_interactions_file = "{}.interactions.txt.gz".format(all_prefix)
        if not os.path.isfile(all_interactions_file):
            get_timepoint_consistent_links(links_files, all_prefix)

    # -------------------------------------------
    # ANALYSIS - different styles of ABC linking
    # input: ABC links
    # output: replicate consistent link sets
    # -------------------------------------------

    # distance based ABC
    logger.info("ANALYSIS: build replicate consistent ABC links - distance-based")
    abc_dist_dir = "abc.distance"
    if not os.path.isdir("{}/{}".format(results_dir, abc_dist_dir)):
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="distance",
            out_dir=abc_dist_dir)

    # average hiC based ABC
    logger.info("ANALYSIS: build replicate consistent ABC links - cell type avg HiC")
    abc_avg_hic_dir = "abc.hic.celltype_avg"
    if not os.path.isdir("{}/{}".format(results_dir, abc_avg_hic_dir)):
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="hic.celltype_avg",
            out_dir=abc_avg_hic_dir)

    # TODO - for each of these link sets, get correlation matrix
    # between ATAC and RNA
    # input: atac clusters, rna clusters/data
    links_dirs = [
        #"abc.distance"
        "abc.hic.celltype_avg"
        #"proximity",
        #"proximity.corr",
        #"proximity.corr.dynamic"
    ]
    for links_dir in links_dirs:
        # run workflow for using links to get from ATAC traj to RNA traj        
        run_traj_linking_workflow(args, prefix, "{}/{}".format(results_dir, links_dir))

    quit()

        
    # TODO get rid of everything down here once clean
        
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

    return args
