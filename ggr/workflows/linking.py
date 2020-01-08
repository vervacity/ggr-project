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
from ggr.analyses.linking import build_confusion_matrix
from ggr.analyses.linking import get_union_links_signal_mat


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

    if "hichip" in link_key:
        add_mirrored_interactions = True
    else:
        add_mirrored_interactions = False
    
    # go through each day
    for link_day in link_days:
        print link_day
        pooled_file = glob.glob(
            "{}/*{}*pooled*.gz".format(link_dir, link_day))
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
            idr_thresh=args.inputs["params"]["linking"][out_dir]["idr_thresh"],
            add_mirrored_interactions=add_mirrored_interactions)

    # and reconcile across timepoints
    all_prefix = "{}/{}.ALL".format(replicated_link_dir, prefix)
    all_interactions_file = "{}.interactions.txt.gz".format(all_prefix)
    if not os.path.isfile(all_interactions_file):
        links_files = sorted(
            glob.glob("{}/*interactions.txt.gz".format(
                replicated_link_dir)))
        get_timepoint_consistent_links(
            links_files, all_prefix)

    # now, for PCA/QC purposes, build a matrix of signals (by rep) across unified links
    if False:
        mat_files = glob.glob("{}/*idr_filt*mat.txt.gz".format(replicated_link_dir))
        get_union_links_signal_mat(all_interactions_file, mat_files)
    
        
    return args


def run_traj_linking_workflow(args, prefix, links_dir):
    """link ATAC trajectories to RNA trajectories with linking file
    """
    logging.info("ANALYSIS: link ATAC traj to downstream genes/gene trajs")

    # setup
    traj_dir = "{}/traj.linking".format(links_dir)
    os.system("mkdir -p {}".format(traj_dir))

    # ATAC traj files
    atac_clusters_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    atac_mat_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.traj.mat"]

    # RNA traj files
    rna_clusters_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
        "clusters.reproducible.hard.reordered.list"]
    rna_mat_file = args.outputs["data"][
        "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
    
    # get linking file
    interactions_file = "{}/{}.ALL.overlap.interactions.txt.gz".format(links_dir, prefix)

    # tss file and all expressed genes file
    tss_file = args.outputs["annotations"]["tss.pc.bed"]
    rna_expressed_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    
    # first, per cluster, get gene sets and enrichments
    # DONT run a correlation filter since we're checking region to gene mapping without
    # assuming anything about genes (outside of dynamic status)
    if "proximity" in links_dir:
        score_filter = 0.5
    elif "hichip" in links_dir:
        score_filter = None
    else:
        score_filter = 2

    if "hichip" in links_dir:
        extend_len = 2500
    else:
        extend_len = 0
        
    region_clusters_to_genes(
        atac_clusters_file,
        interactions_file,
        tss_file,
        traj_dir,
        region_signal_file=atac_mat_file,
        rna_signal_file=rna_mat_file,
        filter_gene_file=rna_mat_file, # only use dynamic genes
        filter_by_score=score_filter,
        run_enrichments=True,
        background_gene_file=rna_expressed_file,
        chromsizes=args.inputs["annot"][args.cluster]["chromsizes"],
        extend_len=extend_len)

    # then, use the gene set files and overlap with rna clusters
    # to get counts per cluster
    confusion_prefix = "{}/{}.traj.atac_to_rna".format(traj_dir, prefix)
    build_confusion_matrix(
        atac_clusters_file, rna_clusters_file, traj_dir, confusion_prefix)
    
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
    all_prefix = "{}/{}.ALL".format(proximity_links_dir, prefix)
    all_interactions_file = "{}.overlap.interactions.txt.gz".format(all_prefix)
    out_results["links.proximity"] = all_interactions_file
    if not os.path.isfile(all_interactions_file):
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
        get_timepoint_consistent_links(links_files, all_prefix)

    # proximity with corr
    proximity_corr_links_dir = "{}/proximity.corr".format(results_dir)
    if not os.path.isdir(proximity_corr_links_dir):
        run_shell_cmd("mkdir -p {}".format(proximity_corr_links_dir))
        method_type = "proximity.corr"
        k_nearest = args.inputs["params"]["linking"][method_type]["k_nearest"]
        max_dist = args.inputs["params"]["linking"][method_type]["max_dist"]
        corr_coeff_thresh = args.inputs["params"]["linking"][method_type]["corr_coeff_thresh"]
        if isinstance(corr_coeff_thresh, basestring):
            corr_coeff_thresh = None
        abs_corr_coeff_thresh = args.inputs["params"]["linking"][method_type]["abs_corr_coeff_thresh"]
        if isinstance(abs_corr_coeff_thresh, basestring):
            abs_corr_coeff_thresh = None
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
                    abs_corr_coeff_thresh=abs_corr_coeff_thresh,
                    is_ggr=True)
            links_files.append(day_links_file)

        # and reconcile across timepoints
        all_prefix = "{}/{}.ALL".format(proximity_corr_links_dir, prefix)
        all_interactions_file = "{}.interactions.txt.gz".format(all_prefix)
        if not os.path.isfile(all_interactions_file):
            get_timepoint_consistent_links(links_files, all_prefix)

    # -------------------------------------------
    # ANALYSIS - different styles of ABC linking
    # input: ABC links
    # output: replicate consistent link sets
    # -------------------------------------------

    # distance based ABC w TSS focus
    logger.info("ANALYSIS: build replicate consistent ABC links - distance-based, my provided TSS links")
    abc_dist_dir = "abc.distance.tss"
    if not os.path.isdir("{}/{}".format(results_dir, abc_dist_dir)):
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="distance.tss",
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

    # average hiC based ABC w TSS focus
    logger.info("ANALYSIS: build replicate consistent ABC links - cell type avg HiC, my provided tss links")
    abc_avg_hic_dir = "abc.hic.celltype_avg.tss"
    if not os.path.isdir("{}/{}".format(results_dir, abc_avg_hic_dir)):
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="hic.celltype_avg.tss",
            out_dir=abc_avg_hic_dir)

    # hichip
    logger.info("ANALYSIS: build replicate consistent ABC links - hichip")
    hichip_dir = "hichip"
    if not os.path.isdir("{}/{}".format(results_dir, hichip_dir)):
        run_replicate_consistent_links_workflow(
            args,
            prefix,
            results_dir,
            link_key="hichip",
            out_dir=hichip_dir)

        # get replicates into union links to be able to get PCA for QC 
        mat_files = sorted(
            glob.glob("{}/{}/*idr_filt*mat.txt.gz".format(results_dir, hichip_dir)))
        union_file = "{}/{}/ggr.linking.ALL.overlap.interactions.txt.gz".format(
            results_dir, hichip_dir)
        union_reps_file = "{}/{}/{}.ALL.reps.mat.txt.gz".format(results_dir, hichip_dir, prefix)
        get_union_links_signal_mat(
            all_interactions_file, mat_files, union_reps_file, "{}/{}".format(results_dir, hichip_dir))
        # PCA
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.pca.R {}/{}/{}.hichip.reps.pca.pdf {}".format(
            results_dir, hichip_dir, prefix, union_reps_file)
        print plot_cmd
        os.system(plot_cmd)
        
    # analyze trajectories (ATAC to RNA)
    links_dirs = [
        "proximity",
        "proximity.corr",
        "abc.distance.tss",
        "abc.hic.celltype_avg",
        "abc.hic.celltype_avg.tss",
        "hichip"
    ]
    
    for links_dir in links_dirs:
        # run workflow for using links to get from ATAC traj to RNA traj
        traj_dir = "{}/{}/traj.linking".format(results_dir, links_dir)
        if not os.path.isdir(traj_dir):
            run_traj_linking_workflow(args, prefix, "{}/{}".format(results_dir, links_dir))

    # viz: filter hichip (loose threshold) links with proximity links
    viz_dir = "{}/viz".format(results_dir)
    if not os.path.isdir(viz_dir):
        os.system("mkdir -p {}".format(viz_dir))
        hichip_interactions = "{}/hichip/ggr.linking.ALL.overlap.interactions.txt.gz".format(results_dir)
        proximity_interactions = "{}/proximity/ggr.linking.ALL.overlap.interactions.txt.gz".format(results_dir)
        out_file = "{}/ggr.linking.proximity.hichip_overlap.interactions.txt.gz".format(viz_dir)
        overlap_cmd = (
            "bedtools slop -i {} -g {} -b {} | "
            "bedtools intersect -u -a {} -b - | "
            "sort -k1,1 -k2,2n | "
            "bgzip -c > {}").format(
                hichip_interactions,
                args.inputs["annot"][args.cluster]["chromsizes"],
                2500,
                proximity_interactions,
                out_file)
        print overlap_cmd
        os.system(overlap_cmd)
        # index
        tabix_cmd = "tabix -p bed {}".format(out_file)
        print tabix_cmd
        os.system(tabix_cmd)
    
    return args
