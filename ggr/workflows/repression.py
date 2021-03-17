# workflow analyzing repression in skin differentiation


import os
import gzip
import glob
import logging

import numpy as np
import pandas as pd

from ggr.util.utils import run_shell_cmd

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.bioinformatics import run_homer
from ggr.analyses.bioinformatics import make_deeptools_heatmap
from ggr.analyses.bioinformatics import run_bioinformatics_on_bed
from ggr.analyses.epigenome import get_distances_to_nearest_region
from ggr.analyses.epigenome import get_promoter_atac
from ggr.analyses.linking import regions_to_genes_w_correlation_filtering


def _plot_profile_heatmaps_workflow(args, tss_file, plot_prefix):
    """plot histone heatmaps for a specific set of tss regions (ordered
    """
    # make a subsample file if necessary
    tss_data = pd.read_csv(tss_file, sep="\t")
    if tss_data.shape[0] > 1000:
        step_size = int(tss_data.shape[0] / 1000)
        keep_indices = np.arange(0, tss_data.shape[0], step=step_size)
        print tss_data.shape
        tss_data = tss_data.iloc[keep_indices]
        print tss_data.shape
        subsample_file = "{}.subsample.tmp.bed.gz".format(plot_prefix)
        tss_data.to_csv(subsample_file, sep="\t",
                        header=False, index=False, compression="gzip")
        tss_file = subsample_file

        # NOTE rowseps won't match if subsampling!!
        print "NOTE rowseps not adjusted!!"
        
    # plot ATAC-seq
    if not os.path.isfile("{}.atac.heatmap.pdf".format(plot_prefix)):
        r_cmd = "~/git/ggr-project/R/plot.tss.timeseries_heatmap.R {} {} {}.atac ATAC-seq".format(
            "{}.promoter_data.mat.txt.gz".format(plot_prefix),
            "{}.row_seps.txt.gz".format(plot_prefix),
            plot_prefix)
        print r_cmd
        os.system(r_cmd)
    
    # plot RNA-seq
    if not os.path.isfile("{}.rna.heatmap.pdf".format(plot_prefix)):
        r_cmd = "~/git/ggr-project/R/plot.tss.timeseries_heatmap.R {} {} {}.rna PAS-seq".format(
            "{}.promoter_rna_data.mat.txt.gz".format(plot_prefix),
            "{}.row_seps.txt.gz".format(plot_prefix),
            plot_prefix)
        print r_cmd
        os.system(r_cmd)
    
    # plot each histone mark
    histones = ["H3K27ac", "H3K4me1", "H3K27me3"]    
    histone_colors = args.inputs["chipseq"][args.cluster]["histones"]["ordered_deeptools_colors"]
    histone_r_colors = args.inputs["chipseq"][args.cluster]["histones"]["ordered_r_colors"]
    
    for histone_idx in range(len(histones)):
        histone = histones[histone_idx]
        histone_color = histone_colors[histone_idx]
        histone_r_color = histone_r_colors[histone_idx]
        histone_bigwigs = sorted(
            glob.glob("{}/{}".format(
                args.inputs["chipseq"][args.cluster]["data_dir"],
                args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_bigwig_glob"])))

        out_prefix = "{}.{}_overlap".format(plot_prefix, histone)
        plot_file = "{}.heatmap.profile.pdf".format(out_prefix)
    
        if not os.path.isfile(plot_file):
            make_deeptools_heatmap(
                tss_file,
                histone_bigwigs,
                out_prefix,
                sort=False,
                referencepoint="center",
                color=histone_color)

        # make profile heatmap in R, with strand oriented TSS
        row_sep_file = "{}.row_seps.txt.gz".format(plot_prefix)
        out_mat_file = "{}.point.mat.gz".format(out_prefix)
        out_r_file = "{}.replot.pdf".format(plot_file.split(".pdf")[0])
        if not os.path.isfile(out_r_file):
            replot = (
                "~/git/ggr-project/R/plot.profile_heatmaps.stranded.R {} {} {} {} {} "
                "1,100 101,200 201,300").format(
                    out_mat_file,
                    histone, 
                    row_sep_file,
                    out_r_file,
                    histone_r_color)
            print replot
            run_shell_cmd(replot)
        
        
    return args




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
    
    # linking file
    links_dir = "{}/linking/proximity".format(args.outputs["results"]["dir"])
    interactions_file = "{}/ggr.linking.ALL.overlap.interactions.txt.gz".format(links_dir)

    # background atac file
    atac_background_file = args.outputs["data"]["atac.master.bed"]
    
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
        
        # setup
        region_dir = "{}/{}".format(
            work_dir, os.path.basename(region_set).split(".bed")[0])
        if not os.path.isdir(region_dir):
            run_shell_cmd("mkdir -p {}".format(region_dir))

        # make sure it has ATAC oriented IDs
        tmp_file = "{}/{}".format(region_dir, os.path.basename(region_set))
        regions = pd.read_csv(region_set, sep="\t", header=None)
        regions[3] = regions[0] + ":" + regions[1].astype(str) + "-" + regions[2].astype(str)
        regions.to_csv(tmp_file, sep="\t", compression="gzip", header=False, index=False)
        
        # get genes
        gene_file = "{}/genes.txt.gz".format(region_dir)
        if not os.path.isfile(gene_file):
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
        gprofiler_output = "{}/genes.go_gprofiler.txt".format(region_dir)
        if not os.path.isfile(gprofiler_output):
            run_gprofiler(
                gene_file,
                args.outputs["data"]["rna.counts.pc.expressed.mat"],
                region_dir, ordered=True, header=True)

        # homer analysis
        homer_dir = "{}/homer".format(region_dir)
        if not os.path.isdir(homer_dir):
            run_shell_cmd("mkdir -p {}".format(homer_dir))
            filt_bed_file = "{}.linked_regions.bed.gz".format(
                tmp_file.split(".bed")[0])
            run_homer(
                filt_bed_file,
                atac_background_file,
                homer_dir)
        
        # clean up
        #run_shell_cmd("rm {}".format(tmp_file))
        
    # ---------------------------
    # ANALYSIS: H3K27me3, ATAC-centric
    # ---------------------------
    
    # work dir
    work_dir = "{}/H3K27me3".format(results_dir)
    if not os.path.isdir(work_dir):
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
        if not os.path.isfile(gene_file):
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
        gprofiler_output = "{}/genes.go_gprofiler.txt".format(region_dir)
        if not os.path.isfile(gprofiler_output):
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
        if not os.path.isdir(region_dir):
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
        if not os.path.isfile(gene_file):
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
        gprofiler_output = "{}/genes.go_gprofiler.txt".format(region_dir)
        if not os.path.isfile(gprofiler_output):
            run_gprofiler(
                gene_file,
                args.outputs["data"]["rna.counts.pc.expressed.mat"],
                region_dir, ordered=True, header=True)

    # ---------------------------
    # ANALYSIS: TSS-centric viewpoint
    # ---------------------------
    work_dir = "{}/tss-centric".format(results_dir)
    if not os.path.isdir(work_dir):
        run_shell_cmd("mkdir -p {}".format(work_dir))
    out_prefix = "{}/{}".format(work_dir, prefix)
    
    # first, look at distances to TSSs from epigenetic feature
    tss_file = args.outputs["annotations"]["tss.pc.bed"]

    # get nearest TSS and distances for ATAC
    atac_file = args.outputs["data"]["atac.master.bed"]
    atac_distances_file = "{}.distance_to_tss.ATAC.txt.gz".format(out_prefix)
    if not os.path.isfile(atac_distances_file):
        get_distances_to_nearest_region(
            atac_file, tss_file, atac_distances_file)

    # nearest for histones
    histones = ["H3K27ac", "H3K4me1", "H3K27me3"]    
    for histone in histones:
        # get nearest TSS and distances
        histone_file = args.outputs["data"]["{}.master.bed".format(histone)]
        histone_distances_file = "{}.distance_to_tss.{}.txt.gz".format(out_prefix, histone)
        if not os.path.isfile(histone_distances_file):
            get_distances_to_nearest_region(
                histone_file, tss_file, histone_distances_file)

    # plot all together
    distance_files = sorted(glob.glob("{}*distance_to_tss*txt.gz".format(out_prefix)))
    plot_file = "{}.pdf".format(out_prefix)
    if not os.path.isfile(plot_file):
        # notes: this plot shows pretty similar splits for distance for all histone
        # marks
        r_cmd = "~/git/ggr-project/R/plot.distance_to_tss.R {} {}".format(
            out_prefix, " ".join(distance_files))
        print r_cmd
        os.system(r_cmd)

    # dynamic genes - how many of them have H3K27me3 at the promoter?
    dynamic_tss_w_H3K27me3_prefix = "{}.tss.dynamic.H3K27me3_present".format(out_prefix)
    dynamic_tss_w_H3K27me3_sorted = "{}.feature_sorted.bed.gz".format(
        dynamic_tss_w_H3K27me3_prefix)
    #if not os.path.isfile(dynamic_tss_w_H3K27me3_sorted):
    if True:
        # slopbed on TSS file only downstream direction to capture downstream H3K27me3 signal
        dynamic_tss_w_H3K27me3_bed = "{}.bed.gz".format(dynamic_tss_w_H3K27me3_prefix)
        downstream_bp_ext = 1000
        bedtools_cmd = (
            "bedtools slop -s -i {0} -g {1} -l 0 -r {2} | " # increase DOWNSTREAM tss region
            "bedtools intersect -u -a stdin -b {3} | " # intersect
            "awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$2+1\"\t\"$4\"\t\"$5\"\t\"$6 }}' | " # return to point
            "gzip -c > {4}").format(
                args.outputs["data"]["tss.dynamic"],
                args.inputs["annot"][args.cluster]["chromsizes"],
                downstream_bp_ext,
                args.outputs["data"]["H3K27me3.master.bed"],
                dynamic_tss_w_H3K27me3_bed)
        print bedtools_cmd
        os.system(bedtools_cmd)

        # now overlap with ATAC data so that we can have promoter info
        get_promoter_atac(
            dynamic_tss_w_H3K27me3_bed,
            args.outputs["data"]["atac.master.bed"],
            dynamic_tss_w_H3K27me3_prefix,
            rna_mat_file=args.outputs["data"][
                "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat"],
            rna_cluster_file=args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"],
            atac_mat_file=args.outputs["data"]["atac.counts.pooled.rlog.mat"],
            atac_cluster_file=args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"])

        # debug: read in hgnc mapping
        if False:
            hgnc_file = "{}.hgnc_ids.mat.txt.gz".format(dynamic_tss_w_H3K27me3.split(".bed")[0])
            mappings = pd.read_csv(args.outputs["annotations"]["geneids.mappings.mat"], sep="\t")
            tss_data = tss_data.merge(mappings, left_on=3, right_on="ensembl_gene_id")
            tss_data.to_csv(hgnc_file, sep="\t", compression="gzip", header=False, index=False)

    # plot data
    args = _plot_profile_heatmaps_workflow(
        args, dynamic_tss_w_H3K27me3_sorted, dynamic_tss_w_H3K27me3_prefix)
        
    # functional enrichments
    dynamic_genes_w_H3K27me3 = "{}.genes.dynamic.H3K27me3_at_tss.bed.gz".format(out_prefix)
    if not os.path.isfile(dynamic_genes_w_H3K27me3):
        tss_data = pd.read_csv(dynamic_tss_w_H3K27me3, sep="\t", header=None)
        tss_data.to_csv(
            dynamic_genes_w_H3K27me3, columns=[3], header=False, index=False, compression="gzip")
    gprofiler_results = "{}.go_gprofiler.txt".format(dynamic_genes_w_H3K27me3)
    if not os.path.isfile(gprofiler_results):
        run_gprofiler(
            dynamic_genes_w_H3K27me3,
            args.outputs["data"]["rna.counts.pc.expressed.mat"],
            work_dir, ordered=False, header=False)


    quit()
    
    # stable genes, subset to those with H3K27me3
    stable_tss_w_H3K27me3 = "{}.tss.stable.H3K27me3_present.bed.gz".format(out_prefix)
    if not os.path.isfile(stable_tss_w_H3K27me3):
        bedtools_cmd = "bedtools intersect -u -a {} -b {} | gzip -c > {}".format(
            args.outputs["data"]["tss.stable"],
            args.outputs["data"]["H3K27me3.master.bed"],
            stable_tss_w_H3K27me3)
        print bedtools_cmd
        os.system(bedtools_cmd)
        
    args = _plot_profile_heatmaps_workflow(args, stable_tss_w_H3K27me3, "{}".format(out_prefix))

    # for non-expressed, figure out how many have H3K27me3 repression on top
    non_expr_tss_w_H3K27me3 = "{}.tss.non_expr.H3K27me3_present.bed.gz".format(out_prefix)
    if not os.path.isfile(non_expr_tss_w_H3K27me3):
        bedtools_cmd = "bedtools intersect -u -a {} -b {} | gzip -c > {}".format(
            args.outputs["data"]["tss.non_expr"],
            args.outputs["data"]["H3K27me3.master.bed"],
            non_expr_tss_w_H3K27me3)
        print bedtools_cmd
        os.system(bedtools_cmd)

    args = _plot_profile_heatmaps_workflow(args, non_expr_tss_w_H3K27me3, "{}".format(out_prefix))


    quit()

    # to consider: analysis on looping across H3K27me3 domains?

            
    return args
