# workflows integrating epigenomic datasets

import os
import gzip
import glob
import logging

import numpy as np
import pandas as pd

from ggr.util.utils import run_shell_cmd

from ggr.analyses.bioinformatics import make_deeptools_heatmap
from ggr.analyses.bioinformatics import run_bioinformatics_on_bed

from ggr.util.bed_utils import id_to_bed

from ggr.analyses.epigenome import get_histone_overlaps
from ggr.analyses.epigenome import cluster_by_chromatin_marks
from ggr.analyses.epigenome import split_stable_atac_by_dynamic_marks
from ggr.analyses.epigenome import get_best_summit

from ggr.analyses.filtering import sort_by_clusters
from ggr.analyses.filtering import get_ordered_subsample

from ggr.analyses.timeseries import plot_clusters


def count_region_nums(file_list, assay, out_file, method="genome_coverage"):
    """count regions and add to file
    """
    if method == "region_count":
        method = "wc -l"
    elif method == "genome_coverage":
        method = "awk -F '\t' '{{ sum += $3-$2 }} END {{ print sum }}'"
    elif method == "genome_fraction":
        method = "awk -F '\t' '{{ sum += $3-$2 }} END {{ print sum/3095693983.0  }}'"
    
    for filename in file_list:

        if "ATAC" in assay:
            awk_cmd = ("awk '{{ print \"{0}\t{1}\t\"$1 }}'").format(
                float(os.path.basename(filename).split(".")[0].split("-")[1].split("d")[1]) / 10,
                assay)
        else:
            awk_cmd = "awk '{{ print \"{0}\t{1}\t\"$1 }}'".format(
                float(os.path.basename(filename).split(".")[0].split("-")[1].split("d")[1]),
                assay)
        
        get_nums = (
            "zcat {0} | "
            "sort -k1,1 -k2,2n | "
            "bedtools merge -i stdin | "
            "{1} | "
            "{2} >> "
            "{3} ").format(
                filename,
                method,
                awk_cmd,
                out_file)
        print get_nums
        os.system(get_nums)

    return None


def get_epigenome_static_metrics_workflow(args, prefix):
    """Per timepoint file (ATAC), determine global stats
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get static metrics across GGR epigenome")

    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "static"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]
    
    # set up output file
    region_nums_key = "epigenome.region_nums"
    out_results[region_nums_key] = "{0}/{1}.region_nums_per_timepoint.txt".format(
        results_dir, prefix)

    if not os.path.isfile(out_results[region_nums_key]):
        # header
        print_header = "echo \"timepoint\ttype\tcount\" > {0}".format(
            out_results[region_nums_key])
        run_shell_cmd(print_header)
        
        # get total counts for ATAC
        atac_timepoint_files = sorted(
            glob.glob("{}/*.narrowPeak.gz".format(
                args.outputs["results"]["atac"]["timepoint_region_dir"])))
        for timepoint_string in args.inputs["params"]["media_timepoints"]:
            atac_timepoint_files = [
                filename for filename in atac_timepoint_files
                if timepoint_string not in filename]
        count_region_nums(
            atac_timepoint_files,
            "ATAC-seq",
            out_results[region_nums_key])
        
        # get counts for histones
        histone_inputs = args.inputs["chipseq"][args.cluster]
        histones = histone_inputs["histones"]["ordered_names"]
        
        for histone in histones:
            histone_overlap_files = sorted(
                glob.glob("{0}/{1}".format(
                    histone_inputs["data_dir"],
                    histone_inputs["histones"][histone]["overlap_glob"])))
            count_region_nums(
                histone_overlap_files,
                "{} ChIP-seq".format(histone),
                out_results[region_nums_key])
            
        # get counts for CTCF
        tf_inputs = args.inputs["chipseq"][args.cluster]
        tfs = tf_inputs["tfs"].keys()
        for tf in tfs:
            peak_files = sorted(
                glob.glob("{0}/{1}".format(
                    tf_inputs["data_dir"],
                    tf_inputs["tfs"][tf]["idr_peak_glob"])))
            count_region_nums(
                peak_files,
                "{} ChIP-seq".format(tf),
                out_results[region_nums_key])

        # and plot
        plot_file = "{}.pdf".format(out_results[region_nums_key].split(".txt")[0])
        plot_counts = "plot.region_nums.R {} {}".format(
            out_results[region_nums_key], plot_file)
        run_shell_cmd(plot_counts)
    
    return args

def get_timepoint_dynamics(deseq_files, region_id_file, assay, out_file):
    """Read in deseq files and calculate signal summary stats
    """
    region_ids = pd.read_table(region_id_file, header=None).iloc[:,0].tolist()

    with open(out_file, "a") as out:
        for deseq_file in deseq_files:
            
            # get prefix
            if "ATAC" in assay:
                prefix = os.path.basename(deseq_file).split(".")[2].split("_results")[0]
                prefix = list(prefix)
                prefix.insert(2, ".")
                prefix.insert(12, ".")
                prefix = "".join(prefix)
            elif "H3" in assay:
                prefix = os.path.basename(deseq_file).split(".")[3].split("_results")[0]
            else:
                prefix = os.path.basename(deseq_file).split(".")[2].split("_results")[0]

            # adjust prefix
            prefix = "{}_to_{}".format(
                prefix.split("_over_")[1],
                prefix.split("_over_")[0])

            # get deseq info
            deseq_results = pd.read_table(deseq_file, index_col=0)
            
            # filter
            deseq_results = deseq_results.loc[region_ids,:]
        
            # calculate. convert into normal signal info, then get fold change, and then return to log2 space
            pos_changes = deseq_results[deseq_results["log2FoldChange"] > 0.0]
            final_signal = np.sum(pos_changes["baseMean"] * 2**pos_changes["log2FoldChange"])
            initial_signal = np.sum(pos_changes["baseMean"])
            pos_FC = final_signal / initial_signal
            posLog2FC = np.log2(pos_FC)
            
            neg_changes = deseq_results[deseq_results["log2FoldChange"] < 0.0]
            final_signal = np.sum(neg_changes["baseMean"] * 2**neg_changes["log2FoldChange"])
            initial_signal = np.sum(neg_changes["baseMean"])
            neg_FC = final_signal / initial_signal
            negLog2FC = np.log2(neg_FC)
            
            # save out
            out.write("{}\t{}\t{}\t{}\n".format(prefix, assay, posLog2FC, negLog2FC))
    
    return None



def get_epigenome_dynamic_metrics_workflow(args, prefix):
    """Get the differential positive and negative signals for 
    each epigenomic dataset
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get static metrics across GGR epigenome")

    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "dynamic"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    # TODO figure out if this is ok
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]
    
    # set up output file
    region_nums_key = "epigenome.timepoint_dynamics"
    out_results[region_nums_key] = "{0}/{1}.region_dynamics_per_timepoint.txt".format(
        results_dir, prefix)

    if not os.path.isfile(out_results[region_nums_key]):
        # header
        print_header = "echo \"timepoint\ttype\tposLog2FC\tnegLog2FC\" > {0}".format(
            out_results[region_nums_key])
        run_shell_cmd(print_header)

        # get dynamic regions handle
        dynamic_regions_file = args.outputs["results"]["atac"]["timeseries"]["dynamic_ids.list"]
        
        # get ATAC deseq files and extract signal
        days = args.inputs["params"]["hi_res_filt_days"]
        atac_deseq2_files_all = sorted(
            glob.glob("{}/deseq2/*resultsAll.txt.gz".format(
                args.outputs["results"]["atac"]["timeseries"]["dir"])))
        timepoint_comparisons = [
            "{}_over_{}".format(days[i+1], days[i])
            for i in xrange(len(days)-1)]
        atac_deseq2_files = []
        for timepoint_string in timepoint_comparisons:
            atac_deseq2_files += [
                filename for filename in atac_deseq2_files_all
                if timepoint_string in filename]
        get_timepoint_dynamics(
            atac_deseq2_files,
            dynamic_regions_file,
            "ATAC-seq",
            out_results[region_nums_key])

        # set up for lo res groups
        days = args.inputs["params"]["lo_res_days"]
        timepoint_comparisons = [
            "{}_over_{}".format(days[i+1], days[i])
            for i in xrange(len(days)-1)]
        
        # get counts for histones
        histone_inputs = args.inputs["chipseq"][args.cluster]
        histones = histone_inputs["histones"]["ordered_names"]
        
        for histone in histones:

            # get dynamic regions handle
            dynamic_regions_file = args.outputs["results"]["histones"][histone]["timeseries"]["ggr.histone.{}.dynamic.ids.list".format(histone)]

            # get histone deseq files
            histone_deseq2_files_all = sorted(
            glob.glob("{}/deseq2/*resultsAll.txt.gz".format(
                args.outputs["results"]["histones"][histone]["timeseries"]["dir"])))
            histone_deseq2_files = []
            for timepoint_string in timepoint_comparisons:
                histone_deseq2_files += [
                    filename for filename in histone_deseq2_files_all
                    if timepoint_string in filename]
            get_timepoint_dynamics(
                histone_deseq2_files,
                dynamic_regions_file,
                "{} ChIP-seq".format(histone),
                out_results[region_nums_key])
            
        # get counts for CTCF
        tf_inputs = args.inputs["chipseq"][args.cluster]
        tfs = tf_inputs["tfs"].keys()

        # for now don't analyze, need to run deseq on tfs
        tfs = []
        
        for tf in tfs:

            # get dynamic regions handle

            # get tf deseq files
            pass

        # and plot
        plot_file = "{}.pdf".format(out_results[region_nums_key].split(".txt")[0])
        plot_counts = "plot.region_dynamics.R {} {}".format(
            out_results[region_nums_key], plot_file)
        run_shell_cmd(plot_counts)

    return args


def make_fake_cluster_file(id_list_file, out_file):
    """Make a fake cluster file
    """
    with gzip.open(id_list_file, "r") as fp:
        with gzip.open(out_file, "w") as out:
            out.write("cluster\tid\n")
            for line in fp:
                if "d00" in line:
                    continue
                out.write("0\t{}\n".format(line.strip().split("\t")[0]))
                
    return None


def run_dynamic_epigenome_workflow(
        args,
        prefix,
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
    logger.info("ANALYSIS: overlap ATAC with histones")
    atac_dynamic_bed_key = "atac.dynamic.bed"
    histone_overlap_dir = "{}/overlap_histone".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(histone_overlap_dir))

    histone_overlap_mat_key = "{}.overlap_histones.mat".format(atac_dynamic_bed_key.split(".bed")[0])
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

    # TODO factor this out to subworkflow?
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
        (atac_cluster_file, [atac_cluster_col]),
        (histone_cluster_file, [histone_cluster_col])]

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

    # produce a subsample BED too
    atac_ordered_subsample_bed_key = "{}.bed".format(
        atac_ordered_subsample_key.split(".list")[0])
    out_results[atac_ordered_subsample_bed_key] = "{}.bed".format(
        out_results[atac_ordered_subsample_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_bed_key]):
        id_to_bed(
            out_results[atac_ordered_subsample_key],
            "{}.gz".format(out_results[atac_ordered_subsample_bed_key]),
            col=2)

        # and unzip
        run_shell_cmd("gunzip {}.gz".format(out_results[atac_ordered_subsample_bed_key]))

    # -----------------------------------------
    # ANALYSIS 2 - separate out cluster files (ids and beds)
    # inputs: cluster file
    # outputs: separated id and bed files
    # -----------------------------------------
    logger.info("ANALYSIS: separate out cluster files")
    mark_dir = "{}/clusters/by_mark".format(results_dir)
    state_dir = "{}/clusters/by_state".format(results_dir)
    if not os.path.isdir(state_dir):
        run_shell_cmd("mkdir -p {0}/ids {0}/bed {1}/ids {1}/bed".format(
            mark_dir, state_dir))
        cluster_by_chromatin_marks(
            out_results[atac_ordered_clusters_key],
            ["H3K27ac", "H3K4me1"],
            "histone_cluster",
            "{}/ids".format(mark_dir),
            "{}/ids".format(state_dir),
            prefix,
            min_region_num=args.inputs["params"]["chrom_state_cluster_min_size"])

        # and now convert files to beds and keep bed dir
        mark_files = glob.glob("{}/ids/*.txt.gz".format(mark_dir))
        for mark_file in mark_files:
            mark_bed_file = "{0}/bed/{1}.bed.gz".format(
                mark_dir, os.path.basename(mark_file).split('.txt')[0])
            id_to_bed(mark_file, mark_bed_file, sort=True)
        
        state_files = glob.glob("{}/ids/*.txt.gz".format(state_dir))
        for state_file in state_files:
            state_bed_file = "{0}/bed/{1}.bed.gz".format(
                state_dir, os.path.basename(state_file).split('.txt')[0])
            id_to_bed(state_file, state_bed_file, sort=True)

    mark_bed_dir = "{}/clusters/by_mark/bed".format(results_dir)
    state_bed_dir = "{}/clusters/by_state/bed".format(results_dir)
            
    args.outputs["results"]["label_dirs"].append(mark_bed_dir)
    args.outputs["results"]["label_dirs"].append(state_bed_dir)
        
    # -----------------------------------------
    # ANALYSIS 3 - plot this out
    # inputs: dynamic ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: plot out the heatmaps with the ordering in the subsample")
    mat_key = "atac.counts.pooled.rlog.dynamic.mat"
    plot_dir = "{}/plots".format(results_dir)

    # TODO - want to generate a file that has the strongest ATAC
    # summit across time, for each region. use that summit
    # to center on for the histone marks.
    # make this an analysis
    # build this on the ordered subsample
    ordered_subsample_summits_key = "{}.summits".format(atac_ordered_subsample_bed_key)
    out_results[ordered_subsample_summits_key] = "{}.summits.bed".format(
        out_results[atac_ordered_subsample_bed_key].split(".bed")[0])
    if not os.path.isfile(out_results[ordered_subsample_summits_key]):

        # get atac timeseries files
        timeseries_dir = args.outputs["results"]["atac"]["timepoint_region_dir"]
        atac_timeseries_files = sorted(
            glob.glob("{}/*narrowPeak.gz".format(timeseries_dir)))
        
        # get best summits
        get_best_summit(
            out_results[atac_ordered_subsample_bed_key],
            atac_timeseries_files,
            out_results[ordered_subsample_summits_key])
    
    # plot ATAC with the summits
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
        histone_r_colors = args.inputs["chipseq"][args.cluster]["histones"]["ordered_r_colors"]
        for histone_idx in range(len(histones)):
            histone = histones[histone_idx]
            histone_color = histone_colors[histone_idx]
            histone_r_color = histone_r_colors[histone_idx]
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
                    out_results[ordered_subsample_summits_key],
                    histone_bigwigs,
                    out_prefix,
                    sort=False,
                    referencepoint="center",
                    extend_dist=1500,
                    bin_total=100,
                    color=histone_color)

            # and then make own heatmap file in R with matrix output
            row_sep_file = "{}/{}.row_seps.txt".format(plot_dir, prefix)
            out_mat_file = "{}.point.mat.gz".format(out_prefix)
            out_r_file = "{}.replot.pdf".format(out_file.split(".pdf")[0])
            if not os.path.isfile(out_r_file):
                replot = (
                    "plot.profile_heatmaps.R {0} {1} {2} {3} "
                    "1,100 101,200 201,300").format(
                        out_mat_file,
                        row_sep_file,
                        out_r_file,
                        histone_r_color)
                run_shell_cmd(replot)
            
    # -----------------------------------------
    # ANALYSIS 4 - bioinformatics
    # inputs: BED dirs
    # outputs: HOMER/GREAT results
    # -----------------------------------------
    logger.info("ANALYSIS: run HOMER and GREAT")
    bioinformatics_bed_dirs = [mark_bed_dir, state_bed_dir]
    for bioinformatics_bed_dir in bioinformatics_bed_dirs:

        if not os.path.isdir("{}/homer".format(bioinformatics_bed_dir)):
            
            bed_files = glob.glob("{}/*.bed.gz".format(bioinformatics_bed_dir))
            
            # make a background bed file: stable + this info?
            background_bed_file = "{}/{}.atac.background.bed.gz".format(
                bioinformatics_bed_dir, prefix)
            make_background = (
                "zcat {0} {1} | "
                "sort -k1,1 -k2,2n | "
                "gzip -c > {2}").format(
                    " ".join(bed_files),
                    out_data["atac.stable.bed"],
                    background_bed_file)
            run_shell_cmd(make_background)
            
            # run files
            for bed_file in bed_files:
                run_bioinformatics_on_bed(
                    bed_file,
                    background_bed_file,
                    bioinformatics_bed_dir)

    return args


def run_stable_epigenome_workflow(
        args,
        prefix,
        histone_assignment="nearest"):
    """Same thing as above, except all ATAC is one big cluster (no splits)
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get stable epigenome")

    # assertions
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "stable"
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
    logger.info("ANALYSIS: overlap ATAC with histones")
    atac_stable_bed_key = "atac.stable.bed"
    histone_overlap_dir = "{}/overlap_histone".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(histone_overlap_dir))

    histone_overlap_mat_key = "{}.overlap_histones.mat".format(atac_stable_bed_key.split(".bed")[0])
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
            out_data[atac_stable_bed_key],
            histone_files,
            histone_overlap_dir,
            out_results[histone_overlap_mat_key],
            args.inputs["annot"][args.cluster]["chromsizes"],
            histone_assignment=histone_assignment)

    # -----------------------------------------
    # ANALYSIS 1 - reorder ATAC with all clusters
    # inputs: stable ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: reorder ATAC regions using histone clusters")

    # make a fake ATAC clusters file
    atac_stable_mat_key = "atac.counts.pooled.rlog.stable.mat"
    atac_fake_clusters_key = "atac.stable.fake_cluster.list"
    out_results[atac_fake_clusters_key] = "{}/{}.atac.clusters.txt.gz".format(results_dir, prefix)
    if not os.path.isfile(out_results[atac_fake_clusters_key]):
        make_fake_cluster_file(
            out_data[atac_stable_mat_key],
            out_results[atac_fake_clusters_key])
    atac_cluster_cols = ["cluster"]
        
    histone_clusters_key = histone_overlap_mat_key
    histone_cluster_file = out_results[histone_clusters_key]
    histone_cluster_cols = ["histone_cluster"]

    cluster_files = [
        (out_results[atac_fake_clusters_key], atac_cluster_cols),
        (histone_cluster_file, histone_cluster_cols)]

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
    # ANALYSIS 2 - separate out cluster files (ids and beds)
    # inputs: cluster file
    # outputs: separated id and bed files
    # -----------------------------------------
    logger.info("ANALYSIS: separate out cluster files")
    mark_dir = "{}/clusters/by_mark".format(results_dir)
    state_dir = "{}/clusters/by_state".format(results_dir)
    mark_bed_dir = "{}/clusters/by_mark/bed".format(results_dir)
    state_bed_dir = "{}/clusters/by_state/bed".format(results_dir)
    if not os.path.isdir(state_dir):
        run_shell_cmd("mkdir -p {0}/ids {0}/bed {1}/ids {1}/bed".format(
            mark_dir, state_dir))
        cluster_by_chromatin_marks(
            out_results[atac_ordered_clusters_key],
            ["H3K27ac", "H3K4me1", "H3K27me3"],
            "histone_cluster",
            "{}/ids".format(mark_dir),
            "{}/ids".format(state_dir),
            prefix,
            min_region_num=args.inputs["params"]["chrom_state_cluster_min_size"])

        # and now convert files to beds
        mark_files = glob.glob("{}/ids/*.txt.gz".format(mark_dir))
        for mark_file in mark_files:
            mark_bed_file = "{0}/bed/{1}.bed.gz".format(
                mark_dir, os.path.basename(mark_file).split('.txt')[0])
            id_to_bed(mark_file, mark_bed_file, sort=True)
        
        state_files = glob.glob("{}/ids/*.txt.gz".format(state_dir))
        for state_file in state_files:
            state_bed_file = "{0}/bed/{1}.bed.gz".format(
                state_dir, os.path.basename(state_file).split('.txt')[0])
            id_to_bed(state_file, state_bed_file, sort=True)

    args.outputs["results"]["label_dirs"].append(mark_bed_dir)
    args.outputs["results"]["label_dirs"].append(state_bed_dir)

    # -----------------------------------------
    # ANALYSIS 3 - separate out dynamic and non-dynamic chromatin mark data,
    # mostly for plotting
    # inputs: clusters
    # outputs: separate id lists
    # -----------------------------------------
    cluster_file_key = atac_ordered_clusters_key
    stable_atac_dynamic_histones_key = "{}.dynamic_histones.clusters".format(
        cluster_file_key.split(".mat")[0])
    out_results[stable_atac_dynamic_histones_key] = "{}/{}.dynamic_histones.clusters.txt".format(
        results_dir, os.path.basename(out_results[cluster_file_key]).split(".clusters")[0])
    stable_atac_stable_histones_key = "{}.stable_histones.clusters".format(
        cluster_file_key.split(".mat")[0])
    out_results[stable_atac_stable_histones_key] = "{}/{}.stable_histones.clusters.txt".format(
        results_dir, os.path.basename(out_results[cluster_file_key]).split(".clusters")[0])
    if not os.path.isfile(out_results[stable_atac_stable_histones_key]):
        # TODO: check that it does not reorder
        split_stable_atac_by_dynamic_marks(
            out_results[cluster_file_key],
            out_results[stable_atac_dynamic_histones_key],
            out_results[stable_atac_stable_histones_key])
        
    # keep keys for future use
    stable_set_keys = [
        stable_atac_dynamic_histones_key,
        stable_atac_stable_histones_key]

    # get subsamples, and subsample bed files for each of these
    subsample_bed_keys = []
    for stable_set_key in stable_set_keys:
        # and generate a subsample for plotting too
        subsample_mat_key = "{}.subsample.mat".format(
            stable_set_key.split(".list")[0])
        out_results[subsample_mat_key] = "{}.subsampled.mat".format(
            out_results[stable_set_key].split(".txt")[0])
        if not os.path.isfile(out_results[subsample_mat_key]):
            get_ordered_subsample(
                out_results[stable_set_key],
                out_results[subsample_mat_key])

        # extract just the id
        subsample_key = "{}.subsample.list".format(
            stable_set_key.split(".list")[0])
        out_results[subsample_key] = "{}.subsampled.txt".format(
            out_results[stable_set_key].split(".txt")[0])
        if not os.path.isfile(out_results[subsample_key]):
            get_ids = "cat {} | awk '{{print $2}}' | grep -v id > {}".format(
                out_results[subsample_mat_key],
                out_results[subsample_key])
            run_shell_cmd(get_ids)
        
        # produce a BED too
        subsample_bed_key = "{}.bed".format(
            subsample_key.split(".list")[0])
        out_results[subsample_bed_key] = "{}.bed".format(
            out_results[subsample_key].split(".txt")[0])
        subsample_bed_keys.append(subsample_bed_key)
        if not os.path.isfile(out_results[subsample_bed_key]):
            id_to_bed(
                out_results[subsample_key],
                "{}.gz".format(out_results[subsample_bed_key]))

            # and unzip
            run_shell_cmd("gunzip {}.gz".format(out_results[subsample_bed_key]))
        
    # -----------------------------------------
    # ANALYSIS 4 - plot this out
    # inputs: stable groups
    # outputs: heatmaps with histone marks
    # -----------------------------------------
    logger.info("ANALYSIS: plot out the heatmaps with the ordering in the subsamples")
    plot_dir = "{}/plots".format(results_dir)
    
    # plot ATAC and histone marks
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))

        for subsample_bed_key in subsample_bed_keys:

            # plot the ATAC alongside, using summits
            ordered_subsample_summits_key = "{}.summits".format(subsample_bed_key)
            out_results[ordered_subsample_summits_key] = "{}.summits.bed".format(
                out_results[subsample_bed_key].split(".bed")[0])
            if not os.path.isfile(out_results[ordered_subsample_summits_key]):
                
                # get atac timeseries files
                timeseries_dir = args.outputs["results"]["atac"]["timepoint_region_dir"]
                atac_timeseries_files = sorted(
                    glob.glob("{}/*narrowPeak.gz".format(timeseries_dir)))
        
                # get best summits
                get_best_summit(
                    out_results[subsample_bed_key],
                    atac_timeseries_files,
                    out_results[ordered_subsample_summits_key])
    
            # plot ATAC with the summits
            if False:
                plot_clusters(
                    out_results[atac_fake_clusters_key],
                    out_results[atac_ordered_subsample_key],
                    out_data[atac_stable_mat_key],
                    plot_dir,
                    prefix)
            
            # TODO need to plot color bars?
            
            # plot the histone signal profiles with deeptools
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

                out_prefix = "{}/{}.{}_overlap.{}".format(
                    plot_dir,
                    prefix,
                    histone,
                    subsample_bed_key.split(".")[-4])
                out_file = "{}.heatmap.profile.pdf".format(out_prefix)
                
                if not os.path.isfile(out_file):
                    make_deeptools_heatmap(
                        out_results[subsample_bed_key],
                        histone_bigwigs,
                        out_prefix,
                        sort=False,
                        referencepoint="center",
                        color=histone_color)

                
                # and then make own heatmap file in R with matrix output
                row_sep_file = "{}/{}.row_seps.txt".format(plot_dir, prefix)
                out_mat_file = "{}.point.mat.gz".format(out_prefix)
                out_r_file = "{}.replot.pdf".format(out_file.split(".pdf")[0])
                if not os.path.isfile(out_r_file):
                    replot = (
                        "plot.profile_heatmaps.R {0} {1} {2} {3} "
                        "1,100 101,200 201,300").format(
                            out_mat_file,
                            row_sep_file,
                            out_r_file,
                            histone_r_color)
                    run_shell_cmd(replot)
<<<<<<< HEAD

=======
>>>>>>> afae7595e389c5624099d49214e13134d397ec2b
    
    # -----------------------------------------
    # ANALYSIS 5 - bioinformatics
    # inputs: BED dirs
    # outputs: HOMER/GREAT results
    # -----------------------------------------
    logger.info("ANALYSIS: run HOMER and GREAT")
    bioinformatics_bed_dirs = [mark_bed_dir, state_bed_dir]
    for bioinformatics_bed_dir in bioinformatics_bed_dirs:

        if not os.path.isdir("{}/homer".format(bioinformatics_bed_dir)):
            
            bed_files = glob.glob("{}/*.bed.gz".format(bioinformatics_bed_dir))
            
            # make a background bed file: stable + this info?
            background_bed_file = "{}/{}.atac.background.bed.gz".format(
                bioinformatics_bed_dir, prefix)
            make_background = (
                "zcat {0} {1} | "
                "sort -k1,1 -k2,2n | "
                "gzip -c > {2}").format(
                    " ".join(bed_files),
                    out_data["atac.stable.bed"],
                    background_bed_file)
            run_shell_cmd(make_background)
            
            # run files
            for bed_file in bed_files:
                run_bioinformatics_on_bed(
                    bed_file,
                    background_bed_file,
                    bioinformatics_bed_dir)

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
    # ANALYSIS - per timepoint, total num of regions and overlap with histones
    # inputs - peak files
    # outputs - plot of region counts
    # -----------------------------------------
    logger.info("ANALYSIS: get region counts for epigenome datasets")
    args = get_epigenome_static_metrics_workflow(
        args, "{}.static".format(prefix))

    # -----------------------------------------
    # ANALYSIS - per timepoint, dynamics of regions
    # inputs - deseq2 files
    # outputs - plot of region dynamics per timepoint
    # -----------------------------------------
    logger.info("ANALYSIS: get region dynamics for epigenome datasets")
    args = get_epigenome_dynamic_metrics_workflow(
        args, "{}.dynamic".format(prefix))
    
    # -----------------------------------------
    # ANALYSIS 0 - integrate histones into dynamic ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    logger.info("ANALYSIS: integrate dynamic ATAC w histones")
    args = run_dynamic_epigenome_workflow(
        args, "{}.dynamic".format(prefix))

    # -----------------------------------------
    # ANALYSIS 1 - integrate histones into stable ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    logger.info("ANALYSIS: integrate stable ATAC w histones")
    args = run_stable_epigenome_workflow(
        args, "{}.stable".format(prefix))

    # TODO: look at the histone marks that are outside of ATAC regions
    

    return args
