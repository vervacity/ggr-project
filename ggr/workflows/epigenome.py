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
from ggr.analyses.epigenome import convert_overlaps_bed_to_id_mappings
from ggr.analyses.epigenome import get_aggregate_chromatin_state_summary
from ggr.analyses.epigenome import plot_region_set_chromatin_summary
from ggr.analyses.epigenome import plot_signal_aggregation_plots

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
        #tfs = tf_inputs["tfs"].keys()
        tfs = ["CTCF"]
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


def _run_region_set_epigenome_plotting(
        args, id_file, plot_dir, signal_matrices, convert_beds, atac_timeseries_files):
    """given region ids, plot out dynamics and agg plots
    """
    plot_prefix = "{}/{}".format(
        plot_dir, os.path.basename(id_file).split(".txt")[0])

    # plot average dynamics across timepoints
    plot_file = "{}.avg_dynamics.pdf".format(plot_prefix)
    if not os.path.isfile(plot_file):
        plot_region_set_chromatin_summary(
            id_file,
            plot_prefix,
            convert_beds,
            signal_matrices)
            
    # TODO - deeptools to get agg plots to check for spreading
    # plot agg plots: ATAC, H3K27ac, H3K4me1
    id_summit_file = "{}/{}.summit.bed".format(
        plot_dir, os.path.basename(id_file).split(".txt")[0])
    if not os.path.isfile(id_summit_file):
        get_best_summit(
            id_file.replace("/ids", "/bed").replace(".txt.gz", ".bed.gz"),
            atac_timeseries_files,
            id_summit_file)
        
    plot_prefix = "{}.agg".format(plot_prefix)

    # plot ATAC agg - look for spreading, landscape
    plot_atac_prefix = "{}.ATAC".format(plot_prefix)
    if not os.path.isfile("{}.agg_plot.pdf".format(plot_atac_prefix)):
        atac_bigwigs = sorted(glob.glob("{}/{}".format(
            args.inputs["atac"][args.cluster]["data_dir"],
            args.inputs["atac"][args.cluster]["bigwig_pooled_glob"])))
        plot_signal_aggregation_plots(
            id_summit_file, atac_bigwigs, plot_atac_prefix,
            extend_dist=2000, bin_total=200)

    # plot histone agg - look for spreading, landscape
    histones_agg = ["H3K27ac", "H3K4me1"]
    for histone_idx in range(len(histones_agg)):
        histone = histones_agg[histone_idx]
        histone_bigwigs = sorted(
            glob.glob("{}/{}".format(
                args.inputs["chipseq"][args.cluster]["data_dir"],
                args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_bigwig_glob"])))
        plot_histone_prefix = "{}.{}".format(plot_prefix, histone)
        if not os.path.isfile("{}.agg_plot.pdf".format(plot_histone_prefix)):
            plot_signal_aggregation_plots(
                id_summit_file, histone_bigwigs, plot_histone_prefix,
                extend_dist=2000, bin_total=200)
    
    return



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

        # TODO run H3K27me3 separately
        # so that don't get a cluster file but have the nearest neighbor mark
        inactive_mark_dir = "{}/overlap_histone.inactive".format(results_dir)
        run_shell_cmd("mkdir -p {}".format(inactive_mark_dir))
        inactive_overlap_mat = "{}/{}.overlaps.mat.txt.gz".format(inactive_mark_dir, prefix)
        get_histone_overlaps(
            out_data[atac_dynamic_bed_key],
            [histone_files[2]],
            inactive_mark_dir,
            inactive_overlap_mat,
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

    # for plotting
    
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

    # for plotting dynamics plots
    signal_matrices = [
        out_data["atac.counts.pooled.rlog.mat"],
        out_data["H3K27ac.counts.pooled.rlog.mat"],
        out_data["H3K4me1.counts.pooled.rlog.mat"]
    ]
    convert_beds = [
        (out_data["atac.master.bed"], out_data["atac.master.bed"]),
        (out_data["atac.master.bed"], out_data["atac.midpoints.slop_H3K27ac.bed"]),
        (out_data["atac.master.bed"], out_data["atac.midpoints.slop_H3K4me1.bed"]),
    ]

    # for plotting agg plots
    timeseries_dir = args.outputs["results"]["atac"]["timepoint_region_dir"]
    atac_timeseries_files = sorted(
        glob.glob("{}/*narrowPeak.gz".format(timeseries_dir)))
    
    # summary plots
    mark_id_dir = "{}/clusters/by_mark/ids".format(results_dir)
    mark_plot_dir = "{}/clusters/by_mark/plots".format(results_dir)
    os.system("mkdir -p {}".format(mark_plot_dir))
    id_files = sorted(glob.glob("{}/*.txt.gz".format(mark_id_dir)))
    for id_file in id_files:

        # ignore non-changing histone info
        if "cluster_4.txt" in id_file:
            continue

        # run plotting
        _run_region_set_epigenome_plotting(
            args, id_file, mark_plot_dir, signal_matrices, convert_beds, atac_timeseries_files)

    state_id_dir = "{}/clusters/by_state/ids".format(results_dir)
    state_plot_dir = "{}/clusters/by_state/plots".format(results_dir)
    os.system("mkdir -p {}".format(state_plot_dir))
    id_files = sorted(glob.glob("{}/*.txt.gz".format(state_id_dir)))
    for id_file in id_files:

        # ignore non-changing histone info
        if "cluster_44.txt" in id_file:
            continue

        # run plotting
        _run_region_set_epigenome_plotting(
            args, id_file, state_plot_dir, signal_matrices, convert_beds, atac_timeseries_files)

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
                    "plot.profile_heatmaps.R {} {} {} {} {} "
                    "1,100 101,200 201,300").format(
                        out_mat_file,
                        histone, 
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

        if not os.path.isdir("{}/homer_HOCOMOCO".format(bioinformatics_bed_dir)):
            
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
                    bioinformatics_bed_dir,
                    mknown=args.outputs["annotations"]["pwms.renamed.nonredundant"],
                    mknown_name="HOCOMOCO")
                
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

    # also make a BED file of accessible stable but NON histone marked regions
    # this is approximate (using the mark groups, which are filtered for size)
    stable_no_histones_bed = "{}/clusters/{}.no_histones.bed.gz".format(
        results_dir, prefix)
    if not os.path.isfile(stable_no_histones_bed):
        get_stable_no_histones = (
            "bedtools subtract "
            "-a {0} "
            "-b {1}/*cluster*bed.gz | "
            "sort -k1,1 -k2,2n | "
            "bedtools merge -i stdin | "
            "gzip -c > {2}").format(
                args.outputs["data"]["atac.master.bed"],
                mark_bed_dir,
                stable_no_histones_bed)
        run_shell_cmd(get_stable_no_histones)
    
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
            if True:
                # adjust fake clusters file here to get proper row seps
                subsample_mat_key = "{}.mat".format(subsample_bed_key.split(".bed")[0])
                plot_atac_ordered_subsample_file = "{}.plot.mat".format(
                    out_results[subsample_mat_key].split(".mat")[0])
                fake_clusters_df = pd.read_table(
                    out_results[subsample_mat_key])
                colnames = list(fake_clusters_df.columns)
                colnames[-1] = "cluster"
                colnames[0] = "fake_cluster"
                fake_clusters_df.columns = colnames
                fake_clusters_df.to_csv(
                    plot_atac_ordered_subsample_file,
                    sep="\t", index=False)
                plot_clusters(
                    out_results[atac_fake_clusters_key],
                    plot_atac_ordered_subsample_file,
                    out_data[atac_stable_mat_key],
                    plot_dir,
                    os.path.basename(plot_atac_ordered_subsample_file), # prefix
                    plot_individual=False)

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
                row_sep_file = "{}/{}.row_seps.txt".format(plot_dir, os.path.basename(plot_atac_ordered_subsample_file))
                out_mat_file = "{}.point.mat.gz".format(out_prefix)
                out_r_file = "{}.replot.pdf".format(out_file.split(".pdf")[0])
                if not os.path.isfile(out_r_file):
                    replot = (
                        "plot.profile_heatmaps.R {} {} {} {} {} "
                        "1,100 101,200 201,300").format(
                            out_mat_file,
                            histone,
                            row_sep_file,
                            out_r_file,
                            histone_r_color)
                    run_shell_cmd(replot)

    # -----------------------------------------
    # ANALYSIS 5 - bioinformatics
    # inputs: BED dirs
    # outputs: HOMER/GREAT results
    # -----------------------------------------
    logger.info("ANALYSIS: run HOMER and GREAT")
    bioinformatics_bed_dirs = [mark_bed_dir, state_bed_dir]
    for bioinformatics_bed_dir in bioinformatics_bed_dirs:

        if not os.path.isdir("{}/homer_HOCOMOCO".format(bioinformatics_bed_dir)):
            
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
                    bioinformatics_bed_dir,
                    mknown=args.outputs["annotations"]["pwms.renamed.nonredundant.homer"],
                    mknown_name="HOCOMOCO")


    return args


def run_chromatin_states_workflow(args, prefix):
    """aggregate information into chromatin states map
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: aggregate chromatin states")

    # assertions
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "summary"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    # out file
    state_summary_key = "chrom_states.mat"
    state_summary_file = "{}/chrom_states.summary.txt".format(results_dir)
    out_results[state_summary_key] = state_summary_file
    if os.path.isfile(out_results[state_summary_key]):
        return args
    
    # load in relevant matrices with data
    # ATAC mat, H3K27ac mat, H3K4me1 mat, H3K27me3 mat, overlaps
    _MIN_REGION_NUM = 500
    atac_rlog_mat = pd.read_csv(out_data["atac.counts.pooled.rlog.mat"], sep="\t")
    H3K27ac_rlog_mat = pd.read_csv(out_data["H3K27ac.counts.pooled.rlog.mat"], sep="\t")
    H3K4me1_rlog_mat = pd.read_csv(out_data["H3K4me1.counts.pooled.rlog.mat"], sep="\t")
    H3K27me3_rlog_mat = pd.read_csv(out_data["H3K27me3.counts.pooled.rlog.mat"], sep="\t")

    # -----------------------------------------
    # ANALYSIS  - pull in dynamic set
    # inputs: BED dirs
    # outputs: data mat of chrom states
    # -----------------------------------------
    atac_dpgp_dir = args.outputs["results"]["atac"]["timeseries"]["dp_gp"]["dir"]
    traj_dir = "{}/reproducible/hard/reordered/".format(atac_dpgp_dir)
    traj_region_id_files = sorted(glob.glob("{}/*cluster*txt.gz".format(traj_dir)))
    
    # load in overlaps
    overlaps_dir = "{}/overlap_histone".format(
        args.outputs["results"]["epigenome"]["dynamic"]["dir"])
    H3K27ac_overlaps = convert_overlaps_bed_to_id_mappings(
        "{}/ggr.atac.ends.counts.pooled.rlog.dynamic.overlap.H3K27ac.tmp.bed.gz".format(
            overlaps_dir),
        extend_len=args.inputs["params"]["histones"]["H3K27ac"]["overlap_extend_len"])
    H3K4me1_overlaps = convert_overlaps_bed_to_id_mappings(
        "{}/ggr.atac.ends.counts.pooled.rlog.dynamic.overlap.H3K4me1.tmp.bed.gz".format(
            overlaps_dir),
        extend_len=args.inputs["params"]["histones"]["H3K4me1"]["overlap_extend_len"])
    overlaps_dir = "{}/overlap_histone.inactive".format(
        args.outputs["results"]["epigenome"]["dynamic"]["dir"])
    H3K27me3_overlaps = convert_overlaps_bed_to_id_mappings(
        "{}/ggr.atac.ends.counts.pooled.rlog.dynamic.overlap.H3K27me3.tmp.bed.gz".format(
            overlaps_dir),
        extend_len=args.inputs["params"]["histones"]["H3K27me3"]["overlap_extend_len"])
    dynamic_histones = [
        ("H3K27ac", H3K27ac_overlaps, H3K27ac_rlog_mat),
        ("H3K4me1", H3K4me1_overlaps, H3K4me1_rlog_mat),
        ("H3K27me3", H3K27me3_overlaps, H3K27me3_rlog_mat)]
    
    # set up chrom states/marks dirs
    epigenome_dynamic_dir = args.outputs["results"]["epigenome"]["dynamic"]["dir"]
    mark_region_id_dir = "{}/clusters/by_mark/ids".format(epigenome_dynamic_dir)
    state_region_id_dir = "{}/clusters/by_state/ids".format(epigenome_dynamic_dir)
    
    total_states = 0
    full_summary = None # track: TRAJ (15), ATAC (10), H3K27ac (3), H3K4me1 (3), H3K27me3 (3)
    trajectories = ["TRAJ.{}".format(val+1) for val in range(len(traj_region_id_files))]
    for traj_region_id_file in traj_region_id_files:

        # get cluster num and matching chrom mark/state files
        cluster_prefix = os.path.basename(traj_region_id_file).split(".")[-3]
        cluster_num = int(cluster_prefix.split("_")[-1])
        mark_files = glob.glob("{}/*atac-{}.*.txt.gz".format(mark_region_id_dir, cluster_prefix))
        state_files = glob.glob("{}/*atac-{}.*.txt.gz".format(state_region_id_dir, cluster_prefix))
        print cluster_num
        print cluster_prefix
        print mark_files,
        print state_files

        # track regions
        traj_regions = pd.read_csv(traj_region_id_file, sep="\t", header=None).iloc[:,0].values
        
        # check chrom states first. if chrom states exist, use those
        for state_file in state_files:
            print state_file

            # get the region ids
            state_name = os.path.basename(state_file).split(".")[-3]
            state_regions = pd.read_csv(state_file, sep="\t", header=None).iloc[:,0].values

            # extract data
            state_summary = get_aggregate_chromatin_state_summary(
                state_regions, trajectories, cluster_num, atac_rlog_mat, dynamic_histones,
                index=state_name)
            if full_summary is None:
                full_summary = state_summary.copy()
            else:
                full_summary = pd.concat([full_summary, state_summary], axis=0)

            # keep track of whatever is left over
            traj_regions = traj_regions[np.isin(traj_regions, state_regions, invert=True)]
            print traj_regions.shape
            
            total_states += 1
                
        # if many regions still left (>500?), go to marks and use those.
        print "after chrom states, num remaining:", traj_regions.shape[0]
        if traj_regions.shape[0] >= _MIN_REGION_NUM:
            used_mark_regions = [] # TODO add to here as going through, and then reduce out at end
            for mark_file in mark_files:

                # load, and only keep if > 500 in traj_regions
                mark_regions = pd.read_csv(mark_file, sep="\t", header=None).iloc[:,0].values
                mark_not_used_indices = np.where(np.isin(mark_regions, traj_regions))[0]
                if mark_not_used_indices.shape[0] < _MIN_REGION_NUM:
                    continue

                # if keeping, extract data
                print mark_file
                mark_name = os.path.basename(mark_file).split(".")[-3]
                mark_not_used_regions = mark_regions[mark_not_used_indices]
                mark_summary = get_aggregate_chromatin_state_summary(
                    mark_not_used_regions, trajectories, cluster_num, atac_rlog_mat, dynamic_histones,
                    index=mark_name)
                if full_summary is None:
                    full_summary = mark_summary.copy()
                else:
                    full_summary = pd.concat([full_summary, mark_summary], axis=0)

                #print mark_not_used_regions.shape
                used_mark_regions += list(mark_not_used_regions)
                print len(used_mark_regions)

            # keep track of whatever is left over
            traj_regions = traj_regions[np.isin(traj_regions, used_mark_regions, invert=True)]

        # if remainder in traj that are NOT in union (chrom states, marks)
        # is greater than 500, include as null state (just ATAC, no chrom marks)
        if traj_regions.shape[0] >= _MIN_REGION_NUM:
            group_name = "Null"
            null_summary = get_aggregate_chromatin_state_summary(
                traj_regions, trajectories, cluster_num, atac_rlog_mat, dynamic_histones,
                index=group_name)
            if full_summary is None:
                full_summary = null_summary.copy()
            else:
                full_summary = pd.concat([full_summary, null_summary], axis=0)
                
    # -----------------------------------------
    # ANALYSIS  - pull in stable set
    # inputs: BED dirs
    # outputs: data mat of chrom states
    # -----------------------------------------
    stable_regions = pd.read_csv(
        out_data["atac.counts.pooled.rlog.stable.mat"],
        sep="\t", index_col=0).index.values
    
    # load in overlaps
    overlaps_dir = "{}/overlap_histone".format(
        args.outputs["results"]["epigenome"]["stable"]["dir"])
    H3K27ac_overlaps = convert_overlaps_bed_to_id_mappings(
        "{}/ggr.atac.ends.counts.pooled.rlog.stable.overlap.H3K27ac.tmp.bed.gz".format(
            overlaps_dir),
        extend_len=args.inputs["params"]["histones"]["H3K27ac"]["overlap_extend_len"])
    H3K4me1_overlaps = convert_overlaps_bed_to_id_mappings(
        "{}/ggr.atac.ends.counts.pooled.rlog.stable.overlap.H3K4me1.tmp.bed.gz".format(
            overlaps_dir),
        extend_len=args.inputs["params"]["histones"]["H3K4me1"]["overlap_extend_len"])
    H3K27me3_overlaps = convert_overlaps_bed_to_id_mappings(
        "{}/ggr.atac.ends.counts.pooled.rlog.stable.overlap.H3K27me3.tmp.bed.gz".format(
            overlaps_dir),
        extend_len=args.inputs["params"]["histones"]["H3K27me3"]["overlap_extend_len"])
    dynamic_histones = [
        ("H3K27ac", H3K27ac_overlaps, H3K27ac_rlog_mat),
        ("H3K4me1", H3K4me1_overlaps, H3K4me1_rlog_mat),
        ("H3K27me3", H3K27me3_overlaps, H3K27me3_rlog_mat)]
    
    # set up chrom states/marks dirs
    epigenome_stable_dir = args.outputs["results"]["epigenome"]["stable"]["dir"]
    mark_region_id_dir = "{}/clusters/by_mark/ids".format(epigenome_stable_dir)
    state_region_id_dir = "{}/clusters/by_state/ids".format(epigenome_stable_dir)
    mark_files = sorted(glob.glob("{}/*.txt.gz".format(mark_region_id_dir)))
    state_files = sorted(glob.glob("{}/*.txt.gz".format(state_region_id_dir)))
    
    # look at chrom states
    for state_file in state_files:
        print state_file
        
        # get the region ids
        state_name = os.path.basename(state_file).split(".")[-3]
        state_regions = pd.read_csv(state_file, sep="\t", header=None).iloc[:,0].values

        # extract data
        state_summary = get_aggregate_chromatin_state_summary(
            state_regions, trajectories, None, atac_rlog_mat, dynamic_histones,
            index=state_name)
        if full_summary is None:
            full_summary = state_summary.copy()
        else:
            full_summary = pd.concat([full_summary, state_summary], axis=0)
            
        # keep track of whatever is left over
        stable_regions = stable_regions[np.isin(stable_regions, state_regions, invert=True)]
        print stable_regions.shape
            
        total_states += 1

    # sort and save out
    sort_columns = trajectories + ["H3K27ac.max", "H3K4me1.max", "H3K27me3.max"]
    full_summary = full_summary.sort_values(sort_columns, ascending=False)
    full_summary.to_csv(out_results[state_summary_key], sep="\t")

    return args


def _agg_deeptools_mat_data(mat_file, assay, extend_dist, bin_total):
    """summarize agg data
    """
    mat = pd.read_csv(mat_file, sep="\t", header=None, comment="@")
    mat = mat.iloc[:,6:]
    mat = mat.mean(axis=0)
    positions = np.arange(-extend_dist, extend_dist, bin_total)
    positions = np.tile(positions, 3)
    days = np.concatenate([
        np.tile(["d0"], bin_total),
        np.tile(["d3"], bin_total),
        np.tile(["d6"], bin_total)])
    agg_mat = pd.DataFrame({
        "value": mat,
        "position": positions,
        "day": days})
    agg_mat["assay"] = assay
    
    return agg_mat


def plot_profile_heatmaps_workflow(args, tss_file, plot_prefix):
    """plot histone heatmaps for a specific set of tss regions (ordered
    """
    # make a subsample file if necessary
    atac_plot_file = "{}.ATAC.heatmap.profile.pdf".format(plot_prefix)
    if not os.path.isfile(atac_plot_file):
        tss_data = pd.read_csv(tss_file, sep="\t")
        if tss_data.shape[0] > 1000:
            step_size = int(tss_data.shape[0] / 1000)
            keep_indices = np.linspace(0, tss_data.shape[0]-1, num=1000).astype(int)
            #keep_indices = np.arange(0, tss_data.shape[0], step=step_size)
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
    #if not os.path.isfile("{}.atac.heatmap.pdf".format(plot_prefix)):
    if False:
        r_cmd = "~/git/ggr-project/R/plot.tss.timeseries_heatmap.R {} {} {}.atac ATAC-seq".format(
            "{}.promoter_data.mat.txt.gz".format(plot_prefix),
            "{}.row_seps.txt.gz".format(plot_prefix),
            plot_prefix)
        print r_cmd
        os.system(r_cmd)

    # plot ATAC seq as profile map
    atac_bigwigs = sorted(
        glob.glob("{}/{}".format(
            args.inputs["atac"][args.cluster]["data_dir"],
            args.inputs["atac"][args.cluster]["bigwig_pooled_glob"])))
    atac_bigwigs = [atac_bigwigs[0], atac_bigwigs[6], atac_bigwigs[12]]

    out_prefix = "{}.ATAC".format(plot_prefix)
    plot_file = "{}.heatmap.profile.pdf".format(out_prefix)
    if not os.path.isfile(plot_file):
        make_deeptools_heatmap(
            tss_file,
            atac_bigwigs,
            out_prefix,
            sort=False,
            referencepoint="center",
            extend_dist=5000, #10000
            bin_total=100, # note: do not adjust
            color="Blues")
        
    # make profile heatmap in R, with strand oriented TSS
    row_sep_file = "{}.row_seps.txt.gz".format(plot_prefix)
    out_mat_file = "{}.point.mat.gz".format(out_prefix)
    out_r_file = "{}.replot.pdf".format(plot_file.split(".pdf")[0])
    if not os.path.isfile(out_r_file):
        # TODO - note - strands already adjusted
        replot = (
            "~/git/ggr-project/R/plot.profile_heatmaps.not_stranded.R {} {} {} {} {} "
            "1,100 101,200 201,300").format(
                out_mat_file,
                "ATAC", 
                row_sep_file,
                out_r_file,
                "Blues")
        print replot
        run_shell_cmd(replot)
        
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
                args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_fc_bigwig_glob"])))

        out_prefix = "{}.{}_overlap".format(plot_prefix, histone)
        plot_file = "{}.heatmap.profile.pdf".format(out_prefix)
        if not os.path.isfile(plot_file):
            make_deeptools_heatmap(
                tss_file,
                histone_bigwigs,
                out_prefix,
                sort=False,
                referencepoint="center",
                extend_dist=5000, #10000, #15000,
                bin_total=100, # if change this, need to change R plot params below
                color=histone_color)

        # make profile heatmap in R, with strand oriented TSS
        row_sep_file = "{}.row_seps.txt.gz".format(plot_prefix)
        out_mat_file = "{}.point.mat.gz".format(out_prefix)
        out_r_file = "{}.replot.pdf".format(plot_file.split(".pdf")[0])
        if not os.path.isfile(out_r_file):
            # TODO - note - strands already adjusted
            replot = (
                "~/git/ggr-project/R/plot.profile_heatmaps.not_stranded.R {} {} {} {} {} "
                "1,100 101,200 201,300").format(
                    out_mat_file,
                    histone, 
                    row_sep_file,
                    out_r_file,
                    histone_r_color)
            print replot
            run_shell_cmd(replot)

    # and plot agg plot with all together in one plot
    extend_dist = 5000
    bin_total = 100
    
    # pull in ATAC, histone data
    agg_plot_prefix = "{}.agg".format(plot_prefix)
    plot_file = "{}.ALL.pdf".format(agg_plot_prefix)
    if not os.path.isfile(plot_file):
        atac_mat_file = "{}.ATAC.point.mat.gz".format(plot_prefix)
        mat_data = _agg_deeptools_mat_data(atac_mat_file, "ATAC", extend_dist, bin_total)
        for histone in histones:
            histone_mat_file = "{}.{}_overlap.point.mat.gz".format(plot_prefix, histone)
            histone_mat_data = _agg_deeptools_mat_data(histone_mat_file, histone, extend_dist, bin_total)
            mat_data = pd.concat([mat_data, histone_mat_data], axis=0)
        summary_mat_file = "{}.agg_data.melted.txt.gz".format(plot_prefix)
        mat_data.to_csv(
            summary_mat_file,
            header=True, index=False, compression="gzip", sep="\t")

        plot_cmd = "~/git/ggr-project/R/plot.profile_agg_plots.R {} {}".format(
            summary_mat_file, agg_plot_prefix)
        print plot_cmd
        os.system(plot_cmd)
    
    return args


def _run_spreading_workflow_unidirectional(
        args, results_dir, prefix, histone, spread_regions_file, peak_file,
        spread_bp=1000, direction="increase"):
    """keep code for each direction here, to make it easy to keep standardized
    """
    # params
    #spread_thresh = 5*spread_bp # TODO there's gotta be a better heuristic? maybe go by histone lengths?
    spread_thresh = 10*spread_bp # TODO there's gotta be a better heuristic? maybe go by histone lengths?
    
    # generate a day spread file (make day specific domain so that
    # when intersecting, won't capture smaller peaks multiple times)
    day_spread_regions_file = "{}/{}.{}.{}.day_spread_merge_1kb.bed.gz".format(
        results_dir, prefix, histone, direction)
    if not os.path.isfile(day_spread_regions_file):
        merge_cmd = (
            "zcat {0} | sort -k1,1 -k2,2n | "
            "bedtools slop -i stdin -g {1} -b {2} | "
            "bedtools merge -i stdin | "
            "bedtools slop -i stdin -g {1} -b -{2} | "
            "gzip -c > {3}").format(
                peak_file,
                args.inputs["annot"][args.cluster]["chromsizes"],
                spread_bp,
                day_spread_regions_file)
        print merge_cmd
        os.system(merge_cmd)

    # intersect for SPREAD: the domain increased (from the initiating region)
    intersected_spread_file = "{}.intersect_spread.bed.gz".format(
        day_spread_regions_file.split(".bed")[0])
    if not os.path.isfile(intersected_spread_file):
        get_dynamic_spread_cmd = (
            "bedtools intersect -wao -a {} -b {} | " # intersect while keeping all data
            "gzip -c > {}").format(
                spread_regions_file,
                peak_file, #day_spread_regions_file,
                intersected_spread_file)
        print get_dynamic_spread_cmd
        os.system(get_dynamic_spread_cmd)

    # read in and filter for significant spread
    # TO consider: build a distribution and then take outliers? or some actual stat test?
    filt_file = "{}_filt.bed.gz".format(
        intersected_spread_file.split(".bed")[0])
    if not os.path.isfile(filt_file):
        data = pd.read_csv(intersected_spread_file, sep="\t", header=None)
        data = data[data[3] != "."]
        data["domain"] = "domain=" + data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)
        data["init_region"] = "init_region=" + data[3] + ":" + data[
            4].map(str) + "-" + data[5].map(str) + ";score=" + data[9].map(str)
        data["start_diff"] = data[4] - data[1]
        data["stop_diff"] = data[2] - data[5]
        data["total_diff"] = data["start_diff"] + data["stop_diff"]
        data = data[data["total_diff"] > spread_thresh]
        print "after spread filters, {} domains".format(len(set(data["domain"].values.tolist())))
        data.to_csv(
            filt_file,
            columns=[0, 1, 2, "init_region"],
            sep="\t", compression="gzip", header=None, index=None)

    # intersect for DYNAMICS: the receiving region
    intersected_dynamics_file = "{}.intersect_dynamics.bed.gz".format(
        filt_file.split(".bed")[0])
    if not os.path.isfile(intersected_dynamics_file):
        get_dynamics_cmd = (
            "bedtools intersect -wao -a {} -b {} | " # intersect while keeping all data
            "gzip -c > {}").format(
                filt_file,
                args.outputs["results"]["histones"][histone]["timeseries"][
                    "ggr.histone.{}.enumerated.bed".format(histone)],
                intersected_dynamics_file)
        print get_dynamics_cmd
        os.system(get_dynamics_cmd)

    # filter for domains with dynamic recipient peaks only
    filt_file = "{}_filt.bed.gz".format(intersected_dynamics_file.split(".bed")[0])
    if not os.path.isfile(filt_file):
        data = pd.read_csv(intersected_dynamics_file, sep="\t", header=None)
        try:
            data = data[data[7] != "."]
            data[7] = data[7].astype(int)
        except TypeError:
            pass
        data["domain"] = "domain=" + data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)
        data["recip_region"] = "recip_region=" + data[4] + ":" + data[5].map(str) + "-" + data[6].map(str)
        # filter for correct DIRECTION of spread
        if direction == "increase" :
            data = data[data[7] < 4]
        elif direction == "decrease":
            data = data[data[7] > 4]
        else:
            raise ValueError, "direction not recognized!!"
        data = data.sort_values([7, "domain", "recip_region", 3])
        data["metadata"] = data["domain"] + ";" + data[
            3] + ";" + data["recip_region"] + ";cluster=" + data[7].map(str)
        data.to_csv(
            filt_file,
            columns=[4, 5, 6, "metadata"],
            sep="\t", compression="gzip", header=False, index=False)
        print "after filter for dynamic recipient regions, {} domains".format(
            len(set(data["domain"].values.tolist())))

    # overlap with TSS
    tss_file = args.outputs["data"]["tss.dynamic"]
    tss_overlap_file = "{}.tss_filt.bed.gz".format(filt_file.split(".bed")[0])
    if not os.path.isfile(tss_overlap_file):
        bedtools_cmd = (
            "bedtools slop -s -i {0} -g {1} -b {2} | " # increase histone range
            "bedtools intersect -wo -a {3} -b stdin | " # intersect
            "gzip -c > {4}").format(
                filt_file,
                args.inputs["annot"][args.cluster]["chromsizes"],
                spread_bp,
                args.outputs["data"]["tss.dynamic"],
                tss_overlap_file)
        print bedtools_cmd
        os.system(bedtools_cmd)
        
    # generate other output files
    plot_prefix = "{}/{}.{}.{}".format(
        results_dir, prefix, histone, direction)
    out_file = "{}.promoter_rna_data.mat.txt.gz".format(plot_prefix)
    tss_file = "{}.sorted.bed.gz".format(tss_overlap_file.split(".bed")[0])
    if not os.path.isfile(out_file):
        data = pd.read_csv(tss_overlap_file, sep="\t", header=None)

        # filter for gene trajectories
        rna_clusters = pd.read_csv(args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"], sep="\t")
        data = data.merge(rna_clusters, left_on=3, right_on="id")
        data = data.drop("id", axis=1)
        if True:
            if direction == "increase":
                data = data[data["cluster"] > 7].sort_values("cluster")
            elif direction == "decrease":
                data = data[data["cluster"] < 4].sort_values("cluster")
            else:
                raise ValueError, "unknown direction!"

        # add in hgnc ids
        mappings = pd.read_csv(args.outputs["annotations"]["geneids.mappings.mat"], sep="\t")
        data = data.merge(mappings, left_on=3, right_on="ensembl_gene_id")
        data = data.drop(["ensembl_gene_id", "entrezgene"], axis=1)

        # get an ordered TSS file (must simplify down to one gene per line for this)
        metadata = data[9].str.split(";", n=4, expand=True)
        metadata["score"] = metadata[2].str.split("=", n=2, expand=True)[1].astype(float)
        metadata["gene_id"] = data[3]
        # use score to reduce
        #metadata = metadata.loc[metadata.groupby([0])["score"].idxmax()] # this is an example/debug
        data = data.loc[metadata.groupby(["gene_id"])["score"].idxmax()]
        print "after grouping, {} tss".format(data.shape[0])
        # now sort by the init region distance from TSS
        metadata = metadata[1].str.split("=", n=2, expand=True)
        metadata[["chrom", "start-stop"]] = metadata[1].str.split(":", n=2, expand=True)
        metadata[["start", "stop"]] = metadata["start-stop"].str.split("-", n=2, expand=True)
        data["init_midpoint"] = metadata["start"].astype(int) + (metadata["stop"].astype(int) - metadata["start"].astype(int)) / 2
        data["init_to_tss"] = data[1] - data["init_midpoint"] # positive means upstream
        # reorient based on strand
        for row_i in range(data.shape[0]):
            index_val = data.index[row_i]
            if data[5].iloc[row_i] == "-":
                data.at[index_val, "init_to_tss"] = -1 * data["init_to_tss"][index_val]
        data = data.sort_values("init_to_tss", ascending=False)
        # also ignore those that are too close?
        #data = data[data["init_to_tss"].abs() > 1000] # TODO think about this
        print "after removing nearby, {} tss".format(data.shape[0])
        data.to_csv(
            tss_file, columns=[0,1,2,3,4,5],
            sep="\t", compression="gzip", header=False, index=False)
        
        # overlap PAS-seq
        pas_data = pd.read_csv(
            args.outputs["data"]["rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat"],
            sep="\t")
        data = data.merge(pas_data, left_on=3, right_index=True)
        out_file = "{}.promoter_rna_data.mat.txt.gz".format(plot_prefix)
        data.to_csv(out_file, sep="\t", compression="gzip", index=False)
    
    # TODO make a rowseps file?
    plot_profile_heatmaps_workflow(args, tss_file, plot_prefix)

    # TODO gene set enrichment
    
    
    return args



def run_histone_spreading_workflow(args, prefix):
    """analyze histone spreading
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: look at histone spreading")

    # assertions
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "histone_spread"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    histone = "H3K27ac"

    # do a merge on extended regions, to capture nearby histone peaks into larger domain
    master_bed_file = args.outputs["data"]["{}.master.bed".format(histone)]
    spread_regions_file = "{}.spread_1kb_merge.bed.gz".format(
        master_bed_file.split(".bed")[0])
    spread_bp = 1000 # this is approx 6 histones if wrapped?
    if not os.path.isfile(spread_regions_file):
        merge_cmd = (
            "bedtools slop -i {0} -g {1} -b {2} | "
            "bedtools merge -i stdin | "
            "bedtools slop -i stdin -g {1} -b -{2} | "
            "gzip -c > {3}").format(
                master_bed_file,
                args.inputs["annot"][args.cluster]["chromsizes"],
                spread_bp,
                spread_regions_file)
        print merge_cmd
        os.system(merge_cmd)

    # filter: only keep the domains that have multiple peaks
    intersected_spread_file = "{}.multi_peak.bed.gz".format(
        spread_regions_file.split(".bed")[0])
    if not os.path.isfile(intersected_spread_file):
        get_dynamic_spread_cmd = (
            "bedtools intersect -wao -a {} -b {} | " # intersect while keeping all data
            "gzip -c > {}").format(
                spread_regions_file,
                master_bed_file,
                intersected_spread_file)
        print get_dynamic_spread_cmd
        os.system(get_dynamic_spread_cmd)

    # read in to filter for the multiple peaks
    peak_num_thresh = 3
    domain_file = "{}.multi_peak_only.bed.gz".format(
        spread_regions_file.split(".bed")[0])
    if not os.path.isfile(domain_file):
        data = pd.read_csv(intersected_spread_file, sep="\t", header=None)
        data["domain"] = "domain=" + data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)
        data["indiv_region"] = "indiv_region=" + data[3] + ":" + data[4].map(str) + "-" + data[5].map(str)
        print "after merging with ext {}, found {} domains".format(
            spread_bp, len(set(data["domain"].values.tolist())))
        data_peak_check = data.copy()
        data_peak_check["peak_found"] = 1
        data_peak_check = data_peak_check[["domain", "peak_found"]].groupby("domain").sum()
        keep_ids = data_peak_check.index[data_peak_check["peak_found"] >= peak_num_thresh].values
        data = data[data["domain"].isin(keep_ids)]
        print "after filtering for multi-peak, found {} domains".format(
            len(set(data["domain"].values.tolist())))
        data = data[[0,1,2]].drop_duplicates()
        data.to_csv(
            domain_file,
            columns=[0,1,2],
            sep="\t", compression="gzip", header=False, index=False)

    # timepoint peak files
    peak_files = sorted(
        glob.glob(
            "{}/{}".format(
                args.inputs["chipseq"][args.cluster]["data_dir"],
                args.inputs["chipseq"][args.cluster]["histones"][histone]["overlap_glob"])))

    # now look for domains that INCREASED in spread over time
    # extend regions in d0 to match domains file - this file is the baseline for SPREAD
    args = _run_spreading_workflow_unidirectional(
        args, results_dir, prefix, histone, domain_file, peak_files[0],
        spread_bp=spread_bp, direction="increase")
    args = _run_spreading_workflow_unidirectional(
        args, results_dir, prefix, histone, domain_file, peak_files[2],
        spread_bp=spread_bp, direction="decrease")

    # TODO maybe also looped regions with HiChIP support?

    return args


def run_distal_v_proximal_workflow(args, prefix):
    """analyze distal v proximal 
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: look at distal v proximal regions")

    # assertions
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "distal_v_proximal"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    # inputs
    atac_master_bed = args.outputs["data"]["atac.master.bed"]
    atac_dynamic_bed = args.outputs["data"]["atac.dynamic.bed"] # note that this file is NOT traj filt
    atac_stable_bed = args.outputs["data"]["atac.stable.bed"]
    tss_all_bed = args.outputs["annotations"]["tss.pc.bed"]
    tss_dynamic_bed = args.outputs["data"]["tss.dynamic"]
    tss_stable_bed = args.outputs["data"]["tss.stable"]
    tss_nonexpr_bed = args.outputs["data"]["tss.non_expr"]

    # adjust the atac dynamic for a fold change thresh
    # NOTE: currently not used
    atac_dynamic_filt_bed = "{}/{}.atac_fc_filt.dynamic.bed.gz".format(
        results_dir, prefix)
    if not os.path.isfile(atac_dynamic_filt_bed):
        atac_traj_mat = pd.read_csv(
            args.outputs["data"]["atac.counts.pooled.rlog.dynamic.traj.mat"],
            sep="\t", index_col=0)
        max_vals = atac_traj_mat.max(axis=1)
        min_vals = atac_traj_mat.min(axis=1)
        atac_traj_mat["max"] = max_vals
        atac_traj_mat["min"] = min_vals
        # keep only those with minimum 2-fold change
        atac_traj_mat = atac_traj_mat[atac_traj_mat["max"] - atac_traj_mat["min"] > 1]
        atac_traj_mat = atac_traj_mat.reset_index()
        atac_traj_mat[["chrom", "start-stop"]] = atac_traj_mat["index"].str.split(":", n=2, expand=True)
        atac_traj_mat[["start", "stop"]] = atac_traj_mat["start-stop"].str.split("-", n=2, expand=True)
        atac_traj_mat.to_csv(
            atac_dynamic_filt_bed,
            columns=["chrom", "start", "stop"],
            sep="\t", compression="gzip", header=False, index=False)

    atac_files = [atac_master_bed, atac_dynamic_bed, atac_stable_bed]
    tss_files = [tss_dynamic_bed, tss_stable_bed, tss_nonexpr_bed]

    for atac_file in atac_files:
        for tss_file in tss_files:
            # overlap file
            overlap_file = "{}/{}.atac_{}.rna_{}.bed.gz".format(
                results_dir, prefix,
                os.path.basename(atac_file).split(".")[-3],
                os.path.basename(tss_file).split(".")[-3])

            if not os.path.isfile(overlap_file):
                # take ATAC bed and overlap with (extended) TSS bed
                bedtools_cmd = (
                    "bedtools slop -s -i {0} -g {1} -b {2} | " # increase TSS range
                    "bedtools intersect -wo -a {3} -b stdin | " # intersect
                    "gzip -c > {4}").format(
                        tss_file,
                        args.inputs["annot"][args.cluster]["chromsizes"],
                        args.inputs["params"]["tss.extend_len"],
                        atac_file,
                        overlap_file)
                print bedtools_cmd
                os.system(bedtools_cmd)

    # ANALYSIS - look at range of fold changes for distal vs proximal ATAC
    foldchanges_by_group_file = "{}/{}.atac_foldchanges.by_tss_group.txt.gz".format(
        results_dir, prefix)
    if not os.path.isfile(foldchanges_by_group_file):

        # first pull in atac data
        atac_mat = pd.read_csv(
            args.outputs["data"]["atac.counts.pooled.rlog.mat"], sep="\t")
        max_val = atac_mat.max(axis=1)
        min_val = atac_mat.min(axis=1)
        atac_mat["max_val"] = max_val
        atac_mat["min_val"] = min_val
        atac_mat = atac_mat[["max_val", "min_val"]]
        atac_mat["foldchange"] = atac_mat["max_val"] - atac_mat["min_val"]

        # for each tss group, get ATAC ids
        for tss_file_idx in range(len(tss_files)):
            tss_file = tss_files[tss_file_idx]
            tss_group = os.path.basename(tss_file).split(".")[-3]

            atac_ids = []
            for atac_file in atac_files:
                # overlap file
                overlap_file = "{}/{}.atac_{}.rna_{}.bed.gz".format(
                    results_dir, prefix,
                    os.path.basename(atac_file).split(".")[-3],
                    tss_group)
                match_data = pd.read_csv(
                    overlap_file, sep="\t", header=None)
                match_data["atac_id"] = match_data[0] + ":" + match_data[1].map(str) + "-" + match_data[2].map(str)
                atac_ids += match_data["atac_id"].values.tolist()
                
            # now select those in atac ids and save out
            tss_group_atac_mat = atac_mat[atac_mat.index.isin(atac_ids)].copy()
            tss_group_atac_mat["group"] = "{} genes (n={})".format(tss_group, len(set(atac_ids)))
            if tss_file_idx == 0:
                tss_atac_mat = tss_group_atac_mat.copy()
            else:
                tss_atac_mat = pd.concat([tss_atac_mat, tss_group_atac_mat], axis=0)
            print tss_atac_mat.shape

        # now all the rest are distal
        distal_atac_mat = atac_mat[~atac_mat.index.isin(tss_atac_mat.index)].copy()
        distal_atac_mat["group"] = "Distal regions (n={})".format(distal_atac_mat.shape[0])
        tss_atac_mat = pd.concat([tss_atac_mat, distal_atac_mat], axis=0)
        tss_atac_mat.to_csv(
            foldchanges_by_group_file,
            columns=["foldchange", "group"],
            sep="\t", header=True, index=False, compression="gzip")
        print tss_atac_mat.shape
        
    # and plot
    plot_file = "{}/{}.atac_fc.by_tss_group.pdf".format(results_dir, prefix)
    if not os.path.isfile(plot_file):
        plot_cmd = "~/git/ggr-project/R/plot.atac_fc.by_tss_group.R {} {}".format(
            foldchanges_by_group_file, plot_file)
        print plot_cmd
        os.system(plot_cmd)
                
    return args


def _run_epigenome_summary_workflow(args, prefix):
    """build a top view summary of the skin epigenome
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: look at distal v proximal regions")

    # assertions
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "summary"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    # FINAL output files
    summary_mat_file = "{}/{}.summary_mat.txt.gz".format(results_dir, prefix)
    plot_mat_file = "{}/{}.summary_mat.plot.txt.gz".format(results_dir, prefix)
    out_results["summary_mat"] = summary_mat_file
    
    # make a union file of non ATAC centric regions to be able to add to plot
    work_dir = "{}/no_atac".format(results_dir)
    if not os.path.isdir(work_dir):
        os.system("mkdir -p {}".format(work_dir))

    no_atac_master_bed = "{}/ggr.no_atac.master.bed.gz".format(work_dir)
    if not os.path.isfile(no_atac_master_bed):
        no_atac_files = []
        
        # get no atac TSSs
        master_no_atac_file = "{}/ggr.tss.no_atac.bed.gz".format(work_dir)
        if not os.path.isfile(master_no_atac_file):
            master_tss_file = args.outputs["annotations"]["tss.pc.bed"]
            data = pd.read_csv(master_tss_file, sep="\t", header=None)
            data["id"] = data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)

            # use atac master to remove atac-proximal marks
            tmp_intersect_file = "{}/ggr.tss.atac_overlap.match.bed.gz".format(work_dir)
            bedtools_cmd = (
                "bedtools slop -i {0} -g {1} -b {2} | "
                "bedtools intersect -u -a {3} -b stdin | "
                "gzip -c > {4}").format(
                    args.outputs["data"]["atac.master.bed"],
                    args.inputs["annot"][args.cluster]["chromsizes"],
                    args.inputs["params"]["tss.extend_len"],
                    master_tss_file,
                    tmp_intersect_file)
            print bedtools_cmd
            os.system(bedtools_cmd)
            remove_data = pd.read_csv(tmp_intersect_file, sep="\t", header=None)
            remove_data["id"] = remove_data[0] + ":" + remove_data[1].map(str) + "-" + remove_data[2].map(str)
            data = data[~data["id"].isin(remove_data["id"])]

            # save out
            data.to_csv(
                master_no_atac_file, columns=[0,1,2],
                sep="\t", compression="gzip", header=False, index=False)
        no_atac_files.append(master_no_atac_file)
            
        # histones
        histones = ["H3K27ac", "H3K4me1", "H3K27me3"]
        for histone in histones:
            master_no_atac_file = "{}/ggr.{}.no_atac.bed.gz".format(work_dir, histone)
            if not os.path.isfile(master_no_atac_file):
            
                # load master ids
                master_bed_file = args.outputs["data"]["{}.master.bed".format(histone)]
                data = pd.read_csv(master_bed_file, sep="\t", header=None)
                data["id"] = data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)

                # use atac master to remove atac-proximal marks
                tmp_intersect_file = "{}/ggr.{}.atac_overlap.match.bed.gz".format(work_dir, histone)
                bedtools_cmd = (
                    "bedtools slop -i {0} -g {1} -b {2} | "
                    "bedtools intersect -u -a {3} -b stdin | "
                    "gzip -c > {4}").format(
                        args.outputs["data"]["atac.master.bed"],
                        args.inputs["annot"][args.cluster]["chromsizes"],
                        1000,
                        master_bed_file,
                        tmp_intersect_file)
                print bedtools_cmd
                os.system(bedtools_cmd)
                remove_data = pd.read_csv(tmp_intersect_file, sep="\t", header=None)
                remove_data["id"] = remove_data[0] + ":" + remove_data[1].map(str) + "-" + remove_data[2].map(str)
                data = data[~data["id"].isin(remove_data["id"])]

                # save out
                data.to_csv(
                    master_no_atac_file, columns=[0,1,2],
                    sep="\t", compression="gzip", header=False, index=False)
            no_atac_files.append(master_no_atac_file)

        # then union across all of these
        bedtools_cmd = "zcat {} | sort -k1,1 -k2,2n | bedtools merge -i stdin | gzip -c > {}".format(
            " ".join(no_atac_files), no_atac_master_bed)
        print bedtools_cmd
        os.system(bedtools_cmd)
        
    # now go through and pull in matching data
    # tss - need to get dynamic, stable, off
    tss_dynamic_bed = args.outputs["data"]["tss.dynamic"]
    tss_stable_bed = args.outputs["data"]["tss.stable"]
    tss_nonexpr_bed = args.outputs["data"]["tss.non_expr"]
    tss_files = [tss_dynamic_bed, tss_stable_bed, tss_nonexpr_bed]

    for tss_file in tss_files:
        overlap_file = "{}/{}.no_atac.rna_{}.bed.gz".format(
            work_dir, prefix, os.path.basename(tss_file).split(".")[-3])
        if not os.path.isfile(overlap_file):
            bedtools_cmd = (
                "bedtools slop -s -i {0} -g {1} -b {2} | " # increase TSS range
                "bedtools intersect -wo -a {3} -b stdin | " # intersect
                "gzip -c > {4}").format(
                    tss_file,
                    args.inputs["annot"][args.cluster]["chromsizes"],
                    args.inputs["params"]["tss.extend_len"],
                    no_atac_master_bed,
                    overlap_file)
            print bedtools_cmd
            os.system(bedtools_cmd)

    # histone marks - use get_histone_overlaps
    histones = ["H3K27ac", "H3K4me1", "H3K27me3"]
    histone_files = []
    for histone in histones:
        # need to make files look like enumerated files
        master_bed_file = args.outputs["data"]["{}.master.bed".format(histone)]
        histone_enum_file = "{}/{}.{}.no_atac.enum.bed.gz".format(work_dir, prefix, histone)
        if not os.path.isfile(histone_enum_file):
            data = pd.read_csv(master_bed_file, sep="\t", header=None)
            data["enum"] = 1
            data.to_csv(histone_enum_file, sep="\t", compression="gzip", header=False, index=False)

        histone_files.append((
            histone, histone_enum_file,
            args.inputs["params"]["histones"][histone]["overlap_extend_len"]))

    histone_overlaps_file = "{}/ggr.histones.no_atac.histone_overlap.mat.txt.gz".format(work_dir, histone)
    if not os.path.isfile(histone_overlaps_file):
        get_histone_overlaps(
            no_atac_master_bed,
            histone_files,
            work_dir,
            histone_overlaps_file,
            args.inputs["annot"][args.cluster]["chromsizes"])

    # now pull it all together for the part of the matrix that is non-ATAC centric
    if not os.path.isfile(out_results["summary_mat"]):

        data = pd.read_csv(no_atac_master_bed, sep="\t", header=None)
        data["atac_id"] = data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)
        data = data[["atac_id"]] # not really atac id, but use dummy id to match the atac-centric part
        data["atac_cluster"] = -1

        # pull in TSS data (dynamic, stable, off)
        tss_dynamic_overlap_file = "{}/ggr.epigenome.summary.no_atac.rna_dynamic.bed.gz".format(work_dir)
        tss_dynamic_overlap = pd.read_csv(tss_dynamic_overlap_file, sep="\t", header=None)
        tss_dynamic_overlap["atac_id_tss"] = tss_dynamic_overlap[
            0] + ":" + tss_dynamic_overlap[1].map(str) + "-" + tss_dynamic_overlap[2].map(str)
        tss_dynamic_overlap = tss_dynamic_overlap[["atac_id_tss", 6]]
        tss_dynamic_overlap.columns = ["atac_id_tss", "rna_id"]
        tss_dynamic_overlap["rna_pattern"] = "dynamic"

        tss_stable_overlap_file = "{}/ggr.epigenome.summary.no_atac.rna_stable.bed.gz".format(work_dir)
        tss_stable_overlap = pd.read_csv(tss_stable_overlap_file, sep="\t", header=None)
        tss_stable_overlap["atac_id_tss"] = tss_stable_overlap[
            0] + ":" + tss_stable_overlap[1].map(str) + "-" + tss_stable_overlap[2].map(str)
        tss_stable_overlap = tss_stable_overlap[["atac_id_tss", 6]]
        tss_stable_overlap.columns = ["atac_id_tss", "rna_id"]
        tss_stable_overlap["rna_pattern"] = "stable"

        tss_nonexpr_overlap_file = "{}/ggr.epigenome.summary.no_atac.rna_non_expr.bed.gz".format(work_dir)
        tss_nonexpr_overlap = pd.read_csv(tss_nonexpr_overlap_file, sep="\t", header=None)
        tss_nonexpr_overlap["atac_id_tss"] = tss_nonexpr_overlap[
            0] + ":" + tss_nonexpr_overlap[1].map(str) + "-" + tss_nonexpr_overlap[2].map(str)
        tss_nonexpr_overlap = tss_nonexpr_overlap[["atac_id_tss", 6]]
        tss_nonexpr_overlap.columns = ["atac_id_tss", "rna_id"]
        tss_nonexpr_overlap["rna_pattern"] = "nonexpr"

        tss_data = pd.concat(
            [tss_dynamic_overlap, tss_stable_overlap, tss_nonexpr_overlap], axis=0)
        data = data.merge(tss_data, how="left", left_on="atac_id", right_on="atac_id_tss")
        data = data.drop(["atac_id_tss"], axis=1)

        # pull in TSS dynamic traj
        tss_traj = pd.read_csv(
            args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"], sep="\t")
        tss_traj.columns = ["rna_cluster", "id"]
        data = data.merge(tss_traj, how="left", left_on="rna_id", right_on="id")
        data["rna_cluster"] = data["rna_cluster"].fillna(0).astype(int)
        data = data.drop("id", axis=1)

        # pull in histone mark data
        histone_data = pd.read_csv(histone_overlaps_file, sep="\t")
        histone_data = histone_data[["id", "H3K27ac", "H3K4me1", "H3K27me3"]]
        data = data.merge(histone_data, how="left", left_on="atac_id", right_on="id")
        data = data.drop("id", axis=1)
        
        data_no_atac = data.copy()
    
    if not os.path.isfile(out_results["summary_mat"]):
        
        # pull in ATAC region ids
        data = pd.read_csv(args.outputs["data"]["atac.master.bed"], sep="\t", header=None)
        data["atac_id"] = data[0] + ":" + data[1].map(str) + "-" + data[2].map(str)
        data = data[["atac_id"]]

        # pull in ATAC dynamic (1-15) vs stable (0)
        atac_dynamic = pd.read_csv(args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"], sep="\t")
        data = data.merge(atac_dynamic, how="left", left_on="atac_id", right_on="id").fillna(0)
        data["atac_cluster"] = data["cluster"].astype(int)
        data = data.drop(["cluster", "id"], axis=1)

        # pull in TSS data (dynamic, stable, off)
        tss_dynamic_overlap_file = "{}/ggr.epigenome.distal_v_proximal.atac_master.rna_dynamic.bed.gz".format(
            args.outputs["results"]["epigenome"]["distal_v_proximal"]["dir"])
        tss_dynamic_overlap = pd.read_csv(tss_dynamic_overlap_file, sep="\t", header=None)
        tss_dynamic_overlap["atac_id_tss"] = tss_dynamic_overlap[
            0] + ":" + tss_dynamic_overlap[1].map(str) + "-" + tss_dynamic_overlap[2].map(str)
        tss_dynamic_overlap = tss_dynamic_overlap[["atac_id_tss", 6]]
        tss_dynamic_overlap.columns = ["atac_id_tss", "rna_id"]
        tss_dynamic_overlap["rna_pattern"] = "dynamic"

        tss_stable_overlap_file = "{}/ggr.epigenome.distal_v_proximal.atac_master.rna_stable.bed.gz".format(
            args.outputs["results"]["epigenome"]["distal_v_proximal"]["dir"])
        tss_stable_overlap = pd.read_csv(tss_stable_overlap_file, sep="\t", header=None)
        tss_stable_overlap["atac_id_tss"] = tss_stable_overlap[
            0] + ":" + tss_stable_overlap[1].map(str) + "-" + tss_stable_overlap[2].map(str)
        tss_stable_overlap = tss_stable_overlap[["atac_id_tss", 6]]
        tss_stable_overlap.columns = ["atac_id_tss", "rna_id"]
        tss_stable_overlap["rna_pattern"] = "stable"

        tss_nonexpr_overlap_file = "{}/ggr.epigenome.distal_v_proximal.atac_master.rna_non_expr.bed.gz".format(
            args.outputs["results"]["epigenome"]["distal_v_proximal"]["dir"])
        tss_nonexpr_overlap = pd.read_csv(tss_nonexpr_overlap_file, sep="\t", header=None)
        tss_nonexpr_overlap["atac_id_tss"] = tss_nonexpr_overlap[
            0] + ":" + tss_nonexpr_overlap[1].map(str) + "-" + tss_nonexpr_overlap[2].map(str)
        tss_nonexpr_overlap = tss_nonexpr_overlap[["atac_id_tss", 6]]
        tss_nonexpr_overlap.columns = ["atac_id_tss", "rna_id"]
        tss_nonexpr_overlap["rna_pattern"] = "nonexpr"

        tss_data = pd.concat(
            [tss_dynamic_overlap, tss_stable_overlap, tss_nonexpr_overlap], axis=0)
        data = data.merge(tss_data, how="left", left_on="atac_id", right_on="atac_id_tss")
        data = data.drop(["atac_id_tss"], axis=1)

        # pull in TSS dynamic traj
        tss_traj = pd.read_csv(
            args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"], sep="\t")
        tss_traj.columns = ["rna_cluster", "id"]
        data = data.merge(tss_traj, how="left", left_on="rna_id", right_on="id")
        data["rna_cluster"] = data["rna_cluster"].fillna(0).astype(int)
        data = data.drop("id", axis=1)

        # pull in histone mark data
        histone_dynamic_active = pd.read_csv(
            args.outputs["results"]["epigenome"]["dynamic"][
                "atac.dynamic.overlap_histones.mat"], sep="\t")
        histone_dynamic_active = histone_dynamic_active[["id", "H3K27ac", "H3K4me1"]]
        histone_dynamic_repressive = pd.read_csv(
            "{}/overlap_histone.inactive/ggr.epigenome.dynamic.overlaps.mat.txt.gz".format(
                args.outputs["results"]["epigenome"]["dynamic"]["dir"]), sep="\t")
        histone_dynamic_repressive = histone_dynamic_repressive[["id", "H3K27me3"]]
        histone_dynamic = histone_dynamic_active.merge(
            histone_dynamic_repressive, how="outer", left_on="id", right_on="id")
        histone_stable = pd.read_csv(
            args.outputs["results"]["epigenome"]["stable"][
                "atac.stable.overlap_histones.mat"], sep="\t")
        histone_stable = histone_stable[["id", "H3K27ac", "H3K4me1", "H3K27me3"]]
        histone_data = pd.concat([histone_dynamic, histone_stable], axis=0)
        data = data.merge(histone_data, how="left", left_on="atac_id", right_on="id")
        data = data.drop("id", axis=1)

        # MERGE: bring together the non-ATAC and ATAC data
        data = pd.concat([data, data_no_atac], axis=0)
        
        # making some into binary
        data["atac_dynamic"] = 0
        data.loc[data["atac_cluster"] > 0, "atac_dynamic"] = 8
        data.loc[data["atac_cluster"] == 0, "atac_dynamic"] = 7
        data["rna_present"] = 0
        data.loc[data["rna_pattern"] == "dynamic", "rna_present"] = 6
        data.loc[data["rna_pattern"] == "stable", "rna_present"] = 5
        data.loc[data["rna_pattern"] == "nonexpr", "rna_present"] = 4
        data["H3K27ac_present"] = (data["H3K27ac"] != 9).astype(int) * 3
        data["H3K4me1_present"] = (data["H3K4me1"] != 9).astype(int) * 2
        data["H3K27me3_present"] = (data["H3K27me3"] != 9).astype(int) * 1

        # save out file
        data.to_csv(
            summary_mat_file, sep="\t", header=True, index=False, compression="gzip")

        # now do another summary file that is aggregated for plotting
        plot_headers = [
            "atac_dynamic",
            "rna_present",
            "H3K27ac_present",
            "H3K4me1_present",
            "H3K27me3_present"]
        data = data[plot_headers]
        data["count"] = 1
        data = data.groupby(plot_headers).sum().reset_index()
        # TODO - how to best sort?
        
        for plot_header_idx in range(len(plot_headers)):
            plot_header = plot_headers[plot_header_idx]

            data_tmp = data[[plot_header, "count"]].copy()
            data_tmp.columns = ["fill", "count"]
            data_tmp["name"] = plot_header
            
            # save together
            if plot_header_idx == 0:
                data_melt = data_tmp.copy()
            else:
                data_melt = pd.concat([data_melt, data_tmp], axis=0)
                
        data_melt.to_csv(
            plot_mat_file, sep="\t", header=True, index=False, compression="gzip")

    # plot
    plot_file = "{}.pdf".format(summary_mat_file.split(".txt")[0])
    if not os.path.isfile(plot_file):
        r_cmd = "~/git/ggr-project/R/plot.epigenome_summary.R {} {}".format(
            plot_mat_file, plot_file)
        print r_cmd
        os.system(r_cmd)
    
    return args


def _run_epigenome_tss_workflow(args, prefix):
    """run tss centric analyses
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: look at distal v proximal regions")

    # assertions
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "proximal_regions"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    # inputs
    work_dir = args.outputs["results"]["epigenome"]["dir"]
    atac_categories = ["atac_dynamic", "atac_stable", "no_atac"]
    tss_categories = ["rna_dynamic", "rna_stable", "rna_non_expr"]
    dirnames = ["distal_v_proximal", "distal_v_proximal", "summary"]

    tss_dynamic_bed = args.outputs["data"]["tss.dynamic"]
    tss_stable_bed = args.outputs["data"]["tss.stable"]
    tss_nonexpr_bed = args.outputs["data"]["tss.non_expr"]
    tss_files = [tss_dynamic_bed, tss_stable_bed, tss_nonexpr_bed]
    
    # summary file to summary bed
    summary_mat_file = args.outputs["results"]["epigenome"]["summary"]["summary_mat"]
    summary_bed = "{}/{}.all_regions.for_matching.bed.gz".format(results_dir, prefix)
    if not os.path.isfile(summary_bed):
        summary_mat = pd.read_csv(summary_mat_file, sep="\t")
        summary_mat[["chrom", "start-stop"]] = summary_mat["atac_id"].str.split(":", n=2, expand=True)
        summary_mat[["start", "stop"]] = summary_mat["start-stop"].str.split("-", n=2, expand=True)
        summary_mat.to_csv(
            summary_bed, columns=["chrom", "start", "stop"],
            sep="\t", compression="gzip", header=False, index=False)
        
    # don't separate by ATAC category, just do TSS categories and plot out data sorted
    for tss_file_idx in range(len(tss_files)):
        tss_file = tss_files[tss_file_idx]
        tss_category = tss_categories[tss_file_idx]
        
        # overlap with summary mat file
        tss_overlap_file = "{}/{}.overlap_summary.mat.txt.gz".format(
            results_dir, os.path.basename(tss_file).split(".bed")[0])
        if not os.path.isfile(tss_overlap_file):
            bedtools_cmd = (
                "bedtools intersect -wao -a {} -b {} | "
                "gzip -c > {}").format(
                    tss_file, summary_bed, tss_overlap_file)
            print bedtools_cmd
            os.system(bedtools_cmd)

        # sort TSS file according to epigenetic/transcriptomic details
        tss_sorted_file = "{}/{}.sorted.bed.gz".format(
            results_dir, os.path.basename(tss_file).split(".bed")[0])
        if not os.path.isfile(tss_sorted_file):
            
            # (clean up) and match over to relevant data in summary mat
            overlap_mat = pd.read_csv(tss_overlap_file, sep="\t", header=None)
            overlap_mat = overlap_mat[overlap_mat[6] != "."]
            overlap_mat["atac_id"] = overlap_mat[6] + ":" + overlap_mat[7].map(str) + "-" + overlap_mat[8].map(str)
            summary_mat = pd.read_csv(summary_mat_file, sep="\t")
            tss_data = overlap_mat.merge(summary_mat, how="left", on="atac_id")

            # add in orig TSS files to get back to 1bp coords
            tss_orig_data = pd.read_csv(
                args.outputs["annotations"]["tss.pc.bed"], sep="\t", header=None)
            tss_columns = [
                "tss.chrom", "tss.start", "tss.stop", "tss.gene_id", "tss.score", "tss.strand"]
            tss_orig_data.columns = tss_columns
            tss_data = tss_data.merge(
                tss_orig_data, how="left", left_on=3, right_on="tss.gene_id")
            tss_data["atac_not_present"] = (tss_data["atac_cluster"] <= 0).astype(int)
            tss_data["H3K27ac_not_present"] = (tss_data["H3K27ac_present"] == 0).astype(int)
            tss_data["H3K4me1_not_present"] = (tss_data["H3K4me1_present"] == 0).astype(int)
            
            # sort and save out
            sort_columns = [
                "H3K27me3_present",
                "rna_cluster",
                "atac_not_present",
                "H3K27ac_not_present",
                "H3K4me1_not_present"
            ]
            tss_data = tss_data.sort_values(sort_columns, axis=0) # note: ascending sort
            tss_sorted = tss_data[tss_columns].drop_duplicates()
            tss_sorted.to_csv(
                tss_sorted_file,
                sep="\t", header=False, index=False, compression="gzip")
            
            # make a promoter data matrix to get the RNA vals
            rna_data_file = "{}/{}.promoter_rna_data.mat.txt.gz".format(
                results_dir, os.path.basename(tss_sorted_file).split(".bed")[0])
            rna_traj_data = pd.read_csv(
                args.outputs["data"]["rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat"],
                sep="\t")
            rna_traj_filt = tss_sorted.merge(
                rna_traj_data, how="left", left_on="tss.gene_id", right_index=True, sort=False)
            rna_traj_filt.to_csv(
                rna_data_file, sep="\t", header=True, index=False, compression="gzip")

            # get subsetted agg plots here

            # if dynamic TSS, split by RNA trajectories/H3K27me3 into 4 files
            if tss_categories[tss_file_idx] == "rna_dynamic":
                # progenitor genes, no H3K27me3
                tss_progenitor_file = "{}/{}.progenitor_genes.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_dynamic = tss_data[tss_data["rna_cluster"] <= 5]
                tss_dynamic = tss_dynamic[tss_dynamic["rna_cluster"] != 0]
                tss_dynamic = tss_dynamic[tss_dynamic["H3K27me3_present"] == 0]
                tss_dynamic = tss_dynamic.sort_values(sort_columns, axis=0)
                tss_dynamic = tss_dynamic[tss_columns].drop_duplicates()
                tss_dynamic.to_csv(
                    tss_progenitor_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_progenitor_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_progenitor_file, plot_prefix)

                # early diff
                tss_earlydiff_file = "{}/{}.earlydiff_genes.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_earlydiff = tss_data[tss_data["rna_cluster"] >=6]
                tss_earlydiff = tss_earlydiff[tss_earlydiff["rna_cluster"] >=7]
                tss_earlydiff = tss_earlydiff[tss_earlydiff["H3K27me3_present"] == 0]
                tss_earlydiff = tss_earlydiff.sort_values(sort_columns, axis=0)
                tss_earlydiff = tss_earlydiff[tss_columns].drop_duplicates()
                tss_earlydiff.to_csv(
                    tss_earlydiff_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_earlydiff_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_earlydiff_file, plot_prefix)

                # late diff
                tss_latediff_file = "{}/{}.latediff_genes.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_latediff = tss_data[tss_data["rna_cluster"] >=8]
                tss_latediff = tss_latediff[tss_latediff["H3K27me3_present"] == 0]
                tss_latediff = tss_latediff.sort_values(sort_columns, axis=0)
                tss_latediff = tss_latediff[tss_columns].drop_duplicates()
                tss_latediff.to_csv(
                    tss_latediff_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_latediff_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_latediff_file, plot_prefix)

                # H3K27me3 group
                tss_repr_file = "{}/{}.H3K27me3_genes.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_repr = tss_data[tss_data["H3K27me3_present"] == 1]
                tss_repr = tss_repr.sort_values(sort_columns, axis=0)
                tss_repr = tss_repr[tss_columns].drop_duplicates()
                tss_repr.to_csv(
                    tss_repr_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_repr_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_repr_file, plot_prefix)
            
            # if stable TSS, split H3K27me3 group and H3K27ac/H3K4me1 group
            if tss_categories[tss_file_idx] == "rna_stable":
                # activating marks
                tss_activemarks_file = "{}/{}.activating_marks.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_activating = tss_data[tss_data["H3K27me3_present"] == 0]
                tss_activating = tss_activating.sort_values(sort_columns, axis=0)
                tss_activating = tss_activating[tss_columns].drop_duplicates()
                tss_activating.to_csv(
                    tss_activemarks_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_activemarks_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_activemarks_file, plot_prefix)

                # repressive mark
                tss_repr_file = "{}/{}.repressing_marks.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_repr = tss_data[tss_data["H3K27me3_present"] == 1]
                tss_repr = tss_repr.sort_values(sort_columns, axis=0)
                tss_repr = tss_repr[tss_columns].drop_duplicates()
                tss_repr.to_csv(
                    tss_repr_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_repr_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_repr_file, plot_prefix)

            # if non_expr TSS, split H3K27me3 group and H3K27ac/H3K4me1 group
            if tss_categories[tss_file_idx] == "rna_non_expr":
                # activating marks
                tss_activemarks_file = "{}/{}.activating_marks.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_activating = tss_data[tss_data["H3K27me3_present"] == 0]
                tss_activating = tss_activating[tss_activating["H3K27ac_not_present"] == 0]
                tss_activating = tss_activating[tss_activating["H3K4me1_not_present"] == 0]
                tss_activating = tss_activating.sort_values(sort_columns, axis=0)
                tss_activating = tss_activating[tss_columns].drop_duplicates()
                tss_activating.to_csv(
                    tss_activemarks_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_activemarks_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_activemarks_file, plot_prefix)

                # repressive mark
                tss_repr_file = "{}/{}.repressing_marks.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_repr = tss_data[tss_data["H3K27me3_present"] == 1]
                tss_repr = tss_repr.sort_values(sort_columns, axis=0)
                tss_repr = tss_repr[tss_columns].drop_duplicates()
                tss_repr.to_csv(
                    tss_repr_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_repr_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_repr_file, plot_prefix)

                # and neither
                tss_no_file = "{}/{}.no_marks.bed.gz".format(
                    results_dir, os.path.basename(tss_file).split(".bed")[0])
                tss_no = tss_data[tss_data["H3K27me3_present"] == 0]
                tss_no = tss_no[tss_no["H3K27ac_not_present"] == 1]
                tss_no = tss_no[tss_no["H3K4me1_not_present"] == 1]
                tss_no = tss_no.sort_values(sort_columns, axis=0)
                tss_no = tss_no[tss_columns].drop_duplicates()
                tss_no.to_csv(
                    tss_no_file,
                    sep="\t", header=False, index=False, compression="gzip")
                plot_prefix = "{}/{}".format(
                    results_dir, os.path.basename(tss_no_file).split(".bed")[0])
                plot_profile_heatmaps_workflow(args, tss_no_file, plot_prefix)
            
        # and plot global ones
        plot_prefix = "{}/{}".format(
            results_dir, os.path.basename(tss_sorted_file).split(".bed")[0])
        plot_profile_heatmaps_workflow(args, tss_sorted_file, plot_prefix)

    # here, get distal regions and plot same plots based on ATAC summits
    summary_mat = pd.read_csv(summary_mat_file, sep="\t")
    summary_mat = summary_mat[summary_mat["atac_cluster"] > 0]
    summary_mat = summary_mat[summary_mat["rna_id"].isna()]
    summary_mat[["chrom", "start-stop"]] = summary_mat["atac_id"].str.split(":", n=2, expand=True)
    summary_mat[["start", "stop"]] = summary_mat["start-stop"].str.split("-", n=2, expand=True)
    distal_regions_file = "{}/{}.atac.distal.master.bed.gz".format(results_dir, prefix)
    summary_mat.to_csv(
        distal_regions_file, columns=["chrom", "start", "stop"],
        header=False, index=False, compression="gzip", sep="\t")

    # also split into groups: progenitor assoc, early diff assoc, late diff assoc
    timeseries_dir = args.outputs["results"]["atac"]["timepoint_region_dir"]
    atac_timeseries_files = sorted(
        glob.glob("{}/*narrowPeak.gz".format(timeseries_dir)))

    # progenitor assoc
    distal_prefix = "{}/{}.atac.distal.progenitor_assoc".format(results_dir, prefix)
    distal_regions_file = "{}.bed.gz".format(distal_prefix)
    if not os.path.isfile(distal_regions_file):
        distal_data = summary_mat[summary_mat["atac_cluster"] <= 6]
        distal_data.to_csv(
            distal_regions_file, columns=["chrom", "start", "stop"],
            header=False, index=False, compression="gzip", sep="\t")
    summit_file = "{}.summits.bed".format(distal_regions_file.split(".bed")[0])
    if not os.path.isfile(summit_file):
        get_best_summit(
            distal_regions_file,
            atac_timeseries_files,
            summit_file)
    plot_profile_heatmaps_workflow(args, summit_file, distal_prefix)

    # early diff assoc
    distal_prefix = "{}/{}.atac.distal.earlydiff_assoc".format(results_dir, prefix)
    distal_regions_file = "{}.bed.gz".format(distal_prefix)
    if not os.path.isfile(distal_regions_file):
        distal_data = summary_mat[summary_mat["atac_cluster"] >= 7]
        distal_data = distal_data[distal_data["atac_cluster"] <= 10]
        distal_data.to_csv(
            distal_regions_file, columns=["chrom", "start", "stop"],
            header=False, index=False, compression="gzip", sep="\t")
    summit_file = "{}.summits.bed".format(distal_regions_file.split(".bed")[0])
    if not os.path.isfile(summit_file):
        get_best_summit(
            distal_regions_file,
            atac_timeseries_files,
            summit_file)
    plot_profile_heatmaps_workflow(args, summit_file, distal_prefix)

    # late diff assoc
    distal_prefix = "{}/{}.atac.distal.latediff_assoc".format(results_dir, prefix)
    distal_regions_file = "{}.bed.gz".format(distal_prefix)
    if not os.path.isfile(distal_regions_file):
        distal_data = summary_mat[summary_mat["atac_cluster"] >= 11]
        distal_data = distal_data[distal_data["atac_cluster"] <= 14]
        distal_data.to_csv(
            distal_regions_file, columns=["chrom", "start", "stop"],
            header=False, index=False, compression="gzip", sep="\t")
    summit_file = "{}.summits.bed".format(distal_regions_file.split(".bed")[0])
    if not os.path.isfile(summit_file):
        get_best_summit(
            distal_regions_file,
            atac_timeseries_files,
            summit_file)
    plot_profile_heatmaps_workflow(args, summit_file, distal_prefix)

    
    quit()

        
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

    # -----------------------------------------
    # ANALYSIS - condense all this information into dynamic chromatin states
    # inputs: epigenome outputs
    # outputs: chromatin states heatmap
    # -----------------------------------------
    args = run_chromatin_states_workflow(
        args, "{}.chromatin_states".format(prefix))

    # -----------------------------------------
    # ANALYSIS - look at distal v proximal regions
    # inputs: epigenome outputs
    # outputs: analysis of prom v enh
    # -----------------------------------------
    args = run_distal_v_proximal_workflow(
        args, "{}.distal_v_proximal".format(prefix))
    
    # -----------------------------------------
    # ANALYSIS - look at distal v proximal regions
    # inputs: epigenome outputs
    # outputs: analysis of prom v enh
    # -----------------------------------------
    args = _run_epigenome_summary_workflow(
        args, "{}.summary".format(prefix))

    # -----------------------------------------
    # ANALYSIS - look at proximal regions only ie TSS
    # inputs: TSS files
    # outputs: plots
    # -----------------------------------------
    args = _run_epigenome_tss_workflow(
        args, "{}.tss".format(prefix))

    quit()

    # TODO also look at distal no-ATAC regions
    # maybe just homer analysis?
    
    
    # -----------------------------------------
    # ANALYSIS - look at signal spreading
    # inputs: epigenome outputs
    # outputs: signal spreading results
    # -----------------------------------------
    args = run_histone_spreading_workflow(
        args, "{}.histone_spread".format(prefix))
    
    return args
