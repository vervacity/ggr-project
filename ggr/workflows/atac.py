"""ggr workflow for atac analyses
"""

import os
import glob
import signal
import logging

#import pandas as pd

from ggr.util.utils import run_shell_cmd
from ggr.util.utils import parallel_copy

#from ggr.util.parallelize import setup_multiprocessing_queue
#from ggr.util.parallelize import run_in_parallel

#from ggr.timeseries import get_consistent_dpgp_trajectories
#from ggr.timeseries import merge_cluster_files

from ggr.util.bed_utils import merge_regions
from ggr.util.bed_utils import id_to_bed
from ggr.util.filtering import filter_for_ids

from ggr.analyses.counting import make_count_matrix
#from ggr.analyses.counting import split_count_matrix_by_replicate

from ggr.workflows.timeseries import run_timeseries_workflow

#from ggr.util.bioinformatics import make_deeptools_heatmap


def runall(args, prefix):
    """all workflows for atac-seq data
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: run atac analyses")

    # set up inputs
    inputs = args.inputs["atac"][args.cluster]
    
    # assertions
    assert inputs.get("data_dir") is not None
    assert inputs.get("idr_peak_glob") is not None

    # set up data
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]
    
    results_dirname = "atac"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]

    # -------------------------------------------
    # ANALYSIS 0 - download timepoint BED files
    # input: peak files
    # output: same peak files, new location
    # -------------------------------------------
    logger.info("ANALYSIS: copy idr peak files to new dir")
    atac_peak_files = sorted(
        glob.glob('{0}/{1}'.format(
            inputs['data_dir'],
            inputs['idr_peak_glob'])))
    timepoint_dir = "{}/peaks.timepoints".format(results_dir)
    out_results["timepoint_region_dir"] = timepoint_dir
    timepoints_files = glob.glob("{}/*.narrowPeak.gz".format(
        timepoint_dir))
    if len(timepoints_files) != len(atac_peak_files):
        parallel_copy(
            atac_peak_files,
            timepoint_dir)
    timepoints_files = glob.glob("{}/*.narrowPeak.gz".format(
        timepoint_dir))
    
    # -------------------------------------------
    # ANALYSIS 1 - generate master regions file
    # input: peak files
    # output: master peak file (BED)
    # -------------------------------------------
    logger.info("ANALYSIS: generate master regions file")
    master_regions_key = "atac.master.bed"
    out_data[master_regions_key] = '{0}/{1}.idr.master.bed.gz'.format(
        data_dir, prefix)
    if not os.path.isfile(out_data[master_regions_key]):
        merge_regions(timepoints_files, out_results[master_regions_key])

    # -------------------------------------------
    # ANALYSIS 2 - get read counts in these regions
    # input: master regions, read files (BED/tagAlign format)
    # output: matrix of read counts per region
    # -------------------------------------------
    logger.info("ANALYSIS: get read counts per region")
    adjustment = "ends"
    counts_key = "atac.counts.mat"
    out_data[counts_key] = '{0}/{1}.{2}.counts.mat.txt.gz'.format(
        data_dir, prefix, adjustment)
    if not os.path.isfile(out_data[counts_key]):
        # glob tagalign files
        atac_tagalign_files = sorted(
            glob.glob('{0}/{1}'.format(
                inputs["data_dir"],
                inputs["tagalign_glob"])))
        # remove files from media influenced timepoints
        for timepoint_string in args.inputs["params"]["media_timepoints"]:
            atac_tagalign_files = [
                filename for filename in atac_tagalign_files
                if timepoint_string not in filename]
        # make count matrix
        make_count_matrix(
            out_data[master_regions_key],
            atac_tagalign_files,
            out_data[counts_key],
            "ATAC",
            adjustment=adjustment,
            tmp_dir=results_dir)

    # -------------------------------------------
    # ANALYSIS 3 - run timeseries analysis on these regions
    # input: count matrix of regions
    # output: region trajectories
    # -------------------------------------------
    args = run_timeseries_workflow(
        args,
        prefix,
        datatype_key="atac",
        mat_key=counts_key)

    quit()
    
    # extract dynamic BED and stable BED files
    # also cluster BED files

    # also set up stable mat for ATAC pooled

    # bioinformatics: GREAT, HOMER



    # now plot everything - pdfs, to put into illustrator

    
    

    
    
    
    quit()

    
    # -------------------------------------------
    # ANALYSIS 3 - run DESeq2 for all pairs of timepoints
    # input: matrix of read counts per region
    # output: dynamic IDs
    # -------------------------------------------
    logger.info("ANALYSIS: run deseq2 in timeseries mode...")
    args.atac['dynamic_ids.list'] = '{0}/{1}.dynamic.ids.txt.gz'.format(
        args.folders['atac_timeseries_dir'], atac_prefix)
    if not os.path.isfile(args.atac['dynamic_ids']):
	run_timeseries_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3} full".format(
            args.atac['counts'],
            '{0}/{1}'.format(
                args.folders['atac_deseq2_dir'],
                atac_prefix),
            args.params['atac_deseq2_fdr'],
            args.atac['dynamic_ids'])
        run_shell_cmd(run_timeseries_deseq2)

    # -------------------------------------------
    # ANALYSIS 4 - normalization of count files for timeseries analysis
    # input: matrix of read counts per region with reps
    # output: normalized matrices (pooled, rep1, rep2)
    # -------------------------------------------
    # 4) trajectories: make normalized rep1, rep2, pooled files (for DP_GP analysis)
    counts_prefix = args.atac["counts"].split(".mat")[0]
    args.atac["counts_rep1"] = "{}.rep1.mat.txt.gz".format(counts_prefix)
    args.atac["counts_rep2"] = "{}.rep2.mat.txt.gz".format(counts_prefix)
    args.atac["counts_pooled"] = "{}.pooled.mat.txt.gz".format(counts_prefix)

    if not os.path.isfile(args.atac["counts_pooled"]):
        logger.info("ATAC: Separating count mat into reps/pooled...")
        split_count_matrix_by_replicate(
            args.atac["counts"],
            args.atac["counts_rep1"],
            args.atac["counts_rep2"],
            args.atac["counts_pooled"])


    if not os.path.isfile("{}.rlog.mat.txt.gz".format(
            args.atac["counts_pooled"].split(".mat")[0])):
        logger.info("ATAC: normalizing count files...")
        run_rlogs = "normalize_count_mats.R {0} {1} {2} {3}".format(
            args.atac["counts"],
            args.atac["counts_rep1"],
            args.atac["counts_rep2"],
            args.atac["counts_pooled"])
        print run_rlogs
        run_shell_cmd(run_rlogs)
    
    # for each of the counts files:
    logger.info("ATAC: filtering dynamic/stable...")
    counts_files_to_normalize = ["counts_rep1", "counts_rep2", "counts_pooled"]
    for counts_handle in counts_files_to_normalize:

        rlog_handle = "{}_rlog".format(counts_handle)
        args.atac[rlog_handle] = "{}.rlog.mat.txt.gz".format(
            args.atac[counts_handle].split('.mat')[0])

        # 6) filter appropriate count matrices for the dynamic and stable ones
        dynamic_handle = "{}_dynamic".format(rlog_handle)
        args.atac[dynamic_handle] = "{}.dynamic.mat.txt.gz".format(
            args.atac[rlog_handle].split(".mat")[0])
        if not os.path.isfile(args.atac[dynamic_handle]):
            filter_for_ids(args.atac[rlog_handle],
                           args.atac["dynamic_ids"],
                           args.atac[dynamic_handle])

        # only filter stable for the POOLED set
        if "pooled" in counts_handle:
            stable_handle = "{}_stable".format(rlog_handle)
            args.atac[stable_handle] = "{}.stable.mat.txt.gz".format(
                args.atac[rlog_handle].split('.mat')[0])
            if not os.path.isfile(args.atac[stable_handle]):
                filter_for_ids(args.atac[rlog_handle],
                               args.atac["dynamic_ids"],
                               args.atac[stable_handle],
                               opposite=True)
                
            # also make a stable BED file for downstream analyses
            stable_bed_handle = "{}_bed".format(stable_handle)
            args.atac[stable_bed_handle] = "{}.bed.gz".format(
                args.atac[stable_handle].split(".mat")[0])
            if not os.path.isfile(args.atac[stable_bed_handle]):
                make_bed = (
                    "zcat {} | "
                    "awk -F '\t' '{{ print $1 }}' | "
                    "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
                    "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
                    "grep -v d0 | "
                    "awk 'NF' | "
                    "gzip -c > {}").format(
                        args.atac[stable_handle],
                        args.atac[stable_bed_handle])
                print make_bed
                run_shell_cmd(make_bed)

            # also make a dynamic BED file for downstream analyses
            dynamic_bed_handle = "{}_bed".format(dynamic_handle)
            args.atac[dynamic_bed_handle] = "{}.bed.gz".format(
                args.atac[dynamic_handle].split(".mat")[0])
            if not os.path.isfile(args.atac[dynamic_bed_handle]):
                make_bed = (
                    "zcat {} | "
                    "awk -F '\t' '{{ print $1 }}' | "
                    "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
                    "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
                    "grep -v d0 | "
                    "awk 'NF' | "
                    "gzip -c > {}").format(
                        args.atac[dynamic_handle],
                        args.atac[dynamic_bed_handle])
                print make_bed
                run_shell_cmd(make_bed)
                
    # 7) Run trajectories (on dynamic set) using DP_GP clustering
    args.atac["consistent_clusters"] = get_consistent_dpgp_trajectories(
        args.atac["counts_rep1_rlog_dynamic"],
        args.atac["counts_rep2_rlog_dynamic"],
        args.atac["counts_pooled_rlog_dynamic"],
        args.folders["atac_dp-gp_dir"],
        atac_prefix,
        raw_cluster_min_size=args.params["raw_cluster_min_size"],
        raw_cluster_reject_null_ci_interval=args.params["raw_cluster_reject_null_ci_interval"],
        rep_to_cluster_ci_interval=args.params["rep_to_cluster_ci_interval"],
        rep_to_cluster_corr_cutoff=args.params["rep_to_cluster_corr_cutoff"],
        epsilon=args.params["rep_to_cluster_epsilon"])

    # 8) reorder (ie renumber) the SOFT clusters (numerically) by order of hclust and make sure to propagate to hard clusters
    args.atac["final_hard_clusters"] = "{}/{}.clusters.hard.renumbered.all.txt.gz".format(
        args.folders["atac_dp-gp_final_dir"],
        atac_prefix)
    if not os.path.isfile(args.atac["final_hard_clusters"]):
        reorder_clusters = "reorder_soft_clusters_w_hclust.R {0} {1}/hard/*hard.all*gz {1}/soft/*soft*gz".format(
            args.folders["atac_dp-gp_final_dir"],
            args.folders["atac_dp-gp_dir"])
        print reorder_clusters
        run_shell_cmd(reorder_clusters)
    args.atac["final_soft_clusters"] = sorted(
        glob.glob("{}/*soft*.gz".format(args.folders["atac_dp-gp_final_dir"])))

    # 9) make BED files from soft cluster files
    atac_cluster_beds = glob.glob("{}/*bed.gz".format(args.folders["atac_dp-gp_final_bed_dir"]))
    if len(atac_cluster_beds) == 0:
        atac_cluster_id_files = glob.glob("{}/*soft.txt.gz".format(args.folders["atac_dp-gp_final_dir"]))
        for atac_cluster_id_file in atac_cluster_id_files:
            atac_cluster_bed_file = "{}/{}.bed.gz".format(
                args.folders["atac_dp-gp_final_bed_dir"],
                os.path.basename(atac_cluster_id_file).split('.txt')[0])
            id_to_bed(atac_cluster_id_file, atac_cluster_bed_file, sort=True)

    return args
