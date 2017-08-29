"""Contains functions to run ATAC-specific analyses
"""

import os
import glob
import signal
import logging

import pandas as pd

from ggr.timeseries import get_consistent_dpgp_trajectories
from ggr.timeseries import merge_cluster_files

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel

from ggr.util.bed_utils import merge_regions
from ggr.util.bed_utils import id_to_bed
from ggr.util.filtering import filter_for_ids
from ggr.util.counting import make_count_matrix
from ggr.util.counting import split_count_matrix_by_replicate

from ggr.util.bioinformatics import make_deeptools_heatmap


def run(args):
    """Run all ATAC analyses
    """
    atac_prefix = 'ggr.atac'

    # 0) download peak files into a folder (timepoint labels)
    timepoints_files = glob.glob("{}/*.narrowPeak.gz".format(
        args.folders["atac_timepoints_bed_dir"]))
    if len(timepoints_files) == 0:
        atac_peak_files = sorted(
            glob.glob('{0}/{1}'.format(
                args.atac['data_dir'],
                args.atac['idr_peak_glob'])))
        for atac_peak_file in atac_peak_files:
            copy = "cp {} {}".format(
                atac_peak_file,
                args.folders["atac_timepoints_bed_dir"])
            print copy
            os.system(copy)
    
    # 1) generate a master regions file
    logging.info("ATAC: Generating master file...")
    args.atac['master_bed'] = '{0}/{1}.idr.master.bed.gz'.format(
        args.folders["data_dir"], atac_prefix)
    if not os.path.isfile(args.atac["master_bed"]):
        atac_peak_files = sorted(
            glob.glob('{0}/{1}'.format(
                args.atac['data_dir'],
                args.atac['idr_peak_glob'])))
        merge_regions(atac_peak_files, args.atac['master_bed'])
    
    # 2) get read counts in these regions (for DESeq2 analysis)
    logging.info("ATAC: Generating read counts per master region...")
    adjustment = "ends"
    args.atac['counts'] = '{0}/{1}.{2}.counts.mat.txt.gz'.format(
        args.folders['data_dir'], atac_prefix, adjustment)
    if not os.path.isfile(args.atac['counts']):
        # glob tagalign files
        atac_tagalign_files = sorted(
            glob.glob('{0}/{1}'.format(
                args.atac['data_dir'],
                args.atac['tagalign_glob'])))
        # remove files from media influenced timepoints
        for timepoint_string in args.misc["media_timepoints"]:
            atac_tagalign_files = [
                filename for filename in atac_tagalign_files
                if timepoint_string not in filename]
        # make count matrix
        make_count_matrix(
            args.atac['master_bed'],
            atac_tagalign_files,
            args.atac['counts'],
            "ATAC",
            adjustment="ends",
            tmp_dir=args.folders["atac_dir"])

    # 3) with read counts run DESeq2 as timeseries mode
    logging.info("ATAC: Running DESeq2 on all pairs of timepoints...")
    args.atac['dynamic_ids'] = '{0}/{1}.dynamic.ids.txt.gz'.format(
        args.folders['atac_timeseries_dir'], atac_prefix)
    if not os.path.isfile(args.atac['dynamic_ids']):
	run_timeseries_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3} full".format(
            args.atac['counts'],
            '{0}/{1}'.format(
                args.folders['atac_deseq2_dir'],
                atac_prefix),
            args.params['atac_deseq2_fdr'],
            args.atac['dynamic_ids'])
        print run_timeseries_deseq2
        os.system(run_timeseries_deseq2)

    # 4) trajectories: make normalized rep1, rep2, pooled files (for DP_GP analysis)
    counts_prefix = args.atac["counts"].split(".mat")[0]
    args.atac["counts_rep1"] = "{}.rep1.mat.txt.gz".format(counts_prefix)
    args.atac["counts_rep2"] = "{}.rep2.mat.txt.gz".format(counts_prefix)
    args.atac["counts_pooled"] = "{}.pooled.mat.txt.gz".format(counts_prefix)
    
    logging.info("ATAC: Separating count mat into reps/pooled...")
    if not os.path.isfile(args.atac["counts_pooled"]):
        split_count_matrix_by_replicate(
            args.atac["counts"],
            args.atac["counts_rep1"],
            args.atac["counts_rep2"],
            args.atac["counts_pooled"])

    logging.info("ATAC: normalizing count files...")
    if not os.path.isfile("{}.rlog.mat.txt.gz".format(
            args.atac["counts_pooled"].split(".mat")[0])):
        run_rlogs = "normalize_count_mats.R {0} {1} {2} {3}".format(
            args.atac["counts"],
            args.atac["counts_rep1"],
            args.atac["counts_rep2"],
            args.atac["counts_pooled"])
        print run_rlogs
        os.system(run_rlogs)
    
    # for each of the counts files:
    logging.info("ATAC: filtering dynamic/stable...")
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
                os.system(make_bed)

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
                os.system(make_bed)
                
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
        rep_to_cluster_corr_cutoff=args.params["rep_to_cluster_corr_cutoff"])

    # 8) reorder (ie renumber) the SOFT clusters (numerically) by order of hclust and make sure to propagate to hard clusters
    args.atac["final_hard_clusters"] = "{}/{}.clusters.hard.renumbered.txt.gz".format(
        args.folders["atac_dp-gp_final_dir"],
        atac_prefix)
    if not os.path.isfile(args.atac["final_hard_clusters"]):
        reorder_clusters = "reorder_soft_clusters_w_hclust.R {0} {1}/soft/*hard*gz {1}/soft/*soft*gz".format(
            args.folders["atac_dp-gp_final_dir"],
            args.folders["atac_dp-gp_dir"])
        print reorder_clusters
        os.system(reorder_clusters)
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
