"""Contains code to run histone analyses
"""

import os
import gzip
import glob
import logging
import signal

import pandas as pd

from ggr.util.utils import add_folder
from ggr.util.counting import make_count_matrix
from ggr.util.counting import split_count_matrix_by_replicate
from ggr.util.bed_utils import merge_regions
from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel

from ggr.util.diff import join_diff_region_lists_to_mat


def setup_master_regions(args, histone, method="atac_midpoint", filter_w_overlap=True):
    """modularize the master region setup for easy switching
    1) use naive overlap set
    2) go atac-centric with fixed extend length from atac midpoint
    """
    # set up naive overlap based master regions
    logging.info("HISTONE: {}: Generating master regions...".format(histone))
    args.chipseq["histones"][histone]["overlap_master_regions"] = "{0}/ggr.{1}.overlap.master.bed.gz".format(
        args.folders["data_dir"], histone)
    if not os.path.isfile(args.chipseq["histones"][histone]["overlap_master_regions"]):
        histone_overlap_files = sorted(
            glob.glob("{0}/{1}".format(
                args.chipseq["data_dir"], args.chipseq["histones"][histone]["overlap_glob"])))
        logging.info("Master regions using: {}".format(" ".join(histone_overlap_files)))
        merge_regions(histone_overlap_files, args.chipseq["histones"][histone]["overlap_master_regions"])
    
    if method == "atac_midpoint":
        # If centered on ATAC region midpoint, then extract the midpoint and then extend out with bedtools slop
        args.atac["master_slop_bed"] = "{}.slop_{}bp.bed.gz".format(
            args.atac["master_bed"].split(".bed")[0],
            args.params["histones"][histone]["overlap_extend_len"])
        if not os.path.isfile(args.atac["master_slop_bed"]):
            slop_bed = (
                "zcat {0} | "
                "awk -F '\t' 'BEGIN{{OFS=\"\t\"}} "
                "{{ midpoint=$2+int(($3-$2)/2); "
                "$2=midpoint; $3=midpoint+1; print }}' | "
                "bedtools slop -i stdin -g {1} -b {2} | "
                "gzip -c > {3}").format(
                    args.atac["master_bed"],
                    args.annot["chromsizes"],
                    args.params["histones"][histone]["overlap_extend_len"],
                    args.atac["master_slop_bed"])
            print slop_bed
            os.system(slop_bed)

        if filter_w_overlap:
            # now intersect ATAC with the naive overlap files and only keep region if has an overlap
            args.chipseq["histones"][histone]["master_slop_marked_bed"] = "{}.{}-marked.bed.gz".format(
                args.atac["master_slop_bed"].split(".bed")[0],
                histone)
            if not os.path.isfile(args.chipseq["histones"][histone]["master_slop_marked_bed"]):
                keep_marked = (
                    "bedtools intersect -u -a {0} -b {1} | "
                    "gzip -c > {2}").format(
                        args.atac["master_slop_bed"],
                        args.chipseq["histones"][histone]["overlap_master_regions"],
                        args.chipseq["histones"][histone]["master_slop_marked_bed"])
                print keep_marked
                os.system(keep_marked)
            master_regions = args.chipseq["histones"][histone]["master_slop_marked_bed"]
            
        else:
            master_regions = args.atac["master_slop_bed"]

    elif method == "naive_overlap":
        # If naive overlap, don't do anything extra - already generated the master file
        master_regions = args.chipseq["histones"][histone]["overlap_master_regions"]

    else:
        raise Exception("non existent master regions method!")
            
    return master_regions


def run(args):
    """Pipeline for histone analyses
    """
    # set up histones
    histones = args.chipseq["histones"]["ordered_names"]
    
    for histone in histones:

        logging.info("HISTONE: Working on {}".format(histone))
        histone_prefix = "ggr.{}".format(histone)
        histone_dir = args.folders["{}_dir".format(histone)]

        args.chipseq["histones"][histone]["master_regions"] = setup_master_regions(
            args, histone, method=args.params["histone_master_region_method"], filter_w_overlap=True)

        # 2) get counts of midpoint of fragments in these regions
        logging.info("HISTONE: {}: Making counts matrix...".format(histone))
        args.chipseq["histones"][histone]["counts"] = "{0}/{1}.midpoints.counts.mat.txt.gz".format(
            args.folders["data_dir"], histone_prefix)
        if not os.path.isfile(args.chipseq["histones"][histone]["counts"]):
            bedpe_file_list = sorted(glob.glob("{}/{}".format(
                args.chipseq["data_dir"],
                args.chipseq["histones"][histone]["bedpe_glob"])))
            assert len(bedpe_file_list) == 6
            # Make count matrix
            make_count_matrix(
                args.chipseq["histones"][histone]["master_regions"],
                bedpe_file_list,
                args.chipseq["histones"][histone]["counts"],
                histone,
                adjustment="midpoints",
                tmp_dir=histone_dir)

        # 3) run DESeq on forward sequential pairs of timepoints
        logging.info("HISTONE: {}: Running DESeq2 analysis...".format(histone))
        histone_diff_dir = args.folders["{}_diff_dir".format(histone)]
        args.chipseq["histones"][histone]["dynamic_ids"] = "{0}/ggr.{1}.dynamic.ids.txt.gz".format(
            histone_diff_dir, histone)
        if not os.path.isfile(args.chipseq["histones"][histone]["dynamic_ids"]):
            run_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3} sequential".format(
                args.chipseq["histones"][histone]["counts"],
                "{0}/{1}".format(histone_diff_dir, histone_prefix),
                args.params["histone_fdr"],
                args.chipseq["histones"][histone]["dynamic_ids"])
            print run_deseq2
            os.system(run_deseq2)
            
        # 4) run rlog to have standardized matrix (to viz and do QC)
        logging.info("HISTONE: {}: splitting count mat into reps/pooled...".format(histone))
        counts_prefix = args.chipseq["histones"][histone]["counts"].split(".mat")[0]
        args.chipseq["histones"][histone]["counts_rep1"] = "{}.rep1.mat.txt.gz".format(
            counts_prefix)
        args.chipseq["histones"][histone]["counts_rep2"] = "{}.rep2.mat.txt.gz".format(
            counts_prefix)
        args.chipseq["histones"][histone]["counts_pooled"] = "{}.pooled.mat.txt.gz".format(
            counts_prefix)
        if not os.path.isfile(args.chipseq["histones"][histone]["counts_pooled"]):
            split_count_matrix_by_replicate(
                args.chipseq["histones"][histone]["counts"],
                args.chipseq["histones"][histone]["counts_rep1"],
                args.chipseq["histones"][histone]["counts_rep2"],
                args.chipseq["histones"][histone]["counts_pooled"])

        logging.info("HISTONE: {}: normalizing count_files...".format(histone))
        args.chipseq["histones"][histone]["counts_pooled_rlog"] = "{}.pooled.rlog.mat.txt.gz".format(
            counts_prefix)
        if not os.path.isfile(args.chipseq["histones"][histone]["counts_pooled_rlog"]):
            run_rlogs = "normalize_count_mats.R {0} {1}".format(
                args.chipseq["histones"][histone]["counts"],
                args.chipseq["histones"][histone]["counts_pooled"])
            print run_rlogs
            os.system(run_rlogs)
            
        # 5) enumerate trajectories (9 groups)
        sig_up_files = sorted(glob.glob("{}/*sigResultsUp.txt.gz".format(histone_diff_dir)))
        sig_down_files = sorted(glob.glob("{}/*sigResultsDown.txt.gz".format(histone_diff_dir)))
        sig_up_down_pairs = zip(sig_up_files, sig_down_files)
        args.chipseq["histones"][histone]["sig_mat"] = "{}/{}.deseq2.mat.txt.gz".format(
            histone_diff_dir, histone_prefix)
        if not os.path.isfile(args.chipseq["histones"][histone]["sig_mat"]):
            join_diff_region_lists_to_mat(
                args.chipseq["histones"][histone]["master_regions"],
                sig_up_down_pairs,
                args.chipseq["histones"][histone]["sig_mat"])
            
        args.chipseq["histones"][histone]["clusters_bed"] = "{}/{}.deseq2.clusters.bed.gz".format(
            histone_diff_dir, histone_prefix)
        if not os.path.isfile(args.chipseq["histones"][histone]["clusters_bed"]):
            make_bed = (
                "zcat {0} | "
                "awk -F '\t' '{{ print $1\"\t\"$4 }}' | "
                "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
                "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
                "grep -v regions | "
                "sort -k1,1 -k2,2n | "
                "gzip -c > {1}").format(
                    args.chipseq["histones"][histone]["sig_mat"],
                    args.chipseq["histones"][histone]["clusters_bed"])
            print make_bed
            os.system(make_bed)
            
    return args
