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



def run_histone_idr(
        args,
        histone,
        master_file,
        bedpe_file_list,
        out_dir,
        prefix,
        fdr=0.10):
    """Do a version of IDR on histone files
    NOTES: do on naive overlap file
    """
    # sort master file
    master_ordered_file = "{}/{}.master.bed.gz".format(out_dir, prefix)
    if not os.path.isfile(master_ordered_file):
        order_bed = (
            "zcat {0} | "
            "awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$3 }}' | "
            "sort -k1,1 -k2,2n | "
            "gzip -c > {1}").format(master_file, master_ordered_file)
        print order_bed
        os.system(order_bed)
    
    # get counts of midpoint of fragments in these regions
    counts_file = "{}/{}.midpoints.counts.mat.txt.gz".format(out_dir, prefix)
    if not os.path.isfile(counts_file):
        make_count_matrix(
            master_ordered_file,
            bedpe_file_list,
            counts_file,
            histone,
            tmp_dir=out_dir)

    # run rlog to normalize
    counts_rlog_file = "{}.rlog.mat.txt.gz".format(counts_file.split(".mat")[0])
    if not os.path.isfile(counts_rlog_file):
        run_rlogs = "normalize_count_mat.R {0} {1}".format(
            counts_file,
            counts_rlog_file)
        print run_rlogs
        os.system(run_rlogs)

    # separate into 2 replicate files
    rep1_out_file = "{}.rep1.bed.gz".format(counts_rlog_file.split('.mat')[0])
    rep2_out_file = "{}.rep2.bed.gz".format(counts_rlog_file.split('.mat')[0])
    if not os.path.isfile(rep2_out_file):
        with gzip.open(counts_rlog_file, 'r') as fp:
            with gzip.open(rep1_out_file, 'w') as out1:
                with gzip.open(rep2_out_file, 'w') as out2:

                    for line in fp:
                        
                        if "b1" in line:
                            continue

                        fields = line.strip().split('\t')
                        chrom = fields[0].split(':')[0]
                        start = fields[0].split(':')[1].split('-')[0]
                        stop = fields[0].split('-')[1]

                        out1.write("{0}\t{1}\t{2}\t{0}:{1}-{2}\t{3}\t.\t{1}\t{2}\t0\n".format(chrom, start, stop, fields[1]))
                        out2.write("{0}\t{1}\t{2}\t{0}:{1}-{2}\t{3}\t.\t{1}\t{2}\t0\n".format(chrom, start, stop, fields[2]))

    # call IDR on these rep files
    # NOTE since idr is python3, need to step out of code to run
    # command is kept here for reference
    args.chipseq["histones"][histone]["idr_bed"] = "{}/{}.idr.bed.gz".format(out_dir, prefix)
    idr_out_log_file = "{}/{}.idr.log".format(out_dir, prefix)

    if not os.path.isfile(args.chipseq["histones"][histone]["idr_bed"]):
        idr = (
            "idr --samples {0} {1} --input-file-type bed "
            "--rank score "
            "--soft-idr-threshold {2} "
            "--output-file {3} "
            "--log-output-file {4} "
            "--plot").format(
                rep1_out_file, rep2_out_file,
                0.05,
                args.chipseq["histones"][histone]["idr_bed"].split(".gz")[0],
                idr_out_log_file)
        print idr

    return args


def run(args):
    """Pipeline for histone analyses
    """
        
    # set up histones
    histones = args.chipseq["histones"].keys()
    days = ["d0", "d3", "d6"]
    
    for histone in histones:

        logging.info("HISTONE: Working on {}".format(histone))
        histone_prefix = "ggr.{}".format(histone)
        histone_dir = args.folders["{}_dir".format(histone)]

        # set up file lists
        histone_overlap_files = sorted(
            glob.glob("{0}/{1}".format(
                args.chipseq["data_dir"], args.chipseq["histones"][histone]["overlap_glob"])))
        bedpe_file_list = sorted(glob.glob("{}/{}".format(
            args.chipseq["data_dir"],
            args.chipseq["histones"][histone]["bedpe_glob"])))
        assert len(bedpe_file_list) == 6


        # make the extended ATAC master regions as key regions for histones
        # TODO go from midpoint outward
        args.atac["master_slop_bed"] = "{}.slop_{}bp.bed.gz".format(
            args.atac["master_bed"].split(".bed")[0],
            args.params["histones"][histone]["atac_extend_len"])
        if not os.path.isfile(args.atac["master_slop_bed"]):
            slop_bed = (
                "zcat {0} | "
                "awk -F '\t' 'BEGIN{{OFS=\"\t\"}} {{ midpoint=$2+int(($3-$2)/2); $2=midpoint; $3=midpoint+1; print }}' | "
                "bedtools slop -i stdin -g {1} -b {2} | "
                "gzip -c > {3}").format(
                    args.atac["master_bed"],
                    args.annot["chromsizes"],
                    args.params["histones"][histone]["atac_extend_len"],
                    args.atac["master_slop_bed"])
            print slop_bed
            os.system(slop_bed)

        # TODO instead of IDR, use ATAC regions as master list.
        # first, only mark a region if it overlaps a naive overlap peak (master regions file)
        # then, with those regions, normalize to {extendlen} size and count fragments
        # and continue on as before
        # make sure to switch the nearest histone to look for nearest
        
        # 1) generate master regions file
        logging.info("HISTONE: {}: Generating master regions...".format(histone))
        args.chipseq["histones"][histone]["master_regions"] = "{0}/ggr.{1}.overlap.master.bed.gz".format(
            args.folders["data_dir"], histone)
        if not os.path.isfile(args.chipseq["histones"][histone]["master_regions"]):
            histone_overlap_files = sorted(
                glob.glob("{0}/{1}".format(
                    args.chipseq["data_dir"], args.chipseq["histones"][histone]["overlap_glob"])))
            logging.info("Master regions using: {}".format(" ".join(histone_overlap_files)))
            merge_regions(histone_overlap_files, args.chipseq["histones"][histone]["master_regions"])

        # now intersect ATAC with the naive overlap files and only keep region if has an overlap
        args.chipseq["histones"][histone]["master_slop_marked_bed"] = "{}.{}-marked.bed.gz".format(
            args.atac["master_slop_bed"].split(".bed")[0],
            histone)
        if not os.path.isfile(args.chipseq["histones"][histone]["master_slop_marked_bed"]):
            keep_marked = (
                "bedtools intersect -u -a {0} -b {1} | "
                "gzip -c > {2}").format(
                    args.atac["master_slop_bed"],
                    args.chipseq["histones"][histone]["master_regions"],
                    args.chipseq["histones"][histone]["master_slop_marked_bed"])
            print keep_marked
            os.system(keep_marked)

        # 2) get counts of midpoint of fragments in these regions
        logging.info("HISTONE: {}: Making counts matrix...".format(histone))
        args.chipseq["histones"][histone]["counts"] = "{0}/{1}.midpoints.counts.mat.txt.gz".format(
            args.folders["data_dir"], histone_prefix)
        if not os.path.isfile(args.chipseq["histones"][histone]["counts"]):
            # Make count matrix
            make_count_matrix(
                args.chipseq["histones"][histone]["master_slop_marked_bed"],
                #args.atac["master_slop_bed"], # NOTE CHANGED HERE
                bedpe_file_list,
                args.chipseq["histones"][histone]["counts"],
                histone,
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
                args.chipseq["histones"][histone]["master_slop_marked_bed"],
                #args.chipseq["histones"][histone]["master_regions"],# fix this
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
