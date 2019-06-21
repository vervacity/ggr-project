#!/usr/bin/env python

# fig 3b preprocess files for allele sensitive ATAC analysis
# description: code to "unfix" bams (ie in the appropriate situations,
# remove the 1bp position shift to realign base pairs and position)
# this is necessary to correctly pull the base pair from the correct position
# data: /mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2016-10-06.fixed_filtered_align

import os
import sys
import glob
import pysam

import numpy as np

_READ_LEN = 76

def unfix_bam(in_bam_file, out_bam_file):
    """read in bam file and unfix
    """
    bam_file = pysam.AlignmentFile(in_bam_file, 'rb')
    to_stdout = pysam.AlignmentFile(out_bam_file, "wb", template=bam_file)

    seen_read_1 = False
    seen_read_2 = False

    num_fixed = 0
    for read in bam_file:

        # process when read 1 AND read 2 are seen
        # only important at the very beginning of file
        if read.is_read1:
            read1 = read
            seen_read_1 = True
        if read.is_read2:
            read2 = read
            seen_read_2 = True
        if not seen_read_1 or not seen_read_2:
            continue

        # check to see if read pair fits the problem (in reverse)
        if read1.query_name == read2.query_name:
            # Check 4 - only fix reads that came from the specific chip in the flowcell
            check_4 = ('H3MCTBBXX' in read1.query_name)

            # Always change the read 1 start point for appropriate reads UNLESS unmapped or will become unmapped (ie 0 or 1 coordinate)
            if (read1.reference_start > 1) and check_4:
                read1.reference_start += 1
                num_fixed += 1

            to_stdout.write(read1)
            to_stdout.write(read2)

    print "adjusted read total:", num_fixed

    return None


def main():

    # server specific
    WORK_DIR = "/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2016-10-06.fixed_filtered_align"
    #OUT_DIR = '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2019-06-04.bams_bp-position-matched'
    OUT_DIR = "/srv/scratch/dskim89/ggr/fig_3-b"
    os.system("mkdir -p {}".format(OUT_DIR))

    # first unfix BAM files
    bam_files = sorted(glob.glob("{}/*/align/rep1/*bam".format(WORK_DIR)))
    print len(bam_files)
    for bam_file in bam_files:

        # read name sort
        tmp_sort_bam_file = "{}/{}.sort.bam".format(
            OUT_DIR, os.path.basename(bam_file).split(".bam")[0])
        sort_cmd = "samtools sort -o {} -n -@ 16 {}".format(
            tmp_sort_bam_file, bam_file)
        if not os.path.isfile(tmp_sort_bam_file):
            os.system(sort_cmd)
            
        # and then unfix
        unsorted_fixed_bam_file = "{}/{}.unfix.tmp.bam".format(OUT_DIR, os.path.basename(bam_file).split(".bam")[0])
        if not os.path.isfile(unsorted_fixed_bam_file):
            unfix_bam(tmp_sort_bam_file, unsorted_fixed_bam_file)

        # and finally sort that file
        out_bam_file = "{}/{}.unfix.bam".format(OUT_DIR, os.path.basename(bam_file).split(".bam")[0])
        sort_cmd = "samtools sort -o {} -@ 16 {}".format(
            out_bam_file, unsorted_fixed_bam_file)
        if not os.path.isfile(out_bam_file):
            os.system(sort_cmd)

    # and clean up
    os.system("rm {}/*sort.bam".format(OUT_DIR))
    os.system("rm {}/*tmp.bam".format(OUT_DIR))

    # merge by bio rep
    donors = ['D1', 'D2']
    timepoints = np.arange(0, 6.5, 0.5)
    for timepoint in timepoints:
        for donor in donors:
            main_day = str(timepoint).split('.')[0]
            half_day = str(timepoint).split('.')[1] == '5'
            donor_num = donor.split('D')[1]

            if not half_day:
                which_files = ('find {0}/ -regex ".*Day{1}.*nodup.unfix.bam" | '
                               'grep {2} | '
                               'grep -v ".*\({1}-5\|{1}5\).*" | '
                               'grep -v "5and"').format(OUT_DIR, main_day, donor)
                pool_files = ('{0} | '
                              'xargs samtools merge {3}/primary_keratinocyte-d{1}0.GGR.Stanford_Greenleaf.ATAC-seq.b{2}.fixedtrim.PE2SE.nodup.unfix.bam').format(
                                  which_files, main_day, donor_num, OUT_DIR)

            else:
                which_files = ('find {0}/ -regex ".*Day{1}.*nodup.unfix.bam" | '
                               'grep {2} | '
                               'grep ".*\({1}-5\|{1}5\).*" | '
                               'grep -v "5and"').format(OUT_DIR, main_day, donor)
                pool_files = ('{0} | '
                              'xargs samtools merge {3}/primary_keratinocyte-d{1}5.GGR.Stanford_Greenleaf.ATAC-seq.b{2}.fixedtrim.PE2SE.nodup.unfix.bam').format(
                                  which_files, main_day, donor_num, OUT_DIR)

            print which_files
            os.system(which_files)
            print pool_files
            os.system(pool_files)
            print "DONE"

    # and do again as pooled (NOTE: do not use pooled, confounds genotyping)
    for timepoint in timepoints:
        main_day = str(timepoint).split('.')[0]
        half_day = str(timepoint).split('.')[1] == '5'

        if not half_day:
            which_files = ('find {0}/ -regex ".*Day{1}.*nodup.unfix.bam" | '
                           #'grep {2} | '
                           'grep -v ".*\({1}-5\|{1}5\).*" | '
                           'grep -v "5and"').format(OUT_DIR, main_day)
            pool_files = ('{0} | '
                          'xargs samtools merge {2}/primary_keratinocyte-d{1}0.GGR.Stanford_Greenleaf.ATAC-seq.pooled.fixedtrim.PE2SE.nodup.unfix.bam').format(
                              which_files, main_day, OUT_DIR)

        else:
            which_files = ('find {0}/ -regex ".*Day{1}.*nodup.unfix.bam" | '
                           #'grep {2} | '
                           'grep ".*\({1}-5\|{1}5\).*" | '
                           'grep -v "5and"').format(OUT_DIR, main_day)
            pool_files = ('{0} | '
                          'xargs samtools merge {2}/primary_keratinocyte-d{1}5.GGR.Stanford_Greenleaf.ATAC-seq.pooled.fixedtrim.PE2SE.nodup.unfix.bam').format(
                              which_files, main_day, OUT_DIR)

        if half_day:
            print which_files
            os.system(which_files)
            print pool_files
            os.system(pool_files)
        print "DONE"

    return None

main()
