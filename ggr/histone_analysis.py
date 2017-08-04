"""Contains code to run histone analyses
"""


import os
import glob
import logging
import signal

import pandas as pd

from ggr.util.utils import add_folder
from ggr.util.bed_utils import merge_regions
from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel



def make_midpoints_file(bedpe_file, out_file):
    """Get the midpoint base pair from a paired ended dataset.
    This is the cleanest signal for where a factor was bound 
    to the fragment.
    """
    assert ".bedpe.gz" in bedpe_file

    # for multiprocessing purposes
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    
    # awk command
    # for BEDPE, if first strand is positive, grab start of 1 and end of 2
    # if first strand is negative, grab end of 1 and start of 2
    make_midpoints = (
        "zcat {0} | "
        "awk -F '\t' 'BEGIN {{ OFS=\"\t\" }} "
        "{{ if ($9==\"+\") {{ $2=$2+int(($6-$2)/2); $3=$2+int(($6-$2)/2)+1 }} "
        "else if ($9==\"-\") {{ $2=$5+int(($3-$5)/2); $3=$5+int(($3-$5)/2)+1 }} "
        "print $1\"\t\"$2\"\t\"$3\"\tN\t1000\t\"$9 }}' | "
        "sort -k1,1 -k2,2n | "
        "gzip -c "
        "> {1}").format(
            bedpe_file,
            out_file)

    print make_midpoints
    os.system(make_midpoints)
    
    return None

def get_counts_from_tagalign(
        master_regions_file,
        tagalign_file,
        out_counts_file):
    """Get counts in regions from tagAlign of reads
    """
    get_coverage = (
        "bedtools coverage "
        "-sorted "
        "-counts "
        "-a {0} -b {1} | "
        "gzip -c > {2}").format(
            master_regions_file,
            tagalign_file,
            out_counts_file)
    print get_coverage
    os.system(get_coverage)
    
    return


def get_counts_in_regions_from_dir(
        reads_dir,
        master_regions_file,
        out_dir,
        unique_string="",
        count_type="cuts"):
    """given a directory with tagAlign files, get overlap
    with master regions and save out to output folder
    Use unique string if multiple assays in same folder and 
    you only want one of them
    """
    count_queue = setup_multiprocessing_queue()

    reads_files = sorted(
        glob.glob("{0}/*{1}*{2}.tagAlign.gz".format(
            reads_dir, unique_string, count_type)))
    
    for reads_file in reads_files:
        count_file = "{0}/{1}.counts.txt.gz".format(
            out_dir,
            os.path.basename(reads_file).split(
                ".{}".format(count_type))[0])
        if not os.path.isfile(count_file):
            print count_file
            count_queue.put([
                get_counts_from_tagalign,
                [master_regions_file, reads_file, count_file]])

    run_in_parallel(count_queue)

    return


def make_counts_mat_from_dir(counts_dir, out_mat_file):
    """given a directory with counts files (from bedtools coverage)
    merge to make a matrix file
    """
    counts_files = sorted(glob.glob('{}/*counts.txt.gz'.format(counts_dir)))
    
    coverage_df = pd.DataFrame()
    if not os.path.isfile(out_mat_file):
        for counts_file in counts_files:
            filename_fields = os.path.basename(counts_file).split('.')

            # some hackyness around names
            header = '{0}_{1}'.format(
                filename_fields[0].split('-')[1],
                filename_fields[4])

            # and now the counts
            counts = pd.read_csv(counts_file, sep='\t', header=None)
            counts.columns = ['chrom', 'start', 'stop', 'counts']

            coverage_df[header] = counts['counts']
        coverage_df = coverage_df.set_index(
            counts['chrom'] + ':' + counts['start'].map(str) + '-' + counts['stop'].map(str))
        coverage_df.to_csv(out_mat_file, sep='\t', index_label=False, compression='gzip')
    
    return


def make_chipseq_count_matrix(
        master_regions_file,
        bedpe_file_list,
        out_mat_file,
        mark_name,
        tmp_dir="."):
    """For regions in master file list, return a matrix file 
    with the counts of midpoints of fragments in these regions
    """
    # generate (sorted) midpoint files from bedpe files
    midpoints_dir = "{}/midpoints".format(tmp_dir)
    os.system("mkdir -p {}".format(midpoints_dir))
    
    midpoints_queue = setup_multiprocessing_queue()
    for bedpe_file in bedpe_file_list:
        midpoints_file = "{0}/{1}.midpoints.tagAlign.gz".format(
            midpoints_dir,
            os.path.basename(bedpe_file).split(".bedpe")[0])
        if not os.path.isfile(midpoints_file):
            print midpoints_file
            midpoints_queue.put([
                make_midpoints_file,
                [bedpe_file, midpoints_file]])
            
    run_in_parallel(midpoints_queue)

    # overlap with master regions and then run bedtools coverage
    counts_dir = "{}/counts".format(tmp_dir)
    os.system("mkdir -p {}".format(counts_dir))
    get_counts_in_regions_from_dir(
        midpoints_dir,
        master_regions_file,
        counts_dir,
        unique_string=mark_name,
        count_type="midpoints")

    # merge all into matrix
    make_counts_mat_from_dir(counts_dir, out_mat_file)
    
    return



def run(args):
    """Pipeline for histone analyses
    """

    # for each histone:
    histones = args.chipseq["histones"].keys()

    for histone in histones:

        logging.info("working on {}".format(histone))
        histone_prefix = "ggr.{}".format(histone)
        histone_dir = args.folders["{}_dir".format(histone)]
        
        # first do an intersect of the raw BED files,
        # the current overlap (1/23/2017) was unstable
        days = ["d0", "d3", "d6"]
        intersect_dir = "{}/intersect".format(histone_dir)
        os.system("mkdir -p {}".format(intersect_dir))
        histone_narrowPeak_files = glob.glob("{0}/{1}".format(
            args.chipseq["data_dir"],
            args.chipseq["histones"][histone]["narrowPeak_glob"]))
        
        for day in days:
            # get right files
            day_histone_files = [peakfile
                                 for peakfile in histone_narrowPeak_files
                                 if day in peakfile]

            assert len(day_histone_files) == 2
            
            # intersect
            new_file = "{0}/{1}.replicated.narrowPeak.gz".format(
                intersect_dir,
                os.path.basename(day_histone_files[0]).split(".narrowPeak")[0])
            if not os.path.isfile(new_file):
                bed_intersect = ("bedtools intersect -a {0} -b {1} | "
                                 "gzip -c > {2}").format(
                                     day_histone_files[0],
                                     day_histone_files[1],
                                     new_file)
                logging.info(bed_intersect)
                os.system(bed_intersect)
        
        # generate master regions file
        args.chipseq["histones"][histone]["master_regions"] = "{0}/ggr.{1}.macs2.master.bed.gz".format(
            args.folders["data_dir"], histone)
        if not os.path.isfile(args.chipseq["histones"][histone]["master_regions"]):
            histone_replicated_files = glob.glob("{0}/*{1}*.replicated.narrowPeak.gz".format(
                intersect_dir, histone))
            logging.info("Master regions using: {}".format(" ".join(histone_replicated_files)))
            merge_regions(histone_replicated_files, args.chipseq["histones"][histone]["master_regions"])
            
        # get midpoints of fragments in these regions
        bedpe_file_list = sorted(glob.glob("{}/{}".format(
            args.chipseq["data_dir"],
            args.chipseq["histones"][histone]["bedpe_glob"])))
        assert len(bedpe_file_list) == 6

        # make count matrix for histone mark
        args.chipseq["histones"][histone]["counts"] = "{0}/ggr.{1}.midpoints.counts.mat.txt.gz".format(
            args.folders["data_dir"], histone)
        make_chipseq_count_matrix(
            args.chipseq["histones"][histone]["master_regions"],
            bedpe_file_list,
            args.chipseq["histones"][histone]["counts"],
            histone,
            tmp_dir=histone_dir)

        # run DESeq on forward sequential pairs of timepoints
        histone_diff_dir = args.folders["{}_diff_dir".format(histone)]
        args.chipseq["histones"][histone]["dynamic_ids"] = "{0}/ggr.{1}.dynamic.ids.txt.gz".format(
            histone_diff_dir, histone)
        if not os.path.isfile(args.chipseq["histones"][histone]["dynamic_ids"]):
            run_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3} sequential".format(
                args.chipseq["histones"][histone]["counts"],
                "{0}/{1}".format(histone_diff_dir, histone_prefix),
                args.params["timeseries_fdr"],
                args.chipseq["histones"][histone]["dynamic_ids"])
            print run_deseq2
            os.system(run_deseq2)

            
        # enumerate trajectories (9 groups) and keep trajectory numbers

        
        # from here, move onto integration with ATAC info
    

    



    

    return None
