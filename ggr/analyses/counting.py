"""Tools to quickly count reads in a region and analyze
"""

import os
import glob
import signal
import logging

import pandas as pd

from ggr.util.utils import run_shell_cmd

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel


def make_ends_file(tagalign_file, out_file, cluster="sherlock"):
    """Take a tagAlign file and get read ends only
    Generally used for ATAC-seq/DNase-seq, where the key
    information in the read is the END of the read (where 
    the cut in the open site occurred)

    Args:
      tagalign_file: BED file of reads
      out_file: BED file of "reads" (ie, 1bp read ends)
    """
    assert os.path.splitext(out_file)[1] == ".gz"
    assert os.path.splitext(tagalign_file)[1] == ".gz"
    
    # Hack for python multiprocessing
    #if not cluster == "sherlock":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    
    # awk is faster than python
    make_ends = (
        "zcat {0} | "
        "awk -F '\t' 'BEGIN {{ OFS=\"\t\" }} "
        "{{ if ($6==\"+\") {{$3=$2+1}} "
        "else if ($6==\"-\") {{$2=$3-1}} "
        "print $0 }}' | "
        "sort -k1,1 -k2,2n | "
        "gzip -c "
        "> {1}").format(
            tagalign_file, out_file)
    
    print make_ends
    logging.info(make_ends)
    run_shell_cmd(make_ends)

    return None


def make_midpoints_file(bedpe_file, out_file, cluster="sherlock"):
    """Get the midpoint base pair from a paired ended dataset.
    This is the cleanest signal for where a factor was bound 
    to the fragment, ie use for paired ended ChIP-seq data.

    Args:
      bedpe_file: BEDPE format file with 2 paired reads
      out_file: BED file of "reads" (ie, 1bp midpoints)
    """
    assert ".bedpe.gz" in bedpe_file
    assert os.path.splitext(bedpe_file)[1] == ".gz"
    assert os.path.splitext(out_file)[1] == ".gz"

    # for multiprocessing purposes
    if not cluster == "sherlock":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    
    # awk is faster than python
    make_midpoints = (
        "zcat {0} | "
        "awk -F '\t' 'BEGIN {{ OFS=\"\t\" }} "
        "{{ if ($9==\"+\") {{ midpoint=$2+int(($6-$2)/2); $2=midpoint; $3=midpoint+1 }} "
        "else if ($9==\"-\") {{ midpoint=$5+int(($3-$5)/2); $2=midpoint; $3=midpoint+1 }} "
        "print $1\"\t\"$2\"\t\"$3\"\tN\t1000\t\"$9 }}' | "
        "sort -k1,1 -k2,2n | "
        "gzip -c "
        "> {1}").format(
            bedpe_file,
            out_file)

    print make_midpoints
    run_shell_cmd(make_midpoints)
    
    return None


def adjust_reads_from_dir(reads_files, out_dir, adjustment="ends"):
    """Use to do reads adjustments (get ends or midpoints)
    by giving an input and output directory

    Args:
      reads_files: list of BED of reads filenames
      out_dir: where to save them
    """
    # set up queue
    adjust_queue = setup_multiprocessing_queue()

    # go through reads files
    for reads_file in reads_files:

        # check adjustment type and set up accordingly
        if adjustment == "ends":
            out_file = "{0}/{1}.{2}.tagAlign.gz".format(
                out_dir,
                os.path.basename(reads_file).split('.tagAlign')[0],
                adjustment)
            adjust_fn = make_ends_file
            
        elif adjustment == "midpoints":
            out_file = "{0}/{1}.{2}.tagAlign.gz".format(
                out_dir,
                os.path.basename(reads_file).split('.bedpe')[0],
                adjustment)
            adjust_fn = make_midpoints_file

        else:
            "Need proper adjustment type"
            return

        # add to queue
        if not os.path.isfile(out_file):
            print out_file
            adjust_queue.put([
                adjust_fn,
                [reads_file, out_file]])

    # run queue
    run_in_parallel(adjust_queue)
    
    return None


def get_counts_from_tagalign(
        master_regions_file,
        tagalign_file,
        out_counts_file):
    """Get counts in regions from tagAlign of reads

    Args:
      master_regions_file: master BED file
      tagalign_file: BED file of reads
      out_counts_file: output file with counts/region
    """
    assert os.path.splitext(out_counts_file)[1] == ".gz"

    # run bedtools coverage
    # TODO allow read extension?
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
    run_shell_cmd(get_coverage)
    
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

    Args:
      reads_dir: directory with BED files of reads
      master_regions_file: BED file to get counts in
      out_dir: where to put counts files
      unique_string: for glob
      count_type: for naming files
    """
    # setup parallelize queue
    count_queue = setup_multiprocessing_queue()

    # get reads files
    reads_files = sorted(
        glob.glob("{0}/*{1}*{2}.tagAlign.gz".format(
            reads_dir, unique_string, count_type)))

    # run reads_files
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

    return None


def make_counts_mat_from_dir(counts_dir, out_mat_file):
    """given a directory with counts files (from bedtools coverage)
    merge to make a matrix file
    
    Args:
      counts_dir: directory of counts files
      out_mat_file: output matrix file of counts
    """
    # get counts files
    counts_files = sorted(
        glob.glob('{}/*counts.txt.gz'.format(
            counts_dir)))

    # join all the counts
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

        # write out
        coverage_df.to_csv(out_mat_file, sep='\t', index_label=False, compression='gzip')
    
    return


def make_count_matrix(
        master_regions_file,
        bedpe_file_list,
        out_mat_file,
        assay_name,
        adjustment="midpoints",
        tmp_dir="."):
    """Given ChIP-seq (paired ended) or ATAC, get count matrix
    Note that if ChIP-seq, it MUST be paired ended, with BEDPE
    files. If ATAC, must have tagAlign files.

    Args:
      master_regions_file: regions you want to count in
      bedpe_file_list: list of BEDPE formatted reads
      out_mat_file: output file
      mark_name: ChIP-seq mark name in file
      tmp_dir: scratch dir as needed
    """

    # first adjust the reads files to get fragment midpoints
    adjust_dir = "{0}/{1}".format(tmp_dir, adjustment)
    os.system("mkdir -p {}".format(adjust_dir))
    adjust_reads_from_dir(
        bedpe_file_list,
        adjust_dir,
        adjustment=adjustment)

    # then overlap the fragment midpoints with the master regions
    counts_dir = "{}/counts".format(tmp_dir)
    os.system("mkdir -p {}".format(counts_dir))
    get_counts_in_regions_from_dir(
        adjust_dir,
        master_regions_file,
        counts_dir,
        unique_string=assay_name,
        count_type=adjustment)

    # finally merge all into matrix
    make_counts_mat_from_dir(counts_dir, out_mat_file)

    return None


def pull_single_replicate(mat_file, out_file, rep="b1"):
    """From a matrix, pull one replicate out
    """
    assert ".mat" in mat_file
    data = pd.read_table(mat_file, index_col=0)
    
    rep_columns = [
        colname for colname in data.columns
        if rep in colname]
    data_rep = data[rep_columns]
    data_rep.to_csv(out_file, sep='\t', compression="gzip")

    return None


def pool_replicates(mat_file, out_file):
    """Pool replicates, assumes naming
    structure is *_b1, *_b2.
    """
    assert ".mat" in mat_file
    data = pd.read_table(mat_file, index_col=0)

    samples = sorted(
        list(set([colname.split("_")[0]
                  for colname in data.columns])))
    for sample in samples:
        data[sample] = data["{}_b1".format(sample)] + data["{}_b2".format(sample)]
    data_pooled = data[samples].astype(int)
    data_pooled.to_csv(out_file, sep='\t', compression="gzip")

    return None


def split_count_matrix_by_replicate(
        count_matrix_file,
        rep1_file,
        rep2_file,
        pooled_file):
    """Takes an input matrix and splits into replicates 
    and a pooled file.
    """
    assert ".mat" in count_matrix_file
    
    # rep1
    pull_single_replicate(count_matrix_file, rep1_file, rep="b1")

    # rep2
    pull_single_replicate(count_matrix_file, rep2_file, rep="b2")

    # pooled
    pool_replicates(count_matrix_file, pooled_file)

    return None

