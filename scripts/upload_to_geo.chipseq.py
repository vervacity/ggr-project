

import os
import glob
import subprocess

import numpy as np
import pandas as pd


def get_md5sum(filename):
    """get md5sum
    """
    md5sum_val = subprocess.check_output(["md5sum", filename])
    md5sum_val = md5sum_val.split()[0]
    
    return md5sum_val


def setup_processed_files(
        processed_file_names,
        file_type,
        upload_dir=None,
        table_name=None):
    """set up table of processed files to easily put into seq_metadata sheet
    """
    file_names = []
    file_types = []
    file_checksums = []

    for file_name in processed_file_names:
        print file_name
        # copy file over for FTP upload dir
        if upload_dir is not None:
            copy_cmd = "rsync -avz --progress {} {}/".format(file_name, upload_dir)
            print copy_cmd
            os.system(copy_cmd)
        # add to metadata
        file_names.append(os.path.basename(file_name))
        file_types.append(file_type)
        file_checksums.append(get_md5sum(file_name))

    processed_table = pd.DataFrame({
        "file_name": file_names,
        "file_types": file_types,
        "file_checksums": file_checksums})
    processed_table = processed_table.sort_values("file_name")
    if table_name is not None:
        processed_table.to_csv(
            table_name, sep="\t", header=True, index=False)
    
    return processed_table


def main():
    """set up data uploads to GEO, centered around
    making it easy to copy/paste into GEO seq_template
    """
    SCRIPT_DIR = "/users/dskim89/git/ggr-project/scripts"
    PROCESSED_DIR = "/mnt/lab_data/kundaje/projects/skin/data/bds"
    
    # =============================
    # ChIP-seq
    # =============================

    #WORK_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/geo_submit/ggr.GEO.chipseq_submit.2021-07-30"
    WORK_DIR = "/srv/scratch/dskim89/ggr.GEO.chipseq_submit.2021-07-30"
    
    # SAMPLES
    # use GEOquery in R to get all GSMs (the curators need the old GSMs to associate)
    # also provide metadata table (link ENCODE portal accessions to GEO accessions)
    # inputs: encode_to_gsm_file (metadata table, downloaded from ENCODE portal)
    # outputs: list of GSMs per GSE, complete list of GSMs
    encode_to_gsm_file = "{}/chipseq.encode_to_geo.tsv".format(WORK_DIR)
    encode_to_gsm_updated_file = "{}.with_GSMs.tsv".format(encode_to_gsm_file.split(".tsv")[0])
    #if not os.path.isfile(encode_to_gsm_updated_file):
    if True:
        get_gsms = "{}/get_gsms_from_encode_metadata.R {} {}".format(
            SCRIPT_DIR, encode_to_gsm_file, encode_to_gsm_updated_file)
        print get_gsms
        os.system(get_gsms)

    quit()
        
    # PROCESSED DATA FILES
    # get files into an upload folder, get md5sums, keep organized to easily copy/paste
    # into seq_metadata
    chipseq_processed_dir = "{}/processed.chipseq.2017-01-23.histones".format(PROCESSED_DIR)
    
    # submit: overlap peaks
    # NOTE: rep narrowPeaks not submitted since never used in paper, only overlap peaks used
    overlap_peaks = sorted(glob.glob("{}/*/peak/macs2/overlap/*overlap*filt.narrowPeak.gz".format(chipseq_processed_dir)))
    overlap_table_file = "{}/overlap_peaks.metadata.tsv".format(WORK_DIR)
    if not os.path.isfile(overlap_table_file):
        overlap_table = setup_processed_files(
            overlap_peaks, "BED",
            upload_dir=WORK_DIR, table_name=overlap_table_file)
    
    # submit: bigwigs
    # NOTE: rep bigwigs not submitted since never used in paper, only pooled bigwigs used
    #bigwigs = sorted(glob.glob("{}/*/signal/macs2/rep*/*bigwig".format(chipseq_processed_dir)))
    bigwigs = sorted(glob.glob("{}/*/signal/macs2/pooled_rep/*bw".format(chipseq_processed_dir)))
    bigwigs = [bigwig for bigwig in bigwigs if "CTCF" not in bigwig]
    bigwig_table_file = "{}/bigwigs.metadata.tsv".format(WORK_DIR)
    if not os.path.isfile(bigwig_table_file):
        bigwig_table = setup_processed_files(
            bigwigs, "BIGWIG",
            upload_dir=WORK_DIR, table_name=bigwig_table_file)
    
    # submit: other useful processed files
    GGR_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a"
    processed_files = []
    histones = ["H3K27ac", "H3K4me1", "H3K27me3"]
    for histone in histones:
        processed_files += ["{}/data/ggr.histone.{}.overlap.master.bed.gz".format(GGR_DIR, histone)] # union peaks
        processed_files += ["{}/data/ggr.histone.{}.overlap.midpoints.counts.mat.txt.gz".format(GGR_DIR, histone)] # count mat
        processed_files += ["{}/data/ggr.histone.{}.overlap.midpoints.counts.pooled.rlog.mat.txt.gz".format(GGR_DIR, histone)] # rlog norm mat
        processed_files += ["{}/data/ggr.histone.{}.midpoints.counts.mat.txt.gz".format(GGR_DIR, histone)] # count mat, ATAC-centric
        processed_files += ["{}/data/ggr.histone.{}.midpoints.counts.pooled.rlog.mat.txt.gz".format(GGR_DIR, histone)] # rlog norm, ATAC-centric      
    processed_table_file = "{}/other_processed.metadata.tsv".format(WORK_DIR)
    processed_table = setup_processed_files(
        processed_files, "TXT",
        upload_dir=WORK_DIR, table_name=processed_table_file)

        
    # RAW FILES: none for ChIP-seq
    # PAIRED-END: EXPTS none for ChIP-seq

    return

main()
