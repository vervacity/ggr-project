

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
    # ATAC-seq
    # =============================

    #WORK_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/geo_submit/ggr.GEO.atac_submit.2021-07-30"
    WORK_DIR = "/srv/scratch/dskim89/ggr.GEO.atac_submit.2021-07-30"
    
    # SAMPLES
    # use GEOquery in R to get all GSMs (the curators need the old GSMs to associate)
    # also provide metadata table (link ENCODE portal accessions to GEO accessions)
    # inputs: encode_to_gsm_file (metadata table, downloaded from ENCODE portal)
    # outputs: list of GSMs per GSE, complete list of GSMs
    encode_to_gsm_file = "{}/atac.encode_to_geo.tsv".format(WORK_DIR)
    encode_to_gsm_updated_file = "{}.with_GSMs.tsv".format(encode_to_gsm_file.split(".tsv")[0])
    if not os.path.isfile(encode_to_gsm_updated_file):
        get_gsms = "{}/get_gsms_from_encode_metadata.R {} {}".format(
            SCRIPT_DIR, encode_to_gsm_file, encode_to_gsm_updated_file)
        print get_gsms
        os.system(get_gsms)

    # get GSMs as list to add to file
    metadata = pd.read_csv(encode_to_gsm_updated_file, sep="\t")
    gsms = ";".join(metadata["GSM"].values)
    print gsms
        
    # PROCESSED DATA FILES
    # get files into an upload folder, get md5sums, keep organized to easily copy/paste
    # into seq_metadata
    
    # submit: biorep tagAligns <- don't submit alignments
    if False:
        atac_tagaligns_dir = "{}/processed.atac.2016-10-09.bioreps".format(PROCESSED_DIR)
        tagalign_files = sorted(glob.glob("{}/*tagAlign.gz".format(atac_tagaligns_dir)))
        tagalign_table_file = "{}/tagaligns.metadata.tsv".format(WORK_DIR)
        if not os.path.isfile(tagalign_table_file):
            tagalign_table = setup_processed_files(
                tagalign_files, "BED", upload_dir=WORK_DIR, table_name=tagalign_table_file)

    # processed dir
    atac_processed_dir = "{}/processed.atac.2016-10-09.bioreps_peaks_signals".format(PROCESSED_DIR)
    
    # submit: idr peaks
    # NOTE: rep narrowPeaks not submitted since never used in paper, only IDR peaks used
    idr_peaks = sorted(glob.glob("{}/*/peak/idr/true_reps/rep1-rep2/*IDR*filt.narrowPeak.gz".format(atac_processed_dir)))
    idr_table_file = "{}/idr_peaks.metadata.tsv".format(WORK_DIR)
    if not os.path.isfile(idr_table_file):
        idr_table = setup_processed_files(
            idr_peaks, "BED",
            upload_dir=WORK_DIR, table_name=idr_table_file)
    
    # submit: pooled bigwigs
    # NOTE: rep bigwigs not submitted since never used in paper, only pooled bigwigs used
    #bigwigs = sorted(glob.glob("{}/*/signal/macs2/rep*/*bigwig".format(atac_processed_dir)))
    bigwigs = sorted(glob.glob("{}/*/signal/macs2/pooled_rep/*bigwig".format(atac_processed_dir)))
    bigwig_table_file = "{}/bigwigs.metadata.tsv".format(WORK_DIR)
    if not os.path.isfile(bigwig_table_file):
        bigwig_table = setup_processed_files(
            bigwigs, "BIGWIG",
            upload_dir=WORK_DIR, table_name=bigwig_table_file)
    
    # submit: other useful processed files
    GGR_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a"
    processed_files = []
    processed_files += ["{}/data/ggr.atac.idr.master.bed.gz".format(GGR_DIR)] # union peaks
    processed_files += ["{}/data/ggr.atac.idr.summits.bed.gz".format(GGR_DIR)] # union summits
    processed_files += ["{}/data/ggr.atac.ends.counts.mat.txt.gz".format(GGR_DIR)] # count mat
    processed_files += ["{}/data/ggr.atac.ends.counts.pooled.rlog.mat.txt.gz".format(GGR_DIR)] # rlog norm mat

    ggr_traj_dir = "{}/results/atac/timeseries/dp_gp/reproducible/hard/reordered".format(GGR_DIR)
    processed_files += ["{}/ggr.atac.reproducible.hard.reordered.clustering.txt".format(ggr_traj_dir)] # clusters

    processed_table_file = "{}/other_processed.metadata.tsv".format(WORK_DIR)
    if not os.path.isfile(processed_table_file):
        processed_table = setup_processed_files(
            processed_files, "TXT",
            upload_dir=WORK_DIR, table_name=processed_table_file)
    
    # submit: epigenome files?
    epigenome_files = []
    
    # RAW FILES: none for ATAC-seq
    # PAIRED-END: EXPTS none for ATAC-seq

    return

main()
