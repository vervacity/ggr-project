#!/usr/bin/env python

# fig 3b processing
# description: set up reads for allele-specific ATAC (quasar) calls
# data: /srv/scratch/dskim89/ggr/ggr.atac.2019-06-06.ase

import os
import glob

import pandas as pd
import numpy as np


def filter_bam_file(in_bam, out_bam):
    """ filter bam file
    """
    # make header
    header_bam = "{}.header.bam".format(out_bam.split(".bam")[0])
    make_header = "samtools view -H {} > {}".format(
        in_bam, header_bam)
    os.system(make_header)

    # filter
    cmd = (
        "samtools view -q 40 {} | " # q40 filter
        "awk '$6 !~ /D|I/' | " # check CIGAR, remove indels
        "awk '$12 ~ /^MD:Z:([0-9]*[A-Z][0-9]*){{0,1}}[0-9]*$/' | " # check MD, only allow 1 mismatch through
        "awk '$12 !~ /^MD:Z:[0-5][A-Z][0-9]*/' | " # just throw out reads with mismatch within 5bp of start
        "cat {} - | "
        "samtools view -bS - > {}").format(
            in_bam, header_bam, out_bam)
    print cmd
    os.system(cmd)

    # index the file
    index_cmd = "samtools index {}".format(out_bam)
    print index_cmd
    os.system(index_cmd)
    
    # clean up
    os.system("rm {}".format(header_bam))
    
    return None


def build_quasar_input_file(
        snp_file,
        bam_file,
        out_file,
        fasta_file,
        quasar_preprocess_script):
    """build quasar input file
    """
    # pileup
    pileup_file = "{}.1KG_filt.pileup.gz".format(
        bam_file.split(".bam")[0])
    pileup_cmd = "samtools mpileup -f {} -l {} {} | gzip -c > {}".format(
        fasta_file, snp_file, bam_file, pileup_file)
    if not os.path.isfile(pileup_file):
        print pileup_cmd
        os.system(pileup_cmd)
    
    # convert to BED with filtering
    pileup_bed_file = "{}.1KG_filt.pileup.bed.gz".format(
        bam_file.split(".bam")[0])
    to_bed_cmd = (
        "zcat {} | "
        "awk -v OFS='\t' "
        "'{{ if ($4>0 && $5 !~ /[^\^][<>]/ && "
        "$5 !~ /\+[0-9]+[ACGTNacgtn]+/ && "
        "$5 !~ /-[0-9]+[ACGTNacgtn]+/ "
        "&& $5 !~ /[^\^]\*/) "
        "print $1,$2-1,$2,$3,$4,$5,$6 }}' | "
        "sortBed -i stdin | "
        "intersectBed -a stdin -b {} -wo | "
        "cut -f 1-7,11-14 | "
        "gzip > {}").format(
            pileup_file, snp_file, pileup_bed_file)
    if not os.path.isfile(pileup_bed_file):
        print to_bed_cmd
        os.system(to_bed_cmd)

    # and then process to quasar input file
    quasar_preprocess_cmd = "R --vanilla --args {} < {}".format(
        pileup_bed_file,
        quasar_preprocess_script)
    if not os.path.isfile(out_file):
        print quasar_preprocess_cmd
        os.system(quasar_preprocess_cmd)

    assert os.path.isfile(out_file)
    
    return None


def main():
    """run ASE on GGR ATAC
    """
    # server specific inputs
    WORK_DIR = "."
    FASTA = "/mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa"
    BAM_DIR = "/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2019-06-04.bams_bp-position-matched"
    SNP_FILE_URL = "http://genome.grid.wayne.edu/centisnps/files/1KG_SNPs_filt.bed.gz"
    quasar_preprocess_script = "/users/dskim89/git/ggr-project/figs/fig_3.motifs_and_tfs/fig_3-b.1b.convertPileupToQuasar.R"
    quasar_script = "/users/dskim89/git/ggr-project/figs/fig_3.motifs_and_tfs/fig_3-b.1c.run_quasar.R"
    
    # pull bam files
    bam_files = sorted(glob.glob("{}/*.b*bam".format(BAM_DIR)))

    # get the snp file
    snp_file = "1KG_SNPs_filt.bed.gz"
    if not os.path.isfile(snp_file):
        download_snp_file = "wget {}".format(SNP_FILE_URL)
        os.system(download_snp_file)
    
    # first filter BAMs
    # up the filter quality
    # adjust positions/base pairs of BAM file (trimming)
    # remove reads with indels, only keep reads that match perfectly or have 1 mismatch
    filt_bams = []
    for bam_file in bam_files:
        out_bam = "{}/{}.ase_filt.bam".format(
            WORK_DIR, os.path.basename(bam_file).split(".bam")[0])
        if not os.path.isfile(out_bam):
            filter_bam_file(bam_file, out_bam)
        filt_bams.append(out_bam)
        
    # then set up files for quasar
    quasar_inputs = []
    for bam_file in filt_bams:
        quasar_input_file = "{}.1KG_filt.quasar.in.gz".format(
            bam_file.split(".bam")[0])
        if not os.path.isfile(quasar_input_file):
            build_quasar_input_file(snp_file, bam_file, quasar_input_file, FASTA, quasar_preprocess_script)
        quasar_inputs.append(quasar_input_file)

    # run each individual through quasar, generate individual results and pooled results
    rep1_results_file = "b1.results_ALL.txt"
    if not os.path.isfile(rep1_results_file):
        rep1_inputs = [filename for filename in quasar_inputs if "b1" in filename]
        run_quasar = "{} b1 {}".format(
            quasar_script,
            " ".join(rep1_inputs))
        print run_quasar
        os.system(run_quasar)
    
    rep2_results_file = "b2.results_ALL.txt"
    if not os.path.isfile(rep2_results_file):
        rep2_inputs = [filename for filename in quasar_inputs if "b2" in filename]
        run_quasar = "Rscript run_quasar.R b2 {}".format(
            " ".join(rep2_inputs))
        print run_quasar
        os.system(run_quasar)

    # merge the two rep files for a final file with sig variants
    # make a VCF file, name column has sig?
    rep1_data = pd.read_csv(rep1_results_file, sep="\t", index_col=0)
    rep2_data = pd.read_csv(rep2_results_file, sep="\t", index_col=0)

    all_data = rep1_data.merge(rep2_data, how="outer", left_index=True, right_index=True)
    all_data = all_data.fillna(0)
    all_data = all_data[sorted(all_data.columns)]
    id_fields = all_data.index.to_series().str.split("_", n=3, expand=True)
    #all_data["rsid"] = id_fields[0]
    #print all_data
    
    # and add column that considers whether sig across any timepoint?
    # need to keep slope and sig effect and max time point of effect?
    fdr = 0.10
    qval_columns = [colname for colname in all_data.columns if "qval" in colname]
    qvals = all_data[qval_columns]
    sig = np.any((qvals < fdr) & (qvals != 0), axis=-1).astype(int)
    all_data["sig"] = sig

    beta_columns = [colname for colname in all_data.columns if "beta" in colname]
    betas = all_data[beta_columns]
    all_data["beta_max_idx"] = np.argmax(np.abs(betas.values), axis=-1) * 2
    betas_max = []
    for i in range(all_data.shape[0]):
        betas_max.append(all_data.values[i,all_data["beta_max_idx"].values[i]])
    all_data["beta_max"] = betas_max
    #print all_data

    # pull in ref/alt from snp file
    all_data["rsid"] = id_fields[0].values
    snps = pd.read_csv(
        snp_file, sep="\t",
        header=None,
        names=["chr", "start", "stop", "rsid", "ref", "alt", "maf"])
    all_data = all_data.merge(snps, how="left", on="rsid")
    #print all_data
    
    # save this out, reduce data
    all_data_file = "all.results_ALL.txt"
    all_data.to_csv(all_data_file, sep="\t")

    all_reduced_file = "all.results.slim.txt"
    all_reduced = all_data[["sig", "beta_max"]]
    all_reduced.to_csv(all_reduced_file, sep="\t")

    # convert to VCF format
    # keep sig and beta val
    vcf_file = "all.1KG_filt.vcf"
    vcf_data = all_data[
        ["chr", "start", "rsid", "ref", "alt"]]
    vcf_data["chr"] = vcf_data["chr"].str.replace("chr", "").values
    vcf_data.columns = ["#CHROM", "POS", "ID", "REF", "ALT"]
    vcf_data["QUAL"] = 1
    vcf_data["FILTER"] = 1
    vcf_data["INFO"] = "sig=" + all_data["sig"].map(str) + ";beta=" + all_data["beta_max"].map(str)
    vcf_data.to_csv(vcf_file, sep="\t", index=False)

    print np.sum(all_data["sig"].values)
    
    return

main()
