# Description: RNA-specific analyses


import os
import glob
import gzip
import pandas as pd
import numpy as np

from ggr.util.utils import run_shell_cmd


def make_rsem_matrix(quant_files, out_file, colname="expected_count"):
    """Given a list of RSEM files, extract the desired column and make matrix
    column for counts: expected_counts
    column for tpm: TPM
    """
    signal_mat = pd.DataFrame()
    for quant_file in quant_files:
        fields = quant_file.split("/")[-2].split(".")
        prefix = '{0}_{1}'.format(fields[0].split('-')[1], fields[4])
        
        data = pd.read_table(quant_file, sep="\t")
        signal_mat[prefix] = data[colname]
        
    signal_mat["gene_id"], _ = data["gene_id"].str.split(".", 1).str
    signal_mat = signal_mat.set_index("gene_id")
    if colname == "expected_count":
        signal_mat = signal_mat.astype(int)
    signal_mat.to_csv(out_file, compression="gzip", sep='\t')
    
    return None


def threshold_empirically_off_genes(mat_file, out_file, thresh=1.0):
    """Given a matrix file, get a distribution and determine empirical threshold
    """
    data = pd.read_table(mat_file, index_col=0)
    data_thresh = data.loc[~(data < thresh).all(axis=1)]

    data_thresh.to_csv(out_file, compression="gzip", sep='\t')
    
    return None


def filter_clusters_for_tfs(
        cluster_file,
        out_cluster_file,
        mapping_file,
        keep_list):
    """Filter cluster file for id strings in list
    """
    keep_list = [str(keep) for keep in keep_list]

    # read in mapping file
    mappings = {}
    with gzip.open(mapping_file, "r") as fp:
        for line in fp:
            fields = line.strip().split("\t")
            mappings[fields[0]] = fields[1]
    
    # open files and write out
    with open(cluster_file, "r") as fp:
        with open(out_cluster_file, "w") as out:

            for line in fp:
                fields = line.strip().split("\t")
                if "cluster" in fields[0]:
                    out.write(line)
                    continue

                current_id_string = fields[1]
                try:
                    mapped_symbol = mappings[current_id_string]
                except:
                    continue
                for keep_id in keep_list:
                    if keep_id in mapped_symbol:
                        out.write("{}\t{}\n".format(line.strip(), mapped_symbol))
                        
    return None


def filter_clusters_for_tss(
        cluster_file,
        out_prefix,
        tss_bed_file):
    """given cluster file, build TSS files for each cluster
    """
    # read in TSS file to dict
    geneid_to_tss = {}
    with gzip.open(tss_bed_file, "r") as fp:
        for line in fp:
            fields = line.strip().split()
            geneid_to_tss[fields[3]] = line
    
    # open files and write out
    with open(cluster_file, "r") as fp:
        for line in fp:
            if line.startswith("cluster"):
                continue
            
            fields = line.strip().split()
            cluster_file = "{}.cluster-{}.tmp.bed.gz".format(out_prefix, fields[0])
            try:
                with gzip.open(cluster_file, "a") as out:
                    out.write(geneid_to_tss[fields[1]])
            except:
                print "did not find {}".format(fields[1])

    # sort 
    cluster_tmp_files = glob.glob("{}*tmp.bed.gz".format(out_prefix))
    for cluster_tmp_file in cluster_tmp_files:
        cluster_file = "{}.bed.gz".format(cluster_tmp_file.split(".tmp")[0])
        sort_file = "zcat {} | sort -k1,1 -k2,2n | gzip -c > {}".format(
            cluster_tmp_file,
            cluster_file)
        run_shell_cmd(sort_file)

    # delete tmp files
    run_shell_cmd("rm {}*tmp*".format(out_prefix))
    
    return None


def get_neighborhood(
        bed_file,
        other_bed_file,
        out_bed_file,
        extend_len,
        chromsizes):
    """given a bed file and extend distance, call back regions
    from other bed file that are in that neighborhood
    """
    neighborhood = (
        "bedtools slop -i {0} -g {1} -b {2} | "
        "bedtools intersect -u -a {3} -b stdin | "
        "gzip -c > {4}").format(
            bed_file,
            chromsizes,
            extend_len,
            other_bed_file,
            out_bed_file)
    run_shell_cmd(neighborhood)

    return None


def get_nearest_regions(
        bed_file,
        other_bed_file,
        out_bed_file,
        opt_string=""):
    """use bedtools closest to get nearest region
    """
    closest = (
        "bedtools closest {0} -a {1} -b {2} | "
        "awk -F '\t' '{{ print $7\"\t\"$8\"\t\"$9\"\t\"$4\"\t\"$10\"\t\"$6 }}' | "
        "gzip -c > {3}").format(
            opt_string,
            bed_file,
            other_bed_file,
            out_bed_file)
    run_shell_cmd(closest)

    return None


def get_highly_expressed_genes(
        rna_mat_file,
        gene_list_file,
        out_file,
        percentile=90):
    """select genes that are most highly expressed and return gene list
    """
    # read in
    data = pd.read_table(rna_mat_file, index_col=0)
    gene_list = pd.read_table(gene_list_file, index_col=0)
    
    # filter
    data = data[data.index.isin(gene_list["id"])]

    # and then reduce
    data["max"] = data.max(axis=1)
    
    # and select top
    k = int(((100 - percentile) / 100.) * data.shape[0])
    data = data.nlargest(k, "max")

    # and save out gene list
    data.to_csv(out_file, columns=[], header=False)
    
    return None


def convert_gene_list_to_tss(gene_list, tss_file, out_file):
    """given a gene list, get TSS
    """
    # read in TSS file to dict
    geneid_to_tss = {}
    with gzip.open(tss_file, "r") as fp:
        for line in fp:
            fields = line.strip().split()
            geneid_to_tss[fields[3]] = line
            
    # open files and write out
    tmp_file = "{}.tmp.bed.gz".format(out_file.split(".bed")[0])
    with open(gene_list, "r") as fp:
        for line in fp:
            fields = line.strip().split()
            try:
                with gzip.open(tmp_file, "a") as out:
                    out.write(geneid_to_tss[fields[0]])
            except:
                print "did not find {}".format(fields[0])
            
    # sort
    sort_file = "zcat {} | sort -k1,1 -k2,2n | gzip -c > {}".format(
        tmp_file, out_file)
    run_shell_cmd(sort_file)

    # delete tmp files
    run_shell_cmd("rm {}".format(tmp_file))

    return None


def run_rdavid(gene_list, background_gene_list, out_dir):
    """Run DAVID using R
    """
    rdavid = ("bioinformatics.go.rdavid.R {} {} {}").format(
        gene_list,
        background_gene_list,
        out_dir)
    run_shell_cmd(rdavid)
    
    return None

