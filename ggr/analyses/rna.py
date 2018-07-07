# Description: RNA-specific analyses


import os
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
            cluster_file = "{}.cluster-{}.bed.gz".format(out_prefix, fields[0])
            try:
                with gzip.open(cluster_file, "a") as out:
                    out.write(geneid_to_tss[fields[1]])
            except:
                print "did not find {}".format(fields[1])
    
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

