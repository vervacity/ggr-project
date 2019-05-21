"""code for enh-gene linking
"""

import os

from collections import Counter

import numpy as np
import pandas as pd

def bed_to_gene_set_by_proximity(
        bed_file,
        tss_file,
        out_file,
        k_nearest=3,
        max_dist=250000):
    """assuming proximal linking captures the majority of linking
    """
    # bedtools closest
    tmp_file = "tss.overlap.tmp.txt"
    closest = "bedtools closest -d -k {} -a {} -b {} > {}".format(
        k_nearest, bed_file, tss_file, tmp_file)
    os.system(closest)

    # load results and use distance cutoff
    data = pd.read_table(tmp_file, header=None)
    data = data[data[9] < max_dist]

    # clean up
    data = data.sort_values([6, 9])
    data = data.drop_duplicates(6, keep="first")

    # split coumn and keep traj, gene name, distance
    gene_set = data[6].str.split(";", n=2, expand=True)
    gene_set["gene_id"] = gene_set[0].str.split("=", n=2, expand=True)[1]
    gene_set["traj"] = gene_set[1].str.split("=", n=2, expand=True)[1].astype(int)
    gene_set["dist"] = data[9]
    gene_set = gene_set[["gene_id", "traj", "dist"]]
    
    # save out
    gene_set.to_csv(out_file, compression="gzip", sep="\t", header=True, index=False)
    
    # cleanup
    os.system("rm {}".format(tmp_file))
    
    return data


def build_confusion_matrix(traj_bed_files, gene_sets, gene_clusters_file, out_file):
    """collect gene set results into a confusion matrix
    """
    # first figure out how many clusters (ATAC)
    num_gene_sets = len(gene_sets)

    # build row normalization vector
    traj_region_counts = np.zeros((num_gene_sets))
    
    # also figure out how many clusters (RNA)
    gene_clusters_data = pd.read_csv(gene_clusters_file, sep="\t")
    cluster_counts = Counter(gene_clusters_data["cluster"])
    num_gene_clusters = len(set(list(gene_clusters_data["cluster"].values)))
    
    # set up results mat
    mat = np.zeros((num_gene_sets, num_gene_clusters))

    # go through gene sets and fill out mat
    for gene_set_idx in range(len(gene_sets)):
        gene_set_file = gene_sets[gene_set_idx]
        traj_bed_file = traj_bed_files[gene_set_idx]
        num_regions = pd.read_csv(traj_bed_file, sep="\t", header=None).shape[0]
        gene_set_idx = int(
            os.path.basename(gene_set_file).split(".cluster_")[1].split(".")[0]) - 1
        gene_set = pd.read_csv(gene_set_file, sep="\t")
        gene_set_counts = Counter(gene_set["traj"])
        #print gene_set_idx, num_regions
        
        for cluster_idx, total in gene_set_counts.most_common():
            if cluster_idx == -1:
                continue
            
            col_idx = cluster_idx - 1
            #print "  ", cluster_idx, col_idx, cluster_counts[cluster_idx], total
            mat[gene_set_idx, col_idx] = total / (float(cluster_counts[cluster_idx]) * float(num_regions))

    # maybe just row normalize?
    
            
    # save out
    results = pd.DataFrame(data=mat)
    results = results.div(results.sum(axis=1), axis=0)
    
    results.to_csv(out_file, sep="\t", compression="gzip")

    return