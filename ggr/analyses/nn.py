"""code for downstream analyses after NN interpretation
"""

import os
import h5py

import numpy as np
import pandas as pd

def get_motif_region_set(
        nn_files,
        motif_name,
        timepoint_idx,
        out_file,
        motif_idx=None,
        motif_key="sequence-weighted.active.pwm-scores.thresh.sum",
        #motif_key="sequence.active.pwm-scores.thresh.sum",
        regions_key="example_metadata"):
    """using a motif hits key, get back all regions that have that motif
    """
    # first, make sure we have motif_idx set up
    with h5py.File(nn_files[0], "r") as hf:
        #motifs = hf[motif_key].attrs["pwm_names"]
        motifs = hf["sequence-weighted.active.pwm-scores.thresh.sum"].attrs["pwm_names"]
        
    for i in range(len(motifs)):
        if motif_name in motifs[i]:
            motif_idx = i
    if motif_idx is None:
        raise ValueError, "No motif by that name found: {}".format(motif_name)

    if "weighted" not in motif_key:
        timepoint_idx = 0
        print "NOT USING WEIGHTED SCORES, ADJUSTED TIMEPOINT IDX"
    
    # then pull for all files
    for nn_file_idx in range(len(nn_files)):
        nn_file = nn_files[nn_file_idx]
        with h5py.File(nn_file, "r") as hf:
            #for key in sorted(hf.keys()): print key, hf[key].shape
            motif_scores = hf[motif_key][:,timepoint_idx,motif_idx] # {N}

            if True:
                print "NOTE USING FWD AND REV"
                num_motifs = hf[motif_key].shape[-1] / 2
                motif_scores_2 = hf[motif_key][:,timepoint_idx,motif_idx+num_motifs]
                motif_scores = motif_scores + motif_scores_2

            # use a threshold
            #thresh = np.sort(motif_scores)[-10000]
            thresh = 0
            print thresh
            
            # and then just pull those indices for regions where active
            region_indices = np.where(motif_scores > thresh)[0]
            print "num regions:", region_indices.shape
            
            # get regions
            regions = hf[regions_key][:][region_indices,0]

            # set it up as pandas dataframe
            regions = pd.DataFrame({"metadata": regions})

            # merge
            if nn_file_idx == 0:
                all_regions = regions
                all_scores = motif_scores[region_indices]
            else:
                all_regions = pd.concat([all_regions, regions])
                all_scores = np.concat([all_scores, motif_scores[region_indices]])
                
    # extract region details to bed format
    all_regions[["region", "active", "features"]] = all_regions["metadata"].str.split(
        ";", expand=True, n=3)
    all_regions["region_id"] = all_regions["region"].str.split("=", expand=True, n=2)[1]
    all_regions["score"] = all_scores

    # reduce to strongest score in region
    all_regions = all_regions[["region_id", "score"]]
    all_regions = all_regions.groupby(["region_id"]).max().reset_index()
    
    all_regions["chrom"] = all_regions["region_id"].str.split(":", expand=True, n=2)[0]
    all_regions["start-stop"] = all_regions["region_id"].str.split(":", expand=True, n=2)[1]
    all_regions[["start", "stop"]] = all_regions["start-stop"].str.split("-", expand=True, n=2)
    print "reduced set of regions:", all_regions.shape

    # TODO test thresholding here
    n_max = 30000
    if all_regions.shape[0] > n_max:
        thresh = np.sort(all_regions["score"].values)[-n_max]
        all_regions = all_regions[all_regions["score"] > thresh]
    
    # save out
    all_regions.to_csv(
        out_file,
        columns=["chrom", "start", "stop", "score"],
        header=False,
        index=False,
        sep="\t",
        compression="gzip")
    
    return


def compare_nn_gene_set_to_diff_expr(genes_file, diff_file, convert_table=None):
    """compare ranked genes list to diff results file
    """
    # ranked genes
    ranked_genes = pd.read_csv(genes_file, sep="\t")

    # add in hgnc IDs and entrez IDs using convert table
    if convert_table is not None:
        conversions = pd.read_csv(convert_table, sep="\t")
        ranked_genes = ranked_genes.merge(
            conversions, how="left", left_on="gene_id", right_on="ensembl_gene_id")
    #print ranked_genes.columns
    #print ranked_genes.shape
    
    # expr genes
    expr_genes = pd.read_csv(diff_file, index_col=0)
    expr_genes = expr_genes[expr_genes["adj.P.Val"] < 0.05]
    expr_genes = expr_genes[expr_genes["logFC"].abs() > 0.1]
    #expr_genes = expr_genes[expr_genes["logFC"] < -0.1]
    try:
        expr_cols = ["adj.P.Val", "logFC", "Gene.ID", "Gene.symbol"]
        expr_genes = expr_genes[expr_cols]
    except KeyError:
        expr_genes["Gene.symbol"] = expr_genes.index.values
    #print expr_genes.columns
    #print expr_genes.shape
    #print expr_genes

    # sort by FC
    if False:
        expr_genes["logFC.abs"] = expr_genes["logFC"].abs()
        expr_genes = expr_genes.sort_values("logFC.abs", ascending=False)
    
    # compare - full list
    #genes_a = ranked_genes["hgnc_symbol"]
    #genes_b = expr_genes["Gene.symbol"]
    #print genes_a.shape, genes_b.shape, len(set(genes_a).intersection(set(genes_b)))

    # compare - pared down to match expr list
    genes_b = expr_genes["Gene.symbol"]
    genes_a = ranked_genes["hgnc_symbol"]
    #genes_a = ranked_genes["hgnc_symbol"].iloc[0:genes_b.shape[0]]
    num_intersect = len(set(genes_a).intersection(set(genes_b)))
    
    print genes_a.shape, genes_b.shape, num_intersect
    print num_intersect / float(genes_a.shape[0])

    # compare - top 1000
    #genes_a = ranked_genes["hgnc_symbol"].iloc[0:1000]
    #genes_b = expr_genes["Gene.symbol"].iloc[0:1000]
    #print genes_a.shape, genes_b.shape, len(set(genes_a).intersection(set(genes_b)))

    # top 100
    #genes_a = ranked_genes["hgnc_symbol"].iloc[0:500]
    #genes_b = expr_genes["Gene.symbol"].iloc[0:500]
    #print genes_a.shape, genes_b.shape, len(set(genes_a).intersection(set(genes_b)))

    
    return
