#!/usr/bin/env python

import os
import sys
import h5py

import numpy as np
import pandas as pd

from tronn.plot.visualization import scale_scores
from tronn.plot.visualization import plot_weights_group
from tronn.util.formats import array_to_bed


idx_to_letter = {
    0: "A",
    1: "C",
    2: "G",
    3: "T"
}


def get_vignettes(
        h5_file,
        bed_file,
        out_dir,
        tmp_dir=".",
        importance_filter=True,
        interpretation_task_indices=[0,1,2,3,4,5,6,9,10,12],
        rsid_to_genes=None,
        ensembl_to_hgnc=None):
    """get the example indices at the locations of interest
    """
    # prefix (scanmotifs)
    dirname = h5_file.split("/")[-2]
    print dirname
    
    # get a bed file from the h5 file and overlap
    metadata_bed_file = "{}/{}.metadata.tmp.bed.gz".format(
        tmp_dir, dirname)
    with h5py.File(h5_file, "r") as hf:
        metadata = hf["example_metadata"][:,0]
    array_to_bed(metadata, metadata_bed_file, name_key="all", merge=False)

    # overlap and read back in
    overlap_file = "{}/{}.overlap.bed.gz".format(
        tmp_dir, os.path.basename(metadata_bed_file).split(".bed")[0])
    overlap_cmd = (
        "zcat {} | "
        "awk -F '\t' 'BEGIN{{OFS=\"\t\"}}{{ $2=$2+20; $3=$3-20; print }}' | " # offset
        "bedtools intersect -wo -a stdin -b {} | "
        "gzip -c > {}").format(metadata_bed_file, bed_file, overlap_file)
    print overlap_cmd
    os.system(overlap_cmd)
    try:
        overlap_data = pd.read_csv(overlap_file, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return 0
    print overlap_data.shape

    # for each, go back in and find the index, then check importance scores
    total = 0
    for overlap_i in range(overlap_data.shape[0]):
        metadata_i = overlap_data[3][overlap_i]
        metadata_string = metadata_i.split("features=")[-1].replace(":", "_")
        h5_i = np.where(metadata == metadata_i)[0][0]
        rsid_i = overlap_data[9][overlap_i]
        variant_pos = overlap_data[7][overlap_i] - overlap_data[1][overlap_i]
        variant_pos -= 1 # offset by 1

        # continue if variant not actually in active region
        if variant_pos < 0:
            continue
        
        # filter here, if using rsid to genes
        if (rsid_to_genes is not None) and (len(rsid_to_genes.keys()) > 0):
            try:
                gene_id = rsid_to_genes[rsid_i]
            except:
                continue
        else:
            gene_id = "UNKNOWN"

        # and get hgnc
        hgnc_id = ensembl_to_hgnc.get(gene_id, "UNKNOWN")
            
        # and get importance score at variant position
        with h5py.File(h5_file, "r") as hf:
            if False:
                variant_impt = hf["sequence-weighted.active"][h5_i,:,variant_pos,:]
            else:
                start_pos = max(variant_pos - 1, 0)
                stop_pos = min(variant_pos + 1, hf["sequence-weighted.active"].shape[2])
                variant_impt = hf["sequence-weighted.active"][h5_i,:,start_pos:stop_pos,:]
            variant_val = hf["sequence.active"][h5_i,:,variant_pos,:]
            variant_ref_bp = idx_to_letter[np.argmax(variant_val)]
            
        # variant scoring
        variant_impt_max = np.max(np.abs(variant_impt))
        if variant_impt_max > 0:
            
            # plot out
            if True:
                # tODO variant pos needs to bne updated!
                plot_prefix = "{}/{}-{}.{}.{}.{}.{}.{}".format(
                    out_dir, gene_id, hgnc_id, dirname, h5_i, metadata_string, rsid_i,
                    variant_ref_bp)
                # remember to add position to plot prefix after adjusments!!
                    #out_dir, dirname, gene_id, hgnc_id, h5_i, rsid_i,
                    #variant_pos, variant_ref_bp, metadata_string)

                with h5py.File(h5_file, "r") as hf:
                    #for key in sorted(hf.keys()): print key, hf[key].shape
                    mut_indices = hf["sequence-weighted.active.pwm-scores.thresh.max.idx.motif_mut"][h5_i,3]
                    
                    # get labels 
                    signals = hf["labels"][h5_i]
                    total_tasks = signals.shape[0]
                    signals = signals[interpretation_task_indices]
                    
                    # get logits
                    logits = hf["logits.motif_mut"][h5_i]
                    logits = logits[:,interpretation_task_indices]
                    
                    # best idx
                    best_task_idx = np.argmax(logits[0])
                    best_orig_task_idx = range(total_tasks)[best_task_idx]
                    
                    # filter for best idx                                        
                    logits_best_task = logits[:,best_task_idx]

                    # check and make sure the outputs is synergy (for good vignettes)
                    if not np.all(logits_best_task[0] >= logits_best_task):
                        continue
                    if not np.all(logits_best_task[3] <= logits_best_task):
                        continue
                    observed = logits_best_task[0] - logits_best_task[3]
                    expected = logits_best_task[1] + logits_best_task[2] - 2*logits_best_task[3]
                    if observed < expected:
                        continue
                    
                    # grab importances and plot
                    orig_importances = hf["sequence-weighted"][h5_i][:,420:580]
                    match_importances = hf["sequence-weighted.active"][h5_i]
                    importances = scale_scores(orig_importances, match_importances)
                    plot_file = "{}.pos-{}.impts-task.pdf".format(plot_prefix, variant_pos)
                    print plot_file
                    #plot_weights_group(importances, plot_file)
                    
                    # grab interaction importances and plot
                    orig_interact_impts = np.mean(
                        hf["sequence-weighted.motif_mut.ci"][h5_i], axis=1)
                    orig_interact_impts = np.expand_dims(orig_interact_impts, axis=-1)
                    mut_seqs = hf["sequence.motif_mut"][h5_i][:,:,420:580]
                    orig_interact_impts = np.multiply(orig_interact_impts, mut_seqs)
                    mut_mask = hf["sequence.motif_mut.mask"][h5_i] == 0
                    orig_interact_impts = np.multiply(orig_interact_impts, mut_mask)
                    orig_interact_impts = orig_interact_impts[:,best_task_idx]

                # test
                orig_interact_impts[0] = importances[best_task_idx]

                # scale appropriately
                if False:
                    match_interact_impts = hf["sequence-weighted.motif_mut"][h5_i]                    
                    match_interact_impts = match_interact_impts[:,best_task_idx]
                    interaction_impts = scale_scores(
                        orig_interact_impts,
                        match_interact_impts)
                else:
                    orig_interact_impts = np.divide(
                        orig_interact_impts,
                        np.sum(abs(orig_interact_impts), axis=(1,2), keepdims=True))
                    interaction_impts = orig_interact_impts

                # To consider - multiply by logits?
                if True:
                    logit_factors = np.reshape(logits_best_task, (-1,1,1))
                    interaction_impts = np.multiply(interaction_impts, logit_factors)

                # clip negative
                flattened_vals = np.abs(interaction_impts).flatten()
                sort_indices = np.argsort(flattened_vals)
                thresh = flattened_vals[sort_indices][-3] # clip worst two
                #thresh = np.percentile(np.abs(interaction_impts), 99.99)
                interaction_impts[np.where(interaction_impts > thresh)] = thresh
                interaction_impts[np.where(interaction_impts < -thresh)] = -thresh

                # clip the importances for better fig
                if np.abs(mut_indices[0] - mut_indices[1]) < 10:
                    continue
                clip_extend = 60
                midpoint = int(np.sum(mut_indices) / 2) - 420
                start = max(midpoint - clip_extend, 0)
                end = start + clip_extend*2
                if end > interaction_impts.shape[1]:
                    end = interaction_impts.shape[1]
                    start = end - clip_extend*2
                
                interaction_impts = interaction_impts[:,start:end]
                print interaction_impts.shape
                
                # plot interaction impts
                plot_file = "{}.pos-{}.impts-interact.task-{}.pdf".format(
                    plot_prefix, variant_pos - start, best_orig_task_idx)
                print plot_file
                plot_weights_group(interaction_impts, plot_file)

                # plot logits
                #logits_df = pd.DataFrame(np.expand_dims(logits_best_task, axis=0))
                logits_tmp_file = "{}.logits.tmp.txt.gz".format(plot_prefix)
                logits_df = pd.DataFrame({
                    "logits": logits_best_task,
                    "mut_type": ["endog", "motif,scr", "scr,motif", "scr,scr"]})
                logits_df.to_csv(logits_tmp_file, sep="\t", compression="gzip", index=False)
                git_dir = "/users/dskim89/git/ggr-project/figs/fig_3.interactions"
                plot_cmd = "{}/plot.interaction.logits.R {} {}.logits.pdf".format(
                    git_dir, logits_tmp_file, plot_prefix)
                print plot_cmd
                os.system(plot_cmd)
                
            total += 1
            
    print "TOTAL:", total
            
    return total



def main():
    """look for interesting vignettes with variants in them
    """
    # args
    gwas_file = sys.argv[1]
    out_dir = sys.argv[2]
    filter_genes_file = sys.argv[3]
    h5_dir = sys.argv[4]
    summary_file = sys.argv[5]
    tmp_dir = "{}/tmp".format(out_dir)
    os.system("mkdir -p {}".format(tmp_dir))

    # set up h5 files
    summary = pd.read_csv(summary_file, sep="\t")
    summary = summary[~summary["interaction"].str.contains("FAILED", regex=False)]
    dirnames = summary["grammar"].values
    h5_files = ["{}/{}/ggr.synergy.h5".format(h5_dir, dirname) for dirname in dirnames]
    
    # filter_genes
    genes_mat = pd.read_csv(filter_genes_file, sep="\t", index_col=0)
    filter_genes = genes_mat.index.values.tolist()

    # keep hgnc ids?
    if "genodermatoses" in filter_genes_file:
        ensembl_to_hgnc = dict(zip(genes_mat.index.values, genes_mat.iloc[:,0]))
    else:
        ensembl_to_hgnc = {}
    
    # rsid to gene
    rsid_to_genes = {}
    if "gtex" in gwas_file:
        signif_variant_gene_pairs_files = [
            "/mnt/lab_data/kundaje/users/dskim89/ggr/variants/GTEx/GTEx_Analysis_v7_eQTL/Skin_Not_Sun_Exposed_Suprapubic.v7.signif_variant_gene_pairs.txt.gz",
            "/mnt/lab_data/kundaje/users/dskim89/ggr/variants/GTEx/GTEx_Analysis_v7_eQTL/Skin_Sun_Exposed_Lower_leg.v7.signif_variant_gene_pairs.txt.gz"]

        for filename in signif_variant_gene_pairs_files:
            print filename
            pairs = pd.read_csv(filename, sep="\t")
            pairs["gene_id"] = pairs["gene_id"].str.split(".").str[0]
            print pairs.shape
            
            if True:
                pairs = pairs[pairs["gene_id"].isin(filter_genes)]
            print pairs.shape
                
            tmp_rsid_to_genes = dict(zip(pairs["variant_id"], pairs["gene_id"]))
            rsid_to_genes.update(tmp_rsid_to_genes)
    
    # go through each file
    total = 0
    for h5_file in h5_files:
        print h5_file
        file_total = get_vignettes(
            h5_file,
            gwas_file,
            out_dir,
            tmp_dir=tmp_dir,
            importance_filter=True,
            rsid_to_genes=rsid_to_genes,
            ensembl_to_hgnc=ensembl_to_hgnc)
        total += file_total
    print total
    
    return

main()
