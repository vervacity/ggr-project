#!/usr/bin/env python

import os
import sys
import h5py

import numpy as np
import pandas as pd

from ggr.analyses.linking import regions_to_genes

from tronn.plot.visualization import scale_scores
from tronn.plot.visualization import plot_weights_group
from tronn.util.formats import array_to_bed
from tronn.util.utils import DataKeys

idx_to_letter = {
    0: "A",
    1: "C",
    2: "G",
    3: "T"
}

SECOND_PASS_GENE_FILT = [
    "LAMB3",
    "LAMC2",
    "KRT5",
    "TP63",
    "EXPH5",
    "CERS3",
    "SPINK5",
    "CLDN1",
    "ELOVL4",
    "ABCA12"
    #"PEX1", # weak rna signal
    #"ERCC6", # infectious
    #"FBN1", # connective tissue (Marfans)
    #"CAST", # rna not good traj
    #"STAT1", # immune
]

THIRD_PASS_GENE_FILT = [
    "LAMB3",
    "LAMC2",
    #"KRT5", # 3rd pass filt
    "TP63",
    "EXPH5",
    "CERS3",
    #"SPINK5", # 3rd pass filt <- BRING BACK?
    "CLDN1",
    #"ELOVL4", # 3rd pass filt
    "ABCA12" # 3rd pass filt
]


def get_vignettes(
        h5_file,
        links_file,
        tss_file,
        filter_gene_set=[],
        ensembl_to_hgnc=None,
        extend_len=0,
        chromsizes=None,
        interpretation_task_indices=[0,1,2,3,4,5,6,9,10,12],
        prefix="./vignettes",
        tmp_dir="./tmp"):
    """get vignettes that are linked to genes of interest
    """
    # get a bed file from the h5 file
    metadata_bed_file = "{}/{}.metadata.tmp.bed.gz".format(
        tmp_dir, os.path.basename(prefix))
    with h5py.File(h5_file, "r") as hf:
        metadata = hf["example_metadata"][:,0]
    array_to_bed(metadata, metadata_bed_file, name_key="all", merge=False)

    # read in metadata bed file as reference
    metadata_df = pd.read_csv(metadata_bed_file, sep="\t", header=None)
    metadata_df["id"] = metadata_df[0].map(str) + ":" + metadata_df[1].map(str) + "-" + metadata_df[2].map(str)

    # regions to genes
    # NOTE since using hichip, need to extend regions 2500
    regions_to_genes_file = "{}.regions_to_genes.txt.gz".format(
        metadata_bed_file.split(".metadata")[0])
    if not os.path.isfile(regions_to_genes_file):
        regions_to_genes(
            metadata_bed_file,
            links_file,
            tss_file,
            regions_to_genes_file,
            tmp_dir=tmp_dir,
            filter_genes=filter_gene_set,
            chromsizes=chromsizes,
            extend_len=extend_len)

    # if any remaining, plot those out
    try:
        results = pd.read_csv(regions_to_genes_file, sep="\t", header=0)
    except (IOError, pd.errors.EmptyDataError, pd.io.common.EmptyDataError):
        # file is missing - no overlaps
        results = np.zeros((0,1))
    print results.shape
    total = 0
    if results.shape[0] > 0:
        print prefix
        for results_row_idx in range(results["gene_id"].shape[0]):
            gene_id = results["gene_id"].iloc[results_row_idx]
            hgnc_id = ensembl_to_hgnc[gene_id]
            bed_region_id = results["region_ids"].iloc[results_row_idx].split(",")[0]
            results_regions = metadata_df[metadata_df["id"] == bed_region_id][3].tolist()

            # check second pass
            if hgnc_id not in THIRD_PASS_GENE_FILT:
                continue

            for region_id in results_regions:
                # find region in h5 file
                h5_i = np.where(metadata == region_id)[0][0]
                metadata_string = region_id.split("features=")[-1].replace(":", "_")
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
                    plot_file = "{}.{}-{}.{}.{}.plot.pdf".format(
                        prefix, gene_id, hgnc_id, h5_i, metadata_string)
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
                    
                    # also grab actual v predicted
                    keep_indices = [0,1,2,3,4,5,6,9,10,12]
                    actual = hf["ATAC_SIGNALS.NORM"][h5_i][keep_indices]
                    predicted = hf["logits.norm"][h5_i][keep_indices]

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

                # clip
                if True:
                    flattened_vals = np.abs(interaction_impts).flatten()
                    sort_indices = np.argsort(flattened_vals)
                    thresh = flattened_vals[sort_indices][-11] # clip worst ten?
                    #thresh = np.percentile(np.abs(interaction_impts), 99.99)
                    interaction_impts[np.where(interaction_impts > thresh)] = thresh
                    interaction_impts[np.where(interaction_impts < -thresh)] = -thresh
                if True:
                    thresh = np.max(interaction_impts)
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
                
                # set up plot prefixes
                dir_name = os.path.dirname(prefix)
                root_prefix = os.path.basename(prefix)
                plot_prefix = "{}/{}.{}.{}.{}.{}.task-{}".format(
                    dir_name, hgnc_id, root_prefix, gene_id, h5_i, metadata_string, best_orig_task_idx)
                    
                # TODO plot out ATAC vs predicted
                actual_v_predicted = pd.DataFrame({
                    "timepoint": ["d0.0", "d0.5", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0"],
                    "ATAC": actual,
                    "predicted": predicted}).set_index("timepoint")
                plot_data_file = "{}.atac.actual_v_predicted.txt".format(plot_prefix)
                actual_v_predicted.to_csv(plot_data_file, header=True, index=True, sep="\t")
                plot_file = "{}.atac.actual_v_predicted.pdf".format(plot_prefix)
                fig2_dir = "/users/dskim89/git/ggr-project/figs/fig_2.modelling"
                plot_cmd = "{}/plot.atac.actual_v_pred.R {} {}".format(
                    fig2_dir, plot_data_file, plot_file)
                print plot_cmd
                os.system(plot_cmd)

                # plot interaction impts
                plot_file = "{}.interactions.plot.pdf".format(plot_prefix)
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

                # plot rna bar
                plot_file = "{}/gene_vals.{}.pdf".format(dir_name, hgnc_id)
                git_dir = "/users/dskim89/git/ggr-project/figs/fig_1.datasets"
                GGR_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a"
                rna_file = "{}/results/rna/timeseries/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz".format(GGR_DIR)
                if not os.path.isfile(plot_file):
                    plot_cmd = "{}/plot.rna_bar.R {} {} {}".format(
                        git_dir, rna_file, gene_id, plot_file)
                    print plot_cmd
                    os.system(plot_cmd)
                
                total += 1
                    
    return total


def main():
    """look for interesting vignettes that:
       - 
    """
    # args
    prefix = sys.argv[1]
    summary_file = sys.argv[2]
    synergy_dir = sys.argv[3]
    links_file = sys.argv[4]
    tss_file = sys.argv[5]
    filter_genes_file = sys.argv[6]
    chromsizes = sys.argv[7]
    out_dir = sys.argv[8]

    # set up tmp dir
    tmp_dir = "{}/tmp".format(out_dir)
    os.system("mkdir -p {}".format(tmp_dir))
    
    # read in summary file
    rules = pd.read_csv(summary_file, sep="\t", index_col=0)
    # NOTE: assumes mpra output file
    rules = rules[rules["interaction"] != "FAILED.NEGATIVE_ENDOG"]
    rules = rules[rules["interaction"] != "FAILED.TRAJ_PATTERN"]
    rules = rules["grammar"]
    print len(rules)

    # read in filter genes file
    genes_mat = pd.read_csv(filter_genes_file, sep="\t", index_col=0)
    filter_genes = genes_mat.index.values.tolist()
    # TODO just pull the mapping file
    if "genodermatoses" in filter_genes_file:
        ensembl_to_hgnc = dict(zip(genes_mat.index.values, genes_mat.iloc[:,0]))
    else:
        ensembl_to_hgnc = None
        
    # go through each rule
    total = 0
    for rule in rules:
        rule_id = os.path.basename(rule).split(".gml")[0]
        synergy_file = "{}/{}/ggr.synergy.h5".format(
            synergy_dir, rule_id)

        # debug
        if False:
            with h5py.File(synergy_file, "r") as hf:
                for key in sorted(hf.keys()): print key, hf[key].shape

        # get vignettes
        num_vignettes = get_vignettes(
            synergy_file,
            links_file,
            tss_file,
            filter_gene_set=filter_genes,
            ensembl_to_hgnc=ensembl_to_hgnc,
            #chromsizes=chromsizes,
            #extend_len=2500,
            prefix="{}/{}.{}".format(out_dir, prefix, rule_id),
            tmp_dir=tmp_dir)
        total += num_vignettes

    print total
    
    return

main()
