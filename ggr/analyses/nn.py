"""code for downstream analyses after NN interpretation
"""

import os
import h5py
import glob

import numpy as np
import pandas as pd

from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve


def get_sig_motif_list(pvals_file, nn_file):
    """quick convenience fn to get motif list
    """
    # first get motifs list
    with h5py.File(nn_file, "r") as hf:
        motif_names = hf["sequence-weighted.active.pwm-scores.thresh.sum"].attrs[
            "pwm_names"]
        
    with h5py.File(pvals_file, "r") as hf:
        keys = sorted(hf["pvals"].keys())

    for key_idx in range(len(keys)):
        key = "pvals/{}/sig".format(keys[key_idx])
        with h5py.File(pvals_file, "r") as hf:
            key_sig = hf[key][:]
        if key_idx == 0:
            sig = key_sig
        else:
            sig = np.logical_or(sig, key_sig)
    motif_names = motif_names[np.where(sig)]
    
    return motif_names



def get_motif_region_set(
        nn_files,
        motif_name,
        out_file,
        timepoint_idx=None,
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
            if timepoint_idx is not None:
                motif_scores = hf[motif_key][:,timepoint_idx,motif_idx] # {N}
            else:
                motif_scores = np.max(hf[motif_key][:,:,motif_idx], axis=1) # {N}
                
            if False:
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
                all_scores = np.concatenate([all_scores, motif_scores[region_indices]])
                
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



def evaluate_enhancer_prediction(eval_files, prefix):
    """defining enhancer as a chrom state definition, figure out
    genome-wide enhancer prediction performance and plot out
    """
    for eval_file_i in range(len(eval_files)):
        eval_file = eval_files[eval_file_i]
        print eval_file
        
        with h5py.File(eval_file, "r") as hf:

            # make enhancer label - strong enhancer state of ATAC + H3K27ac + H3K4me1
            peaks_ATAC = np.any(hf["ATAC_LABELS"][:] != 0, axis=1).astype(int)
            peaks_H3K27ac = np.any(hf["H3K27ac_LABELS"][:] != 0, axis=1).astype(int)
            peaks_H3K4me1 = np.any(hf["H3K4me1_LABELS"][:] != 0, axis=1).astype(int)
            #y = np.all(
            #    [peaks_ATAC, peaks_H3K27ac, peaks_H3K4me1], axis=0)
            y = np.any(
                [peaks_H3K27ac, peaks_H3K4me1], axis=0)
            #print y.shape
            #print "positives:", np.sum(y)
            
            # baseline - just using ATAC
            signals_ATAC = hf["ATAC_SIGNALS.NORM"][:]
            #signals_H3K27ac = hf["H3K27ac_SIGNALS.NORM"][:]
            #signals_H3K4me1 = hf["H3K4me1_SIGNALS.NORM"][:]
            #X_baseline = np.concatenate([signals_ATAC, signals_H3K27ac, signals_H3K4me1], axis=1)
            X_baseline = signals_ATAC
            
            # get features
            X = hf["final_hidden"][:]

        # train test splits
        X_train_baseline, X_test_baseline, y_train, y_test = train_test_split(
            X_baseline, y, test_size=0.20, random_state=42)
        #print X_train_baseline.shape, y_train.shape, np.sum(y_train)
        #print X_test_baseline.shape, y_test.shape, np.sum(y_test)

        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.20, random_state=42)
        #print X_train.shape, y_train.shape, np.sum(y_train)
        #print X_test.shape, y_test.shape, np.sum(y_test)

        # baseline
        clf_baseline = linear_model.LogisticRegression()
        clf_baseline.fit(X_train_baseline, y_train)
        y_pred_baseline = clf_baseline.predict_proba(X_test_baseline)[:,1]
        
        pr_curve = precision_recall_curve(y_test, y_pred_baseline)
        precision, recall = pr_curve[:2]
        print auc(recall, precision)
        
        # hidden layer
        clf = linear_model.LogisticRegression()
        clf.fit(X_train, y_train)
        y_pred = clf.predict_proba(X_test)[:,1]

        pr_curve = precision_recall_curve(y_test, y_pred)
        precision, recall = pr_curve[:2]
        print auc(recall, precision)

        results = pd.DataFrame({
            "precision": precision,
            "recall": recall})
        results["fold"] = eval_file_i

        # subsample the PR curve
        if results.shape[0] > 100:
            step_size = int(results.shape[0] / 100.)
            keep_indices = np.arange(0, results.shape[0], step_size)
            keep_indices = np.append(keep_indices, results.shape[0]-1)
            results = results.iloc[keep_indices]
            results = results.append(
                pd.DataFrame(
                    [[0, 1, eval_file_i]],
                    columns=["precision", "recall", "fold"]))

        if eval_file_i == 0:
            all_results = results
        else:
            all_results = pd.concat([all_results, results], axis=0)
            
    # save out results
    out_data_file = "{}.pr_curve_data.mat.txt.gz".format(prefix)
    all_results.to_csv(out_data_file, sep="\t", compression="gzip", header=True, index=False)

    # plot with R
    r_cmd = "~/git/ggr-project/R/plot.pr_curves.R {} {}".format(out_data_file, prefix)
    print r_cmd
    os.system(r_cmd)
        
    return



def run_low_affinity_workflow(
        motifs_file,
        motif_name,
        timepoint_idx,
        chipseq_key,
        chipseq_idx,
        out_dir,
        prefix):
    """run low affinity analysis
    """
    # get motif idx
    with h5py.File(motifs_file, "r") as hf:
        motifs = hf["sequence-weighted.active.pwm-scores.thresh.sum"].attrs["pwm_names"]
        #for key in sorted(hf.keys()): print key, hf[key].shape
    motif_idx = None
    for i in range(len(motifs)):
        if motif_name in motifs[i]:
            motif_idx = i
    if motif_idx is None:
        raise ValueError, "No motif by that name found: {}".format(motif_name)
    print motifs[motif_idx]

    # then pull from the file
    with h5py.File(motifs_file, "r") as hf:
        #for key in sorted(hf.keys()): print key, hf[key].shape
        nn_scores = hf["sequence-weighted.active.pwm-scores.thresh.sum"][
            :,timepoint_idx,motif_idx]
        pwm_scores = hf["sequence.active.pwm-scores.thresh.sum"][
            :,0,motif_idx]
        labels = hf[chipseq_key][:,chipseq_idx]

    # dataframe, save these out to plot cumulative results
    results = pd.DataFrame({
        "NN": nn_scores,
        "PWM": pwm_scores,
        "labels": labels})
    scores_file = "{}/{}.scores.mat.txt.gz".format(out_dir, prefix)
    results.to_csv(
        scores_file,
        sep="\t", compression="gzip", header=True, index=False)

    # plot by ranked score
    plot_file = "{}/{}.chipseq_by_ranked_scores.pdf".format(out_dir, prefix)
    r_cmd = "~/git/ggr-project/R/plot.chipseq_by_ranked_scores.R {} {} {}".format(
        scores_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)
    
    # AUPRC
    pr_curve = precision_recall_curve(
        results["labels"].values.astype(int),        
        results["NN"].values)
    precision, recall = pr_curve[:2]
    print auc(recall, precision)

    plot_data = pd.DataFrame({
        "precision": precision,
        "recall": recall})
    plot_data["group"] = "NN"
    
    pr_curve = precision_recall_curve(
        results["labels"].values.astype(int),        
        results["PWM"].values)
    precision, recall = pr_curve[:2]
    print auc(recall, precision)
    
    plot_data_pwm = pd.DataFrame({
        "precision": precision,
        "recall": recall})
    plot_data_pwm["group"] = "PWM"

    plot_data = pd.concat([plot_data, plot_data_pwm], axis=0)
    pr_file = "{}/{}.pr_data.txt.gz".format(out_dir, prefix)
    plot_data.to_csv(pr_file, sep="\t", compression="gzip", header=True, index=False)

    # plot PR curves
    plot_file = "{}/{}.chipseq_pr_curves.pdf".format(out_dir, prefix)
    r_cmd = "~/git/ggr-project/R/plot.chipseq_pr_curves.R {} {}".format(
        pr_file, plot_file)
    print r_cmd
    os.system(r_cmd)

    return


def get_motif_scores(
        motif_name,
        nn_motif_files,
        score_key="sequence-weighted.active.pwm-scores.thresh",
        extra_key=None,
        extra_idx=0):
    """using motif name, pull out vectors
    """
    # first, make sure we have motif_idx set up
    with h5py.File(nn_motif_files[0], "r") as hf:
        motifs = hf["sequence-weighted.active.pwm-scores.thresh.sum"].attrs["pwm_names"]
    for i in range(len(motifs)):
        if motif_name in motifs[i]:
            motif_idx = i
    if motif_idx is None:
        raise ValueError, "No motif by that name found: {}".format(motif_name)
    
    # pull in metadata and motif vectors
    # use reduced motif vectors to filter down first before pulling in full motif vector
    for nn_file_idx in range(len(nn_motif_files)):
        nn_file = nn_motif_files[nn_file_idx]
        
        # get motif presence to figure out which subset to pull
        with h5py.File(nn_file, "r") as hf:
            #for key in sorted(hf.keys()): print key, hf[key].shape
            motif_present = hf["sequence.active.pwm-scores.thresh.sum"][:,0,motif_idx]
        keep_indices = np.where(motif_present != 0)[0]
        print motif_present.shape
        print keep_indices.shape
    
        # pull metadata, only keep orig region and start position of pwm hits vector
        with h5py.File(nn_file, "r") as hf:
            metadata = hf["example_metadata"][keep_indices,0]
        metadata_df = pd.DataFrame({"metadata": metadata})
        metadata_df[["region", "active", "features"]] = metadata_df["metadata"].str.split(
            ";", n=3, expand=True)
        metadata_df["region"] = metadata_df["region"].str.split("=", n=2, expand=True)[1]
        metadata_df["active"] = metadata_df["active"].str.split("=", n=2, expand=True)[1]
        metadata_df["start"] = metadata_df["active"].str.split(
            ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[0].astype(int)
        metadata_df["start"] = metadata_df["start"] + 32 # this is for the edge trim, GGR
        metadata_df = metadata_df[["region", "start"]]

        # if also looking at extra key, can pull that as well
        if extra_key is not None:
            with h5py.File(nn_file, "r") as hf:
                extra_data = hf[extra_key][keep_indices,extra_idx]
            extra_data_df = pd.DataFrame({"extra": extra_data})
            
        # pull motif vectors and reduce to max scores
        with h5py.File(nn_file, "r") as hf:
            motif_scores = hf[score_key][keep_indices,:,:,motif_idx]
        if False:
            if "weighted" in score_key:
                motif_scores = motif_scores[:,9]
            else:
                motif_scores = np.max(motif_scores, axis=1)
        else:
            motif_scores = np.max(motif_scores, axis=1)
        motif_scores_df = pd.DataFrame(motif_scores)

        # concat to dataframe
        file_motif_data = pd.concat([metadata_df, motif_scores_df], axis=1)
        if extra_key is not None:
            file_motif_data = pd.concat([file_motif_data, extra_data_df], axis=1)
        
        # merge
        if nn_file_idx == 0:
            motif_data = file_motif_data.copy()
        else:
            motif_data = pd.concat([motif_data, file_motif_data])

    return motif_data


def align_scores_to_summits(
        motif_data, master_bed_file, summit_file, aligned_file,
        seq_len=136,
        index_cols=["region", "start"]):
    """
    """
    # read in master bed file and summit file
    master_regions = pd.read_csv(master_bed_file, sep="\t", header=None)
    master_regions["region"] = master_regions[0] + ":" + master_regions[1].map(str) + "-" + master_regions[2].map(str)
    master_regions = master_regions[["region"]]

    summits = pd.read_csv(
        summit_file, sep="\t", header=None)
    summits.columns = ["chrom_summit", "start_summit", "stop_summit"]
    mappings = pd.concat([master_regions, summits], axis=1)

    # merge the mappings into the scores
    motif_data = mappings.merge(
        motif_data, how="right", left_on="region", right_on="region")
    
    # first throw out all regions that DONT contain the summit
    # TODO revisit this later
    motif_data = motif_data[
        (motif_data["start"].lt(motif_data["start_summit"])) &
        ((motif_data["start"] + seq_len).gt(motif_data["start_summit"]))
    ]
    motif_data = motif_data.reset_index(drop=True)
    
    # now align
    motif_data["start_to_summit"] = motif_data["start_summit"] - motif_data["start"]
    aligned_data = np.zeros((motif_data.shape[0], seq_len*2)) - 1
    metadata_cols = index_cols + ["chrom_summit", "start_summit", "stop_summit", "start_to_summit"]
    for row_idx in range(motif_data.shape[0]):
        # get insert indices
        adjust_num = motif_data.iloc[row_idx]["start_to_summit"]
        start_idx = seq_len - adjust_num
        stop_idx = start_idx + seq_len
        
        # place it into a background of zeros
        aligned_data[row_idx,start_idx:stop_idx] = motif_data.iloc[row_idx].drop(metadata_cols)

    # adjust metadata cols
    index_cols.remove("start")
    index_cols.append("start_summit")
        
    # now save out as dataframe
    aligned_data = pd.DataFrame(aligned_data)
    aligned_data = pd.concat([
        motif_data[index_cols],
        aligned_data], axis=1)
    aligned_data.to_csv(
        aligned_file,
        sep="\t", compression="gzip", header=True, index=False)

    return aligned_data, index_cols


def get_motif_distances(motif_data, metadata_cols, out_file):
    """go through motif data and format by distance
    """
    motif_data = motif_data.set_index(metadata_cols)
    motif_data = motif_data[motif_data.sum(axis=1) != 0]
    
    # for each row
    start_motif = []
    stop_motif = []
    distances = []
    
    for row_idx in range(motif_data.shape[0]):
        motif_vector = motif_data.iloc[row_idx,:].values
        motif_vector[motif_vector == -1] = 0

        # tmp trackers
        row_motifs = []
        row_distances = []
        
        start_idx = 0
        start_score = 0
        for pos_idx in range(1, motif_vector.shape[0]):
            motif_score = motif_vector[pos_idx]

            # save out distance and scores if actually a hit
            if motif_score != 0:

                # first also check if at beginning
                if start_idx != 0:
                    start_motif.append(start_score)
                    stop_motif.append(motif_score)
                    distances.append(pos_idx - start_idx)
                    
                # and increment
                start_idx = pos_idx
                start_score = motif_score

                # and put into tmp trackers
                row_motifs.append(motif_score)
                row_distances.append(pos_idx)
                
        # TODO - need to go through again for the longer distances
        for i in range(len(row_motifs)):
            for j in range(len(row_motifs)):
                if i > j + 1:
                    # add in the motif and distance
                    start_motif.append(row_motifs[i])
                    stop_motif.append(row_motifs[j])
                    distances.append(row_distances[i] - row_distances[j])
        
    results = pd.DataFrame({
        "start_score": start_motif,
        "stop_score": stop_motif,
        "distance": distances})
    results.to_csv(
        out_file,
        header=True, index=False, sep="\t", compression="gzip")

    return results


def get_summit_centric_motif_maps(
        motif_name, nn_motif_files, master_bed_file, summit_file, work_dir):
    """for a motif, plot the following:
        - motif hits: num hits in region vs affinity
        - motif hits: distance (between motif instances) vs affinity
        - motif hits: positions relative to ATAC summit

        - NN motif: num hits in region vs affinity
        - NN motif: distance (between motif instances) vs affinity
        - NN motif: positions relative to ATAC summit
    """
    # check for extra key: ie, if ChIP-seq data is available
    motif_name_to_key = {
        "TP53": ("TP63_LABELS", 1),
        "KLF12": ("KLF4_LABELS", 1),
        "TFAP2A": ("ZNF750_LABELS", 0)}
    extra_key = None
    extra_idx = 0
    for extra_name in motif_name_to_key.keys():
        if extra_name in motif_name:
            extra_key = motif_name_to_key[extra_name][0]
            extra_idx = motif_name_to_key[extra_name][1]

    # ==================================
    # motif HITS
    # ==================================
    motif_data_raw = get_motif_scores(
        motif_name, nn_motif_files,
        score_key="sequence.active.pwm-scores.thresh",
        extra_key=extra_key,
        extra_idx=extra_idx)
    print "Num regions (HITS): {}".format(motif_data_raw.shape[0])
    
    # adjust index cols as needed
    index_cols = ["region", "start"]
    if extra_key is not None:
        index_cols += ["extra"]
    
    # plot HITS: positions relative to ATAC summit
    aligned_file = "{}/{}.aligned.mat.txt.gz".format(work_dir, motif_name)
    motif_data_raw, index_cols = align_scores_to_summits(
        motif_data_raw, master_bed_file, summit_file, aligned_file,
        index_cols=index_cols)
    print "Num regions (HITS) after aligned to summit: {}".format(motif_data_raw.shape[0])
    plot_file = "{}/{}.aligned.hits.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.summit_aligned_motif_maps.R {} {} {}".format(
        aligned_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)
    
    # plot HITS: num hits v affinity
    data_file = "{}/{}.hits.affinities.mat.txt.gz".format(work_dir, motif_name)
    motif_data_raw.to_csv(
        data_file,
        sep="\t", header=True, index=False, compression="gzip")
    plot_file = "{}/{}.num_v_affinity.hits.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.num_motifs_vs_affinity.R {} {} {}".format(
        data_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)

    # plot HITS: distances v affinity
    distances_hits_file = "{}/{}.hits.distances_v_affinity.txt.gz".format(work_dir, motif_name)
    get_motif_distances(motif_data_raw, index_cols, distances_hits_file)
    plot_file = "{}/{}.distance_v_affinity.hits.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.distance_vs_affinity.R {} {} {}".format(
        distances_hits_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)

    if False:
        # plot HITS: positions relative to ATAC summit
        aligned_file = "{}/{}.aligned.mat.txt.gz".format(work_dir, motif_name)
        aligned_data = align_scores_to_summits(motif_data_raw, master_bed_file, summit_file, aligned_file)
        plot_file = "{}/{}.hits.aligned.pdf".format(work_dir, motif_name)
        r_cmd = "~/git/ggr-project/R/plot.summit_aligned_motif_maps.R {} {}".format(
            aligned_file, plot_file)
        print r_cmd
        #os.system(r_cmd)
    
    # ==================================
    # NN motifs
    # ==================================
    motif_data_nn = get_motif_scores(
        motif_name, nn_motif_files,
        score_key="sequence-weighted.active.pwm-scores.thresh",
        extra_key=extra_key,
        extra_idx=extra_idx)
    print "Num regions (NN): {}".format(motif_data_nn.shape[0])
    
    # adjust index cols as needed
    index_cols = ["region", "start"]
    if extra_key is not None:
        index_cols += ["extra"]

    # plot NN: positions relative to ATAC summit
    aligned_file = "{}/{}.aligned.mat.txt.gz".format(work_dir, motif_name)
    motif_data_nn, index_cols = align_scores_to_summits(
        motif_data_nn, master_bed_file, summit_file, aligned_file,
        index_cols=index_cols)
    print "Num regions (NN) after aligned to summit: {}".format(motif_data_nn.shape[0])
    plot_file = "{}/{}.aligned.nn.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.summit_aligned_motif_maps.R {} {} {}".format(
        aligned_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)

    # get the affinities (raw pwm scores) for JUST the NN-active hits
    motif_raw_tmp = motif_data_raw.set_index(index_cols)
    motif_nn_tmp = motif_data_nn.set_index(index_cols)
    motif_raw_tmp[:] = np.multiply(motif_raw_tmp.values, motif_nn_tmp.values != 0)
    motif_data_filt = motif_raw_tmp.reset_index()

    # plot NN: num hits v affinity
    data_file = "{}/{}.nn_filt_affinities.mat.txt.gz".format(work_dir, motif_name)
    motif_data_filt.to_csv(
        data_file,
        sep="\t", header=True, index=False, compression="gzip")
    plot_file = "{}/{}.num_v_affinity.nn.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.num_motifs_vs_affinity.R {} {} {}".format(
        data_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)

    # plot NN: distances vs affinity
    distances_nn_file = "{}/{}.nn.distances_v_affinity.txt.gz".format(work_dir, motif_name)
    get_motif_distances(motif_data_filt, index_cols, distances_nn_file)
    plot_file = "{}/{}.distance_v_affinity.nn.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.distance_vs_affinity.R {} {} {}".format(
        distances_nn_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)

    if False:
        # plot NN: positions relative to ATAC summit
        aligned_file = "{}/{}.aligned.mat.txt.gz".format(work_dir, motif_name)
        align_scores_to_summits(motif_data_nn, master_bed_file, summit_file, aligned_file)
        plot_file = "{}/{}.nn.aligned.pdf".format(work_dir, motif_name)
        r_cmd = "~/git/ggr-project/R/plot.summit_aligned_motif_maps.R {} {}".format(
            aligned_file, plot_file)
        print r_cmd
        #os.system(r_cmd)

    # ==================================
    # joint analysis
    # ==================================
    # TODO - look at distances from hits vs NN - where did NN results remove motifs?
    # TODO also consider using NN weighted scores here
    if True:
        distances_nn_file = "{}/{}.nn.distances_v_affinity.nn_scores.txt.gz".format(work_dir, motif_name)
        get_motif_distances(motif_data_nn, index_cols, distances_nn_file)
    
    plot_file = "{}/{}.distance_v_affinity.nn_v_hits.pdf".format(work_dir, motif_name)
    r_cmd = "~/git/ggr-project/R/plot.distance_vs_affinity.hits_v_nn.R {} {} {} {}".format(
        distances_hits_file, distances_nn_file, plot_file, motif_name)
    print r_cmd
    os.system(r_cmd)

    quit()
    
    # TODO save out a reduced bed file of filtered regions
    bed_file = "{}/{}.summit_center_filt.bed.gz".format(work_dir, motif_name)
    region_data = motif_data_nn["region"].str.split(":", n=2, expand=True)
    region_data[["start", "stop"]] = region_data[1].str.split("-", n=2, expand=True)
    region_data.to_csv(
        bed_file,
        columns=[0, "start", "stop"],
        sep="\t", header=False, index=False, compression="gzip")
    
    return


# manually curated functional terms
GOOD_GO_TERMS = [
    "stem cell differentiation",
    "hemidesmosome",
    "hair",
    "cell migration",
    "skin",
    "keratinocyte",
    "cell cycle",
    "epiderm",
    "cell junction",
    "cell proliferation",
    "adhesion",
    "lipase activity",
    "fatty acid",
    "sphingolipid",
    "glycerolipid",
    "cornif"]

REMOVE_GO_TERMS = [
    "ameboidal",
    "calcium-independent cell-cell adhesion via plasma membrane cell-adhesion molecules",
    "negative regulation of wound healing,",
    "spindle",
    "via plasma membrane",
    "adhesion molecules",
    "cell spreading",
    "cell-cell adhesion",
    "cell-cell junction",
    "endothelial",
    "angiogenesis",
    "muscle",
    "leukocyte",
    "mononuclear",
    "T cell",
    "forebrain"]

REMOVE_EXACT_TERMS = [
    "biological adhesion",
    "positive regulation of cell migration",
    "positive regulation of cell junction assembly",
    "positive regulation of epithelial cell proliferation",
    "regulation of epithelial cell proliferation",
    "calcium-independent cell-matrix adhesion",
    "positive regulation of cell proliferation",
    "regulation of cell proliferation",
    "positive regulation of cell adhesion",
    "regulation of cell adhesion",
    "regulation of cell migration",
    "regulation of cell substrate adhesion",
    "regulation of epidermal growth factor receptor signaling pathway",
    "regulation of epidermal growth factor-activated receptor activity"
]

def summarize_functional_enrichments(work_dir, out_dir):
    """take in a gprofiler outputs and make a functional map
    note: it's not the same as the combinatorial rules situation
    since can't similarly plot overlying ATAC as cleanly
    can maybe plot downstream genes?
    """
    gprofiler_files = sorted(
        glob.glob("{}/*gprofiler.txt".format(work_dir)))
    print len(gprofiler_files)
    
    for file_name_idx in range(len(gprofiler_files)):
        file_name = gprofiler_files[file_name_idx]
        motif_name = os.path.basename(file_name).split(".UNK")[0].split("_")[1]
        print motif_name

        if "SMARCA" in motif_name:
            continue
        
        results = pd.read_csv(file_name, sep="\t")
        results = results[results["domain"] == "BP"]
        
        # filtering: remove terms
        remove_indices = []
        for row_idx in range(results.shape[0]):
            index_id = results.index[row_idx]
            for remove_term in REMOVE_GO_TERMS:
                if remove_term in results.iloc[row_idx]["term.name"]:
                    remove_indices.append(index_id)
        results = results.drop(remove_indices)

        # filtering: keep terms
        keep_indices = []
        for row_idx in range(results.shape[0]):
            index_id = results.index[row_idx]
            for keep_term in GOOD_GO_TERMS:
                if keep_term in results.iloc[row_idx]["term.name"]:
                    keep_indices.append(index_id)
        results = results.loc[keep_indices,:]

        # filtering: remove exact
        remove_indices = []
        for row_idx in range(results.shape[0]):
            index_id = results.index[row_idx]
            for remove_term in REMOVE_EXACT_TERMS:
                if remove_term in results.iloc[row_idx]["term.name"]:
                    remove_indices.append(index_id)
        results = results.drop(remove_indices)
        results["-log10pval"] = -np.log10(results["p.value"].values)
        # now just keep term and pval
        results = results[["term.name", "-log10pval"]]
        results["motif"] = motif_name

        # and merge
        if file_name_idx == 0:
            all_results = results.copy()
        else:
            all_results = pd.concat([all_results, results])

    # save out all results and plot
    all_results_file = "{}/motifs.summary.go_terms.melt.txt.gz".format(out_dir)
    all_results.to_csv(
        all_results_file,
        sep="\t", compression="gzip", header=True, index=False)
    plot_file = "{}/motifs.summary.go_terms.pdf".format(out_dir)
    r_cmd = "~/git/ggr-project/R/plot.homotypic.go_terms.R {} {}".format(
        all_results_file, plot_file)
    print r_cmd
    os.system(r_cmd)
    
    return
