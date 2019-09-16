"""code for enh-gene linking
"""

import os
import gzip

import numpy as np
import pandas as pd

from collections import Counter
from scipy.stats import pearsonr


def link_by_distance(
        bed_file,
        tss_file,
        out_file,
        k_nearest=3,
        max_dist=25000):
    """bed file to tss file
    """
    # bedtools closest
    tmp_file = "tss.overlap.tmp.txt"
    closest = "bedtools closest -d -k {} -a {} -b {} > {}".format(
        k_nearest, bed_file, tss_file, tmp_file)
    os.system(closest)

    # load results and use distance cutoff
    data = pd.read_csv(tmp_file, sep="\t", header=None)
    data = data[data[9] < max_dist]

    # clean up: remove chrom/start/stop of TSS
    remove_cols = [3, 4, 5]
    data = data.drop(data.columns[remove_cols], axis=1)
    
    # groupby to collapse on region
    data[9] = data[9].astype(str)
    data = data.groupby([0,1,2,7,8])[6].apply(",".join).reset_index()
    data = data[[0,1,2,6,7,8]]
    data[6] = data[0].map(str) + ":" + data[1].map(str) + "-" + data[2].map(str) + ";" + data[6]

    # save out
    data.to_csv(out_file, compression="gzip", sep="\t", header=False, index=False)
    
    # cleanup
    os.system("rm {}".format(tmp_file))
    
    return


def setup_proximity_links(
        ref_region_file,
        tss_file,
        out_links_file,
        region_signal_file=None,
        rna_signal_file=None,
        k_nearest=3,
        max_dist=25000,
        corr_coeff_thresh=0,
        pval_thresh=0.10,
        is_ggr=True):
    """actually, just set up links here?
    then can have a BED file with gene IDs
    """
    assert out_links_file.endswith(".gz")
    
    # first setup proximal links
    link_by_distance(
        ref_region_file, tss_file, out_links_file,
        k_nearest=k_nearest, max_dist=max_dist)
    
    # then read in proximal links, region signal mat, rna signal mat
    # if signal files are given
    if region_signal_file is not None:
        assert rna_signal_file is not None

        # read in links
        links = pd.read_csv(out_links_file, sep="\t", header=None)
        mappings = links[3].str.split(";", n=2, expand=True)
        links["region_id"] = mappings[0]
        links["genes"] = mappings[1]
        print links.shape
        
        # read in region signals
        region_signals = pd.read_csv(region_signal_file, sep="\t", header=0)
        if is_ggr:
            region_signals = region_signals.drop("d05", axis=1)
            
        # read in rna signals
        rna_signals = pd.read_csv(rna_signal_file, sep="\t", header=0)

        # then for each row, correlate gene to atac
        filtered_links = []
        for region_idx in range(links.shape[0]):

            if region_idx % 5000 == 0:
                print region_idx
            
            region_id = links["region_id"].iloc[region_idx]
            region_linked_genes = links["genes"].iloc[region_idx].split(",")
            region_signal = region_signals.loc[region_id,:]
            region_signal_str = ",".join(
                ["{:0.3f}".format(val) for val in region_signal.values])
            
            #filtered_genes = []
            filtered_genes = {}
            for gene_idx in range(len(region_linked_genes)):
                gene_id = region_linked_genes[gene_idx]
                try:
                    gene_signal = rna_signals.loc[gene_id,:]
                    gene_signal_str = ",".join(
                        ["{:0.3f}".format(val) for val in gene_signal.values])
                except KeyError:
                    continue
                
                # correlate
                corr_coeff, pval = pearsonr(
                    region_signal.values, gene_signal.values)

                # filters
                if corr_coeff < corr_coeff_thresh:
                    continue

                if pval > pval_thresh:
                    continue

                # append if passed filter
                filtered_genes[gene_id] = gene_signal_str
                
            # first confirm there are any genes left
            if len(filtered_genes.keys()) == 0:
                continue
            
            # pull links, replace genes, and save out
            link = links.iloc[region_idx].copy()
            link["region_signal"] = region_signal_str
            link["genes"] = ",".join([
                "{}:[{}]".format(gene_id, filtered_genes[gene_id])
                for gene_id in filtered_genes.keys()])
            filtered_links.append(link)

        # concat
        filtered_links = pd.concat(filtered_links, axis=1).transpose()
        filtered_links[3] = "region_id=" + filtered_links["region_id"].map(str) + ";region_signal=[" + filtered_links["region_signal"].map(str) + "];genes=" + filtered_links["genes"].map(str)
        filtered_links = filtered_links.drop(["region_id", "region_signal", "genes"], axis=1)
        
        # and save this out
        filtered_links.to_csv(out_links_file, sep="\t", compression="gzip", header=False, index=False)

    return


def setup_interaction_links(interaction_file, tss_file):
    """takes links from interaction format and associates them with genes
    """
    # TODO once we have ABC distance links
    
    return


def regions_to_genes_through_links(region_file, links_file, out_file):
    """intersect regions with link file to get to gene set
    """
    assert out_file.endswith(".gz")
    
    # intersect
    tmp_out_file = "links.overlap.tmp.txt.gz"
    intersect_cmd = (
        "bedtools intersect -u -a {} -b {} | "
        "gzip -c > {}").format(
            links_file, region_file, tmp_out_file)
    print intersect_cmd
    os.system(intersect_cmd)

    # summaries
    region_signals = []
    rna_signals = []
    gene_ids = []
    
    # read in results
    results = pd.read_csv(tmp_out_file, sep="\t", header=None)
    print results

    # extract data from metadata field
    for result_idx in range(results.shape[0]):
        # make a metadata dict
        metadata = results[3][result_idx].split(";")
        metadata = dict([key_val.split("=") for key_val in metadata])
        
        # adjust region signals
        metadata["region_signal"] = [
            float(val)
            for val in metadata["region_signal"].strip("[").strip("]").split(",")]
        
        # adjust gene signals
        metadata["genes"] = metadata["genes"].split("],")
        genes = {}
        for gene_str_idx in range(len(metadata["genes"])):
            gene_id, signal = metadata["genes"][gene_str_idx].split(":")
            signal = [
                float(val)
                for val in signal.strip("[").strip("]").split(",")]
            genes[gene_id] = signal
        metadata["genes"] = genes

        # add to summaries
        region_signals.append(np.array(metadata["region_signal"]))
        for gene_id in metadata["genes"].keys():
            rna_signals.append(np.array(metadata["genes"][gene_id]))

        # keep genes
        gene_ids.append(metadata["genes"].keys())
            
    # summarize
    region_signal_mean = np.mean(np.stack(region_signals, axis=0), axis=0)
    rna_signal_mean = np.mean(np.stack(rna_signals, axis=0), axis=0)
    gene_ids = np.concatenate(gene_ids, axis=0)
    gene_ids = sorted(list(set(gene_ids.tolist())))
    results = pd.DataFrame({"gene_id": gene_ids})
    
    # save out
    with gzip.open(out_file, "w") as out:
        comment_line = "#region_signal={};rna_signal={}\n".format(
            ",".join(["{:0.3f}".format(val) for val in region_signal_mean.tolist()]),
            ",".join(["{:0.3f}".format(val) for val in rna_signal_mean.tolist()]))
        out.write(comment_line)
    results.to_csv(out_file, mode="a", sep="\t", compression="gzip", header=True, index=False)

    # clean up
    os.system("rm {}".format(tmp_out_file))

    return


def traj_bed_to_gene_set_by_proximity(
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


def build_correlation_matrix(
        atac_cluster_ids,
        atac_mat_file,
        rna_cluster_ids,
        rna_mat_file,
        out_file):
    """using the mean trajectory for ATAC and RNA, make a correlation
    matrix on those patterns
    """
    # to consider - do anything with proximity?
    
    # for each cluster id set, pull mean traj out
    atac_data = pd.read_csv(atac_mat_file, sep="\t", index_col=0)
    # for GGR, need to remove d05
    atac_data = atac_data.drop("d05", axis=1)
    atac_clusters = pd.read_csv(atac_cluster_ids, sep="\t")
    atac_unique_clusters = np.sort(
        np.unique(atac_clusters["cluster"].values))
    atac_cluster_means = []
    for atac_cluster_id in atac_unique_clusters:
        # get example ids in cluster
        example_ids = atac_clusters[
            atac_clusters["cluster"] == atac_cluster_id]["id"].values
        
        # extract set from mat file
        atac_cluster_set = atac_data[atac_data.index.isin(example_ids)]
        
        # and get means
        cluster_mean = np.mean(atac_cluster_set, axis=0)
        atac_cluster_means.append(cluster_mean)

    # do the same for RNA
    rna_data = pd.read_csv(rna_mat_file, sep="\t", index_col=0)
    rna_clusters = pd.read_csv(rna_cluster_ids, sep="\t")
    rna_unique_clusters = np.sort(
        np.unique(rna_clusters["cluster"].values))
    rna_cluster_means = []
    for rna_cluster_id in rna_unique_clusters:
        # get example ids in cluster
        example_ids = rna_clusters[
            rna_clusters["cluster"] == rna_cluster_id]["id"].values
        
        # extract set from mat file
        rna_cluster_set = rna_data[rna_data.index.isin(example_ids)]
        
        # and get means
        cluster_mean = np.mean(rna_cluster_set, axis=0)
        rna_cluster_means.append(cluster_mean)
    
    # and now build out matrix
    correlation_mat = np.zeros((
        len(atac_unique_clusters),
        len(rna_unique_clusters)))
    for atac_i in range(len(atac_unique_clusters)):
        for rna_j in range(len(rna_unique_clusters)):
            # do a pearson corr
            coeff, pval = pearsonr(
                atac_cluster_means[atac_i],
                rna_cluster_means[rna_j])

            # save to mat
            correlation_mat[atac_i, rna_j] = coeff

    # and save out
    correlation_df = pd.DataFrame(
        data=correlation_mat,
        columns=rna_unique_clusters,
        index=atac_unique_clusters)
    correlation_df.to_csv(out_file, sep="\t", compression="gzip")

    # and plot
    plot_file = "{}.pdf".format(out_file.split(".txt")[0])
    plot_cmd = "plot.corr.R {} {}".format(
        out_file, plot_file)
    print plot_cmd
    
    return


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

    # and plot
    plot_file = "{}.pdf".format(out_file.split(".txt")[0])
    plot_cmd = "plot.confusion_matrix.R {} {}".format(
        out_file, plot_file)
    print plot_cmd
    
    return



def get_replicate_consistent_links(links_files, out_prefix):
    """get replicate consistent links
    """
    # format link files
    rep1_tmp_file = "{}.b1_tmp.txt.gz".format(out_prefix)
    if not os.path.isfile(rep1_tmp_file):
        reformat = (
            "zcat {} | "
            "awk -F ',' '{{ print $1\"\t\"$2 }}' | "
            "awk -F '\t' '{{ print $1\":\"$2\"-\"$3\",\"$4\"\t\"$5 }}' | "
            "gzip -c > {}").format(links_files[0], rep1_tmp_file)
        print reformat
        os.system(reformat)

    rep2_tmp_file = "{}.b2_tmp.txt.gz".format(out_prefix)
    if not os.path.isfile(rep2_tmp_file):
        reformat = (
            "zcat {} | "
            "awk -F ',' '{{ print $1\"\t\"$2 }}' | "
            "awk -F '\t' '{{ print $1\":\"$2\"-\"$3\",\"$4\"\t\"$5 }}' | "
            "gzip -c > {}").format(links_files[1], rep2_tmp_file)
        print reformat
        os.system(reformat)
    
    # load in reps and merge to make a table
    out_mat_file = "{}.reps.mat.txt.gz".format(out_prefix)
    if False:
        rep1 = pd.read_csv(rep1_tmp_file, sep="\t", header=None, index_col=0)
        rep2 = pd.read_csv(rep2_tmp_file, sep="\t", header=None, index_col=0)

        # normalize by num contacts?
        
        
        replicates = rep1.merge(rep2, how="inner", left_index=True, right_index=True)
        # normalize and filter?

        
        replicates.to_csv(out_mat_file, sep="\t", compression="gzip", header=False)

    # only keep those with good ABC scores?
        
        
    # plot in R
    plot_file = "{}.reps.scatter.pdf".format(out_prefix)
    plot_cmd = "Rscript /users/dskim89/git/ggr-project/R/plot.links.rep_consistency.R {} {}".format(
        out_mat_file, plot_file)
    print plot_cmd
    os.system(plot_cmd)
    
    return
