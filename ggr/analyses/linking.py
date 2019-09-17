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
    region_signal_mean = np.median(np.stack(region_signals, axis=0), axis=0)
    rna_signal_mean = np.median(np.stack(rna_signals, axis=0), axis=0)
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



def reformat_interaction_file(interaction_file, out_file, remove_dups=True):
    """make it easier to work with interaction file
    """
    assert interaction_file.endswith(".gz")
    assert out_file.endswith(".gz")
    
    # reformat
    reformat = (
        "zcat {} | "
        "awk -F ',' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '\t' '{{ print $1\":\"$2\"-\"$3\"\t\"$4\"\t\"$5 }}' | "
        "gzip -c > {}").format(interaction_file, out_file)
    print reformat
    os.system(reformat)

    # and reduce
    if remove_dups:
        tmp_file = "{}.tmp.txt.gz".format(out_file.split(".txt.gz")[0])
        os.system("mv {} {}".format(out_file, tmp_file))
        seen_ids = set()
        with gzip.open(tmp_file, "r") as fp:
            with gzip.open(out_file, "w") as out:
                for line in fp:
                    fields = line.strip().split("\t")
                    reciprocal = "{}\t{}\t{}\n".format(
                        fields[1], fields[0], fields[2])
                    if reciprocal not in seen_ids:
                        # check that first is always smaller than second
                        try:
                            first_start = int(fields[0].split(":")[1].split("-")[0])
                            second_start = int(fields[1].split(":")[1].split("-")[0])
                        except:
                            print line
                        #assert first_start <= second_start, line
                        # and write out
                        out.write(line)
                        seen_ids.add(line)
        # remove tmp
        os.system("rm {}".format(tmp_file))
        
    return None


def _get_chrom_start_stop(region_id):
    """quick wrapper to get region id
    """
    chrom = region_id.split(":")[0]
    start = region_id.split(":")[1].split("-")[0]
    stop = region_id.split(":")[1].split("-")[1]
    split_id = [chrom, start, stop]
    
    return split_id

def _is_overlap(region_a, region_b):
    """check for overlap
    """
    region_a_fields = _get_chrom_start_stop(region_a)
    region_b_fields = _get_chrom_start_stop(region_b)

    print region_a_fields
    print region_b_fields
    
    if region_a_fields[0] != region_b_fields[0]:
        # chrom must match
        return False
    elif (region_b_fields[1] <= region_a_fields[1]) and (region_b_fields[2] >= region_a_fields[1]):
        # if region b starts earlier, then must end later than region a start
        return True
    elif (region_b_fields[1] >= region_a_fields[1]) and (region_b_fields[1] <= region_a_fields[2]):
        # if region b starts later, then must start earlier than region a end
        return True
    else:
        return False


def interaction_merge(data):
    """
    """
    

    
    return
    

def reconcile_interaction_mat(data):
    """determine which interactions are actually same and condense
    """
    # sort and split ids
    data = data.sort_index()
    interactions = data.index.to_series().str.split("_x_", n=2, expand=True)
    data["interaction_start"] = interactions[0]
    data["interaction_end"] = interactions[1]

    # go through each row
    current_rows = []
    current_start = None
    for example_idx in range(data.shape[0]):
        example_start = data["interaction_start"].iloc[example_idx]
        example_row = data.iloc[example_idx]
        print example_row
        
        # if current rows is empty, just put it in there
        if current_start is None:
            current_start = example_start
            current_rows.append(example_row)
            continue

        # if start region is in range of current rows, add
        if _is_overlap(current_start, example_start):
            current_rows.append(example_row)
        else:
            # new row, save out
            # and also check to make sure end parts match up too
            current_rows = pd.concat(current_rows, axis=1).transpose()
            
            print current_rows
            quit()
            
            # and replace
            current_start = example_start
            current_rows = [example_row]
    
    return



def reconcile_interactions(links_files, out_dir, method="pooled"):
    """reconcile interactions so that they'll have the same IDs
    """
    # adjust for just regions, to intersect in next step
    adj_links_files = []
    for links_file in links_files:
        adj_links_file = "{}/{}.unique.bed.gz".format(
            out_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        adjust_cmd = (
            "zcat {} | "
            "awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$3 }}' | "
            "sort -k1,1 -k2,2n | "
            "uniq | "
            "gzip -c > {}").format(
                links_file, adj_links_file)
        print adjust_cmd
        os.system(adjust_cmd)
        adj_links_files.append(adj_links_file)

    # reconcile
    reconciled_files = []
    if method == "pooled":
        # if matching to pooled links, intersect each rep
        # with pooled file to get mapping,
        # then adjust IDs
        pooled_file = adj_links_files[0]
        reconciled_files.append(links_files[0])
        for link_file_idx in range(1, len(adj_links_files)):
            # intersect w pooled file to get mapping
            links_file = adj_links_files[link_file_idx]
            tmp_file = "{}/{}.mapping.txt.gz".format(
                out_dir,
                os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
            intersect_cmd = (
                "bedtools intersect -wao -a {} -b {} | "
                "awk -F '\t' '{{ print $1\":\"$2\"-\"$3\"\t\"$4\":\"$5\"-\"$6\"\t\"$7 }}' | "
                "gzip -c > {}").format(
                    links_file, pooled_file, tmp_file)
            print intersect_cmd
            os.system(intersect_cmd)
            
            # read in the mapping
            mapping = {}
            best_overlap_lens = {}
            with gzip.open(tmp_file, "r") as fp:
                for line in fp:
                    fields = line.strip().split("\t")
                    # check if actually null mapping here?
                    try:
                        best_overlap_len = best_overlap_lens[fields[0]]
                        if int(fields[2]) > best_overlap_len:
                            mapping[fields[0]] = fields[1]
                            best_overlap_lens[fields[0]] = int(fields[2])
                    except:
                        mapping[fields[0]] = fields[1]
                        best_overlap_lens[fields[0]] = int(fields[2])

            reconcile_file = "{}/{}.matched.txt.gz".format(
                out_dir,
                os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
            with gzip.open(links_files[link_file_idx], "r") as fp:
                with gzip.open(reconcile_file, "w") as out:
                    for line in fp:
                        fields = line.strip().split("\t")
                        start_id = "{}:{}-{}".format(
                            fields[0], fields[1], fields[2])
                        end_id = fields[3].split(",")[0]
                        signal_val = fields[3].split(",")[1]
                        
                        # adjust start end of interaction
                        # but only adjust if there is a mapping val
                        if mapping[start_id] != ".:-1--1":
                            start_id = mapping[start_id]
                        
                        # adjust stop end of interaction
                        if mapping[end_id] != ".:-1--1":
                            end_id = mapping[end_id]
                            
                        # new line
                        adj_line = "{}\t{}\t{}\t{},{}\t{}\n".format(
                            start_id.split(":")[0],
                            start_id.split(":")[1].split("-")[0],
                            start_id.split(":")[1].split("-")[1],
                            end_id,
                            signal_val,
                            fields[4],
                            fields[5])
                        
                        # write out
                        out.write(adj_line)

            # save to list
            reconciled_files.append(reconcile_file)
                        
            # clean up
            os.system("rm {} {}".format(links_file, tmp_file))

        # clean up
        os.system("rm {}".format(pooled_file))
            
    return reconciled_files



def get_replicate_consistent_links(
        pooled_file, links_files, out_prefix,
        colnames=["pooled", "rep1", "rep2"]):
    """get replicate consistent links
    """
    # out dir
    out_dir = os.path.dirname(out_prefix)
    
    # map to pooled space
    matched_links_files = reconcile_interactions(
        [pooled_file] + links_files, out_dir, method="pooled")
    
    # format link files
    interaction_files = []
    for links_file in matched_links_files:
        print links_file
        tmp_file = "{}/{}.interaction_ids.tmp.txt.gz".format(
            out_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        reformat_interaction_file(links_file, tmp_file)
        interaction_files.append(tmp_file)
    print interaction_files
        
    # merge to make a table
    out_mat_file = "{}.reps.mat.txt.gz".format(out_prefix)
    summary = None
    for interaction_file_idx in range(len(interaction_files)):
        interaction_file = interaction_files[interaction_file_idx]
        data = pd.read_csv(interaction_file, sep="\t", header=None, index_col=None)
        data.index = data[0].map(str) + "_x_" + data[1].map(str)
        data[colnames[interaction_file_idx]] = data[2]
        data = data.drop([0,1,2], axis=1)

        if summary is None:
            summary = data.copy()
        else:
            summary = summary.merge(data, how="outer", left_index=True, right_index=True)

    # save out
    summary.to_csv(out_mat_file, sep="\t", compression="gzip", header=True)

    #quit()
    # only keep those with good ABC scores?
        
    # plot in R
    plot_file = "{}.reps.scatter.pdf".format(out_prefix)
    plot_cmd = "Rscript /users/dskim89/git/ggr-project/R/plot.links.rep_consistency.R {} {}".format(
        out_mat_file, plot_file)
    print plot_cmd
    os.system(plot_cmd)
    
    return
