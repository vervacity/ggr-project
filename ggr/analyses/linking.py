"""code for enh-gene linking
"""

import os
import gzip

import numpy as np
import pandas as pd

from collections import Counter
from scipy.stats import pearsonr
from scipy.stats import fisher_exact

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.utils import quantile_norm

from ggr.util.utils import run_shell_cmd
from ggr.util.bed_utils import id_to_bed


def _add_mirrored_interactions(interaction_file, out_file):
    """if have an interaction file with just one side of the interaction,
    add in the mirrored version
    """
    # read in data
    data = pd.read_csv(interaction_file, sep="\t", header=None)

    # pull out relevant info
    interactions = data[3].str.split(",", n=2, expand=True)
    chrom = interactions[0].str.split(":", n=2, expand=True)[0]
    start = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[0]
    stop = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[1]
    score = interactions[1]

    # flip
    data_flip = data.copy()
    data_flip[3] = data[0].map(str) + ":" + data[1].map(str) + "-" + data[2].map(str) + "," + score.map(str)
    data_flip[0] = chrom
    data_flip[1] = start
    data_flip[2] = stop
    
    # concat
    data = pd.concat([data, data_flip], axis=0)

    # save out, sort, finalize
    tmp_file = "{}.unsorted.txt.gz".format(out_file.split(".txt")[0])
    data.to_csv(tmp_file, sep="\t", compression="gzip", header=False, index=False)
    sort_cmd = (
        "zcat {} | "
        "sort -k1,1 -k2,2n | "
        "uniq | "
        "gzip -c > {}").format(
            tmp_file, out_file)
    print sort_cmd
    os.system(sort_cmd)

    # clean up
    os.system("rm {}".format(tmp_file))
    
    return None


def _split_interaction_name_field(interactions):
    """take the interaction field "chr:start-stop,score"
    and return component parts as pandas Series
    """
    interactions = interactions.str.split(",", n=2, expand=True)
    chrom = interactions[0].str.split(":", n=2, expand=True)[0]
    start = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[0]
    stop = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[1]
    score = interactions[1]

    return chrom, start, stop, score


def link_by_distance(
        regions_file,
        target_regions_file,
        out_file,
        k_nearest=3,
        max_dist=25000,
        annotate_gene_ids=True,
        add_mirrored_interactions=True):
    """normally use for regulatory regions to TSS
    if annotate_gene_ids, keep track of which genes are linked
    """
    assert out_file.endswith(".gz")
    
    # bedtools closest
    tmp_file = "tss.overlap.tmp.txt"
    closest = "bedtools closest -d -k {} -a {} -b {} > {}".format(
        k_nearest, regions_file, target_regions_file, tmp_file)
    os.system(closest)

    # load results and use distance cutoff
    data = pd.read_csv(tmp_file, sep="\t", header=None)
    data = data[data[9] < max_dist]

    # build a score on the distance vals with an exponential decay
    # such that the median (or half life) of the scoring fn is 25000
    # see Gasperini Cell 2019
    decay_factor = 25000. / np.log(2) # median = decay_factor * ln(2)
    data["score"] = np.exp(-data[9] / decay_factor)
    
    # adjust to interaction format
    data["name"] = data[3].map(str) + ":" + data[4].map(str) + "-" + data[5].map(str) + "," + data["score"].map(str)
    if annotate_gene_ids:
        data["name"] = data["name"].map(str) + "," + data[6]
    data["unique_id"] = range(data.shape[0])
    keep_cols = [0, 1, 2, "name", "unique_id", 8]
    data = data[keep_cols]

    # save out
    data.to_csv(out_file, compression="gzip", sep="\t", header=False, index=False)

    # adjust if mirroring interactions
    if add_mirrored_interactions:
        tmp_mirror_file = "links.tmp.txt.gz"
        os.system("mv {} {}".format(out_file, tmp_mirror_file))
        _add_mirrored_interactions(tmp_mirror_file, out_file)
        os.system("rm {}".format(tmp_mirror_file))
        
    # cleanup
    os.system("rm {}".format(tmp_file))
    
    return


def _filter_links_for_corr(
        links,
        region_signals,
        target_region_signals,
        corr_coeff_thresh=0,
        abs_corr_coeff_thresh=None):
    """given links and signals, correlate
    """
    filtered_links = []
    for region_idx in range(links.shape[0]):
        
        if region_idx % 5000 == 0:
            print region_idx

        # get interaction, region id, target id
        interaction = links.iloc[region_idx].copy()
        region_id = "{}:{}-{}".format(
            interaction[0],
            interaction[1],
            interaction[2])
        target_id = interaction["target_id"]

        # get signals
        try:
            region_signal = region_signals.loc[region_id,:]
            target_region_signal = target_region_signals.loc[target_id]
        except KeyError:
            continue

        # correlate (note pval not reliable at low N)
        corr_coeff, _ = pearsonr(
            region_signal.values, target_region_signal.values)

        # filters
        if corr_coeff_thresh is not None:
            if corr_coeff < corr_coeff_thresh:
                continue
        if abs_corr_coeff_thresh is not None:
            if abs(corr_coeff) < abs_corr_coeff_thresh:
                continue
        
        # append if passed filter
        filtered_links.append(interaction)

    # concat, clean up
    filtered_links = pd.concat(filtered_links, axis=1).transpose()
    filtered_links[3] = filtered_links[3].str.split(",ENSG", expand=True).iloc[:,0]
    filtered_links = filtered_links.drop("target_id", axis=1)

    return filtered_links


def link_by_distance_and_corr(
        regions_file,
        target_regions_file,
        out_file,
        regions_signal_file,
        target_regions_signal_file,
        k_nearest=3,
        max_dist=100000,
        corr_coeff_thresh=0,
        abs_corr_coeff_thresh=None,
        annotate_gene_ids=True,
        add_mirrored_interactions=True,
        is_ggr=False):
    """link by distance and then filter for correlation
    """
    # link by distance
    tmp_file = "{}.tmp.txt.gz".format(out_file.split(".txt")[0])
    link_by_distance(
        regions_file,
        target_regions_file,
        tmp_file,
        k_nearest=k_nearest,
        max_dist=max_dist,
        annotate_gene_ids=True)

    # read in signals
    region_signals = pd.read_csv(
        regions_signal_file, sep="\t", header=0, index_col=0)
    if is_ggr:
        region_signals = region_signals.drop("d05", axis=1)
    target_region_signals = pd.read_csv(
        target_regions_signal_file, sep="\t", header=0, index_col=0)
    
    # read in links
    links = pd.read_csv(tmp_file, sep="\t", header=None)
    end_link_info = links[3].str.split(",", expand=True)
    if annotate_gene_ids:
        links["target_id"] = end_link_info[2]
    else:
        links["target_id"] = end_link_info[0]

    # filter
    filtered_links = _filter_links_for_corr(
        links,
        region_signals,
        target_region_signals,
        corr_coeff_thresh=corr_coeff_thresh,
        abs_corr_coeff_thresh=abs_corr_coeff_thresh)

    # and save this out
    filtered_links.to_csv(out_file, sep="\t", compression="gzip", header=False, index=False)

    # adjust if mirroring interactions
    if add_mirrored_interactions:
        tmp_mirror_file = "links.tmp.txt.gz"
        os.system("mv {} {}".format(out_file, tmp_mirror_file))
        _add_mirrored_interactions(tmp_mirror_file, out_file)
        os.system("rm {}".format(tmp_mirror_file))
    
    # clean up
    os.system("rm {}".format(tmp_file))
    
    return


def _fix_invalid_records(interaction_file, out_file):
    """adjust records for start/stop problems
    """
    with gzip.open(interaction_file, "r") as fp:
        with gzip.open(out_file, "w") as out:
            for line in fp:
                # check for broken lines
                assert line.startswith("chr")
                
                # get fields
                fields = line.strip().split("\t")

                # check for broken lines
                assert len(fields) == 6
                
                # check interaction start
                try:
                    start = int(fields[1])
                    stop = int(fields[2])
                except:
                    continue
                if start > stop:
                    fields[1] = str(stop)
                    fields[2] = str(start)

                # check interaction end
                name_fields = fields[3].split(",")
                start = int(name_fields[0].split(":")[1].split("-")[0])
                stop = int(name_fields[0].split("-")[1])
                if start > stop:
                    name_fields[0] = "{}:{}-{}".format(
                        name_fields[0].split(":")[0], stop, start)
                fields[3] = ",".join(name_fields)
                    
                # save out
                out.write("{}\n".format("\t".join(fields)))
                
    return


def reconcile_interactions(
        links_files, out_dir,
        method="pooled",
        union_min_overlap=200):
    """reconcile interactions so that they'll have the same IDs
    """
    # fix ids (start/stop) to make sure BED format compatible
    # note: interaction format
    fixed_links_files = []
    for links_file in links_files:
        fixed_links_file = "{}/{}.clean.txt.gz".format(
            out_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        _fix_invalid_records(links_file, fixed_links_file)
        fixed_links_files.append(fixed_links_file)
    
    # adjust for just regions, to intersect in next step
    # note: bed format
    adj_links_bed_files = []
    for links_file in fixed_links_files:
        adj_links_bed_file = "{}/{}.unique.bed.gz".format(
            out_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        adjust_cmd = (
            "zcat {} | "
            "awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$3 }}' | "
            "sort -k1,1 -k2,2n | "
            "uniq | "
            "gzip -c > {}").format(
                links_file, adj_links_bed_file)
        run_shell_cmd(adjust_cmd)
        adj_links_bed_files.append(adj_links_bed_file)

    # list of files with reconciled coordinates
    # note: mapping format
    reconciled_files = []
    
    # build master regions for reconciling files
    if method == "pooled":
        # pooled regions is master file
        master_regions_file = adj_links_bed_files[0]
        reconciled_files.append(fixed_links_files[0])
        adj_links_bed_files = adj_links_bed_files[1:]
        fixed_links_files = fixed_links_files[1:]
    elif method == "union":
        # make a union file of all seen regions with an overlap of 100 bp minimum
        # TODO adjust this?
        master_regions_file = "{}/master.union.bed.gz".format(out_dir)
        union_cmd = (
            "zcat {} | "
            "sort -k1,1 -k2,2n | "
            "bedtools merge -d -{} -i stdin | "
            "gzip -c > {}").format(
                " ".join(adj_links_bed_files),
                union_min_overlap,
                master_regions_file)
        run_shell_cmd(union_cmd)

    # for each file, map regions
    for link_file_idx in range(len(adj_links_bed_files)):

        # in file, out file
        links_file = adj_links_bed_files[link_file_idx]
        reconcile_file = "{}/{}.reconciled.txt.gz".format(
            out_dir,
            os.path.basename(links_file).split(".clean")[0])

        # intersect to master regions        
        tmp_file = "{}/{}.mapping.txt.gz".format(
            out_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        intersect_cmd = (
            "bedtools intersect -wao -a {} -b {} | "
            "awk -F '\t' '{{ print $1\":\"$2\"-\"$3\"\t\"$4\":\"$5\"-\"$6\"\t\"$7 }}' | "
            "gzip -c > {}").format(
                links_file, master_regions_file, tmp_file)
        run_shell_cmd(intersect_cmd)

        # read in the mapping, keep the MOST overlapping hit
        mapping = {}
        best_overlap_lens = {}
        with gzip.open(tmp_file, "r") as fp:
            for line in fp:
                fields = line.strip().split("\t")
                try:
                    best_overlap_len = best_overlap_lens[fields[0]]
                    if int(fields[2]) > best_overlap_len:
                        mapping[fields[0]] = fields[1]
                        best_overlap_lens[fields[0]] = int(fields[2])
                except:
                    mapping[fields[0]] = fields[1]
                    best_overlap_lens[fields[0]] = int(fields[2])

        # map from original coords to reconciled coords
        with gzip.open(fixed_links_files[link_file_idx], "r") as fp:
            with gzip.open(reconcile_file, "w") as out:
                for line in fp:
                    # split to fields
                    fields = line.strip().split("\t")
                    start_id = "{}:{}-{}".format(
                        fields[0], fields[1], fields[2])
                    end_id = fields[3].split(",")[0]
                    signal_val = fields[3].split(",")[1]

                    # map to reconciled coords
                    # but only adjust if there is a mapping val
                    if mapping[start_id] != ".:-1--1":
                        start_id = mapping[start_id]
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
        
    return reconciled_files


def run_interaction_idr_test(mat_file, out_file, idr_thresh=0.10):
    """convert to narrowpeak format files, run IDR,
    filter for consistent peaks
    """
    # make tmp files
    rep1_tmp_file = "{}.rep1.bed.gz".format(mat_file.split(".txt")[0])
    rep2_tmp_file = "{}.rep2.bed.gz".format(mat_file.split(".txt")[0])
    
    # read in data and adjust
    data = pd.read_csv(mat_file, sep="\t", index_col=0)
    data["region_id"] = data.index
    interactions = data["region_id"].str.split("_x_", n=2, expand=True)
    # CHROM HACK: to ensure uniqueness, use to map back also
    data["chrom"] = ["chr{}".format(str(val)) for val in range(data.shape[0])] 
    data["start"] = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[0]
    data["stop"] = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[1]
    data["interaction"] = interactions[1].map(str) + "," + data["pooled"].map(str)

    # filter for unique set to put into IDR
    region_ids = data.index.values.tolist()
    seen_ids = set()
    for region_id in region_ids:
        # flip and check if in seen_ids
        fields = region_id.split("_x_")
        flipped_id = "{}_x_{}".format(fields[1], fields[0])
        if flipped_id in seen_ids:
            continue
        seen_ids.add(region_id)
    data_idr = data[data.index.isin(seen_ids)]
    
    # rep1
    rep1_cols = ["chrom", "start", "stop", "interaction", "rep1"]
    rep1 = data_idr[rep1_cols].copy()
    rep1["strand"] = "."
    rep1["plotStart"] = rep1["start"]
    rep1["plotEnd"] = rep1["stop"]
    rep1["rgb"] = 100
    rep1.to_csv(rep1_tmp_file, sep="\t", header=None, index=None, compression="gzip")

    # rep2
    rep2_cols = ["chrom", "start", "stop", "interaction", "rep2"]
    rep2 = data_idr[rep2_cols].copy()
    rep2["strand"] = "."
    rep2["plotStart"] = rep2["start"]
    rep2["plotEnd"] = rep2["stop"]
    rep2["rgb"] = 100
    rep2.to_csv(rep2_tmp_file, sep="\t", header=None, index=None, compression="gzip")

    # IDR
    idr_tmp_file = "{}.idr_results.tmp".format(out_file.split(".gz")[0])
    idr_cmd = (
        "idr --samples {} {} "
        "--input-file-type bed "
        "--rank 5 "
        "-o {} "
        "--idr-threshold {} "
        "--plot").format(
            rep1_tmp_file, rep2_tmp_file,
            idr_tmp_file,
            idr_thresh)
    run_shell_cmd(idr_cmd)

    # read in IDR file and use IDs to filter mat file
    idr_results = pd.read_csv(idr_tmp_file, sep="\t", header=None)
    keep_ids = data[data["chrom"].isin(idr_results.iloc[:,0])].index.to_series()
    flip_ids = keep_ids.str.split("_x_", n=2, expand=True)
    flip_ids["flipped"] = flip_ids[1].map(str) + "_x_" + flip_ids[0].map(str)
    keep_ids = pd.concat([keep_ids, flip_ids["flipped"]], axis=0)
    data_filt = data[data.index.isin(keep_ids)]
    
    # clean up and save out
    drop_cols = ["chrom", "start", "stop", "region_id", "interaction"]
    data_filt = data_filt.drop(drop_cols, axis=1)
    data_filt.to_csv(out_file, sep="\t", compression="gzip", header=True, index=True)

    # plot
    plot_file = "{}.after_idr.scatter.png".format(out_file.split(".txt")[0])
    plot_cmd = "Rscript /users/dskim89/git/ggr-project/R/plot.links.rep_consistency.R {} {}".format(
        out_file, plot_file)
    print plot_cmd
    os.system(plot_cmd)

    # clean up (delete tmp)
    os.system("rm {} {}".format(rep1_tmp_file, rep2_tmp_file))
    
    return None


def _convert_to_coord_format(interaction_file, out_file):
    """make it easier to work with interaction file,
    produces start_id, end_id, score
    """
    assert interaction_file.endswith(".gz")
    assert out_file.endswith(".gz")
    
    # reformat
    reformat = (
        "zcat {} | "
        "awk -F ',' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '\t' '{{ print $1\":\"$2\"-\"$3\"\t\"$4\"\t\"$5 }}' | "
        "gzip -c > {}").format(interaction_file, out_file)
    run_shell_cmd(reformat)
        
    return None


def _convert_to_interaction_format(mat_file, out_file, score_val="pooled"):
    """convert to interaction format
    """
    # read in data
    data = pd.read_csv(mat_file, sep="\t", index_col=0)

    # which method to score with
    if score_val == "pooled":
        pass
    elif score_val == "max":
        data = data.fillna(0)
        data["pooled"] = data.values.max(axis=1)
    else:
        raise ValueError, "not implemented!"

    # adjust
    data["region_id"] = data.index
    interactions = data["region_id"].str.split("_x_", n=2, expand=True)
    data["chrom"] = interactions[0].str.split(
        ":", n=2, expand=True)[0]
    data["start"] = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[0]
    data["stop"] = interactions[0].str.split(
        ":", n=2, expand=True)[1].str.split("-", n=2, expand=True)[1]
    data["interaction"] = interactions[1].map(str) + "," + data["pooled"].map(str)
        
    # clean up
    keep_cols = ["chrom", "start", "stop", "interaction"]
    data = data[keep_cols]
    data["unique_id"] = range(data.shape[0])
    data["strand"] = "."
    data = data.drop_duplicates()
    
    # save out
    data.to_csv(out_file, sep="\t", compression="gzip", header=False, index=False)
    
    return


def get_replicate_consistent_links(
        pooled_file, links_files, out_prefix,
        idr_thresh=0.10,
        colnames=["pooled", "rep1", "rep2"],
        add_mirrored_interactions=False):
    """get replicate consistent links, using pooled results as a guide
    note: includes flipped interactions (start_x_end and end_x_start)
    """
    # out dir, tmp_dir, results file
    out_dir = os.path.dirname(out_prefix)
    tmp_dir = "{}/tmp".format(out_dir)
    os.system("mkdir -p {}".format(tmp_dir))
    out_mat_file = "{}.overlap.mat.txt.gz".format(out_prefix)

    # mirror here if needed
    if add_mirrored_interactions:
        orig_pooled_file = pooled_file
        pooled_file = "{}/{}.add_mirrored.txt.gz".format(
            tmp_dir,
            os.path.basename(orig_pooled_file).split(".bed")[0].split(".txt")[0])
        _add_mirrored_interactions(orig_pooled_file, pooled_file)

        for links_file_idx in range(len(links_files)):
            orig_link_file = links_files[links_file_idx]
            link_file = "{}/{}.add_mirrored.txt.gz".format(
                tmp_dir,
                os.path.basename(orig_link_file).split(".bed")[0].split(".txt")[0])
            _add_mirrored_interactions(orig_link_file, link_file)
            links_files[links_file_idx] = link_file
    
    # first reconcile interaction coords, mapping to pooled links
    matched_links_files = reconcile_interactions(
        [pooled_file] + links_files, tmp_dir, method="pooled")

    # format link files to "coordinate format" (start_id, end_id, score)
    interaction_files = []
    for links_file in matched_links_files:
        tmp_file = "{}/{}.coord_formatted.tmp.txt.gz".format(
            tmp_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        _convert_to_coord_format(links_file, tmp_file)
        interaction_files.append(tmp_file)

    # merge to make a table
    summary = None
    for interaction_file_idx in range(len(interaction_files)):
        # adjust data
        interaction_file = interaction_files[interaction_file_idx]
        data = pd.read_csv(interaction_file, sep="\t", header=None, index_col=None)
        data.index = data[0].map(str) + "_x_" + data[1].map(str)
        data[colnames[interaction_file_idx]] = data[2]
        data = data.drop([0,1,2], axis=1)

        # merge. use INTERSECTION since the peak needs to be present
        # in pooled AND rep1 AND rep2
        if summary is None:
            summary = data.copy()
        else:
            summary = summary.merge(
                data, how="inner", left_index=True, right_index=True)

    # clean up duplicates (multiple matching hits per link)
    # first flip all to make sure both sides are represented
    summary_flipped = summary.copy()
    flip_regions = summary_flipped.index.to_series().str.split("_x_", n=2, expand=True)
    flip_regions["flip"] = flip_regions[1].map(str) + "_x_" + flip_regions[0].map(str)
    summary_flipped.index = flip_regions["flip"]
    summary = pd.concat([summary, summary_flipped], axis=0)

    # then choose max pooled value to keep
    # note max can produce multiple hits, so need to drop dups in index as final step
    summary = summary.reset_index()
    keep_rows = summary.groupby(["index"])["pooled"].transform(max) == summary["pooled"]
    summary = summary[keep_rows]
    summary = summary.set_index("index")
    summary = summary.loc[~summary.index.duplicated(keep="first")]

    # save out    
    summary.to_csv(out_mat_file, sep="\t", compression="gzip", header=True)

    # plot in R
    plot_file = "{}.reps.scatter.png".format(out_prefix)
    plot_cmd = "Rscript /users/dskim89/git/ggr-project/R/plot.links.rep_consistency.R {} {}".format(
        out_mat_file, plot_file)
    run_shell_cmd(plot_cmd)

    # and set threshold for consistency - use IDR framework
    mat_filt_file = "{}.idr_filt.mat.txt.gz".format(out_mat_file.split(".mat")[0])
    run_interaction_idr_test(out_mat_file, mat_filt_file, idr_thresh=idr_thresh)
        
    # and save out as interaction format
    interaction_file = "{}.interactions.txt.gz".format(mat_filt_file.split(".mat")[0])
    _convert_to_interaction_format(mat_filt_file, interaction_file)
    
    return


def get_timepoint_consistent_links(
        links_files, out_prefix, method="union",
        run_quantile_norm=True,
        colnames=["d0", "d3", "d6"]):
    """aggregate links across timepoints
    """
    # out dir, tmp_dir
    out_dir = os.path.dirname(out_prefix)
    tmp_dir = "{}/tmp".format(out_dir)
    os.system("mkdir -p {}".format(tmp_dir))
    out_mat_file = "{}.overlap.mat.txt.gz".format(out_prefix)
    
    # map to union
    matched_links_files = reconcile_interactions(
        links_files, tmp_dir, method=method)

    # format link files
    interaction_files = []
    for links_file in matched_links_files:
        tmp_file = "{}/{}.coord_formatted.tmp.txt.gz".format(
            tmp_dir,
            os.path.basename(links_file).split(".bed")[0].split(".txt")[0])
        _convert_to_coord_format(links_file, tmp_file)
        interaction_files.append(tmp_file)

    # merge to make a table
    summary = None
    for interaction_file_idx in range(len(interaction_files)):
        # adjust data
        interaction_file = interaction_files[interaction_file_idx]
        data = pd.read_csv(interaction_file, sep="\t", header=None, index_col=None)
        data.index = data[0].map(str) + "_x_" + data[1].map(str)
        data[colnames[interaction_file_idx]] = data[2]
        data = data.drop([0,1,2], axis=1)

        # merge. use UNION across timepoints
        if summary is None:
            summary = data.copy()
        else:
            summary = summary.merge(
                data, how="outer", left_index=True, right_index=True)

    # fill na
    summary = summary.fillna(0)
    
    # quantile norm
    if run_quantile_norm:
        summary[:] = quantile_norm(summary.values)
    
    # clean up duplicates (multiple matching hits per link)
    # first flip all to make sure both sides are represented

    summary["max"] = summary.values.max(axis=1)
    summary_flipped = summary.copy()
    flip_regions = summary_flipped.index.to_series().str.split("_x_", n=2, expand=True)
    flip_regions["flip"] = flip_regions[1].map(str) + "_x_" + flip_regions[0].map(str)
    summary_flipped.index = flip_regions["flip"]
    summary = pd.concat([summary, summary_flipped], axis=0)

    # then choose max pooled value to keep
    # note max can produce multiple hits, so need to drop dups in index as final step
    summary = summary.reset_index()
    keep_rows = summary.groupby(["index"])["max"].transform(max) == summary["max"]
    summary = summary[keep_rows]
    summary = summary.set_index("index")
    summary = summary.loc[~summary.index.duplicated(keep="first")]
    summary = summary.drop("max", axis=1)
    
    # save out
    summary.to_csv(out_mat_file, sep="\t", compression="gzip", header=True)
    
    # and save out as interaction format
    interaction_file = "{}.interactions.txt.gz".format(out_mat_file.split(".mat")[0])
    _convert_to_interaction_format(out_mat_file, interaction_file, score_val="max")
    
    return


def _get_aggregate_score(region_ids_string):
    """get back aggregate score
    """
    region_ids_w_score = region_ids_string.split(";")
    scores = [float(val.split(",")[1]) for val in region_ids_w_score]

    return np.sum(scores)


def regions_to_genes(
        region_file, links_file, tss_file, out_file,
        tmp_dir=None,
        filter_by_score=None,
        filter_genes=[],
        chromsizes=None,
        extend_len=0):
    """goes from a region set to a linked gene set
    note that region id file must match ids in region signals if using
    """
    if extend_len > 0:
        assert chromsizes is not None
    
    if tmp_dir is None:
        tmp_dir = os.path.dirname(out_file)
    
    # overlap regions with links
    # keep original region info to pass through
    tmp_file = "{}/{}.region_to_region.tmp.txt.gz".format(
        tmp_dir,
        os.path.basename(out_file).split(".txt")[0])
    intersect_cmd = (
        "bedtools intersect -wo -a {} -b {} | "
        "gzip -c > {}").format(
            links_file, region_file, tmp_file)
    run_shell_cmd(intersect_cmd)
    
    # reformat linked regions to bed file
    # keep original region info to pass through
    tmp_bed_file = "{}.bed.gz".format(tmp_file.split(".txt")[0])
    data = pd.read_csv(tmp_file, sep="\t", header=None)
    chrom, start, stop, score = _split_interaction_name_field(data[3])
    region_id = data[6].map(str) + ":" + data[7].map(str) + "-" + data[8].map(str) + "," + score.map(str)
    bed_data = pd.DataFrame(
        {"chrom": chrom,
         "start": start,
         "stop": stop,
         "region_id": region_id})

    # groupby region ids linking to that position
    groupby_cols = ["chrom", "start", "stop"]
    bed_data = bed_data.groupby(groupby_cols)["region_id"].apply(list)
    bed_data = bed_data.reset_index()
    bed_data["region_id"] = bed_data["region_id"].apply(";".join)
    bed_data["unique_id"] = range(bed_data.shape[0])
    bed_data["strand"] = "."    
    bed_data = bed_data[["chrom", "start", "stop", "region_id", "unique_id", "strand"]]
    bed_data.to_csv(
        tmp_bed_file, sep="\t", compression="gzip",
        header=False, index=False)
    
    # overlap linked regions with TSS file (generally, all expressed)
    # extend overlap distance as needed
    tmp_genes_file = "{}/{}.region_to_tss.tmp.txt.gz".format(
        tmp_dir,
        os.path.basename(out_file).split(".txt")[0])
    if extend_len > 0:
        intersect_cmd = (
            "bedtools slop -i {0} -g {1} -b {2} | "
            "bedtools intersect -wo -a {3} -b - | "
            "gzip -c > {4}").format(
                tmp_bed_file, chromsizes, extend_len, tss_file, tmp_genes_file)
    else:
        intersect_cmd = (
            "bedtools intersect -wo -a {} -b {} | "
            "gzip -c > {}").format(
                tss_file, tmp_bed_file, tmp_genes_file)
    run_shell_cmd(intersect_cmd)
    
    # read in and clean up
    try:
        genes = pd.read_csv(tmp_genes_file, sep="\t", header=None)
    except pd.errors.EmptyDataError:
        return None
    genes = genes[[3,9]]
    genes.columns = ["gene_id", "region_ids"]
    genes = genes.groupby("gene_id")["region_ids"].apply(list).reset_index()
    genes["region_ids"] = genes["region_ids"].apply(";".join)
    genes["score"] = genes["region_ids"].apply(_get_aggregate_score)
    genes = genes.sort_values("score", ascending=False)

    # filtering (as needed) here
    # could filter for dynamic genes only
    if len(filter_genes) != 0:
        print "before filter genes:", genes.shape[0]
        genes = genes[genes["gene_id"].isin(filter_genes)]
        print "after filter genes:", genes.shape[0]

    # filter for score as needed
    if filter_by_score is not None:
        print "before filter score:", genes.shape[0]
        genes = genes[genes["score"] > filter_by_score]
        print "after filter score:", genes.shape[0]
        
    # save out
    print "genes linked: ", genes.shape[0]
    genes.to_csv(out_file, sep="\t", compression="gzip", header=True, index=False)

    # clean up
    os.system("rm {} {} {}".format(tmp_file, tmp_bed_file, tmp_genes_file))
    
    return None


def regions_to_genes_w_correlation_filtering(
        region_file, links_file, tss_file, out_file,
        region_signal_mat,
        rna_signal_mat,
        tmp_dir=None,
        filter_by_score=None,
        filter_genes=[],
        chromsizes=None,
        extend_len=0,
        corr_thresh=0,
        pval_thresh=1):
    """extends regions to genes with correlation filtering
    """
    # set up
    if tmp_dir is None:
        tmp_dir = os.path.dirname(out_file)

    # get region ids
    region_ids = list(set(
        pd.read_csv(region_file, sep="\t", header=None)[3].values.tolist()))
        
    # first run regions_to_genes and read in
    tmp_out_file = "{}/{}.linked_genes.txt.gz".format(
        tmp_dir, os.path.basename(region_file).split(".bed")[0].split(".txt")[0])
    regions_to_genes(
        region_file, links_file, tss_file, tmp_out_file,
        tmp_dir=tmp_dir,
        filter_by_score=filter_by_score,
        filter_genes=filter_genes,
        chromsizes=chromsizes,
        extend_len=extend_len)

    # read in - bug fix here for bedtools error
    if not os.path.isfile(tmp_out_file):
        return None, None
    linked_genes = pd.read_csv(tmp_out_file, sep="\t", header=0, index_col=0)

    # filter vs atac signal pattern
    region_signal = region_signal_mat.loc[region_ids]
    rna_signal = rna_signal_mat.loc[linked_genes.index.values.tolist()]

    # get corr and pval
    keep_genes = pd.DataFrame(
        rna_signal.apply(
            lambda x: pearsonr(x, region_signal.mean(axis=0)),
            axis=1).values.tolist(),
        columns=["corr", "pval"],
        index=rna_signal.index)

    # threshold
    keep_genes = keep_genes[keep_genes["corr"] > corr_thresh]
    keep_genes = keep_genes[keep_genes["pval"] <= pval_thresh]

    # filter
    rna_signal = rna_signal[rna_signal.index.isin(keep_genes.index)]
    linked_genes = linked_genes[linked_genes.index.isin(keep_genes.index)]
    linked_genes = linked_genes.sort_values("score", ascending=False)

    # save out
    linked_genes.reset_index().to_csv(
        out_file, columns=["gene_id"], index=False, compression="gzip")

    # clean up
    os.system("rm {}".format(tmp_out_file))
    
    return linked_genes, rna_signal
    

def region_clusters_to_genes(
        cluster_file,
        links_file,
        tss_file,
        out_dir,
        region_signal_file=None,
        rna_signal_file=None,
        filter_gene_file=None,
        filter_by_signal_corr=False,
        filter_by_score=None,
        run_enrichments=True,
        background_gene_file=None,
        chromsizes=None,
        extend_len=0,
        is_ggr=True):
    """take a cluster file, split to present clusters and get gene sets
    including signal files will output mean signals for regions and signals

    possible filters:
    - filter by genes in gene file (must be first col)
    - filter also for correlation
    """
    # set up tmp dir
    tmp_dir = "{}/tmp".format(out_dir)
    run_shell_cmd("mkdir -p {}".format(tmp_dir))
    
    # read in cluster file
    data = pd.read_csv(cluster_file, sep="\t")
    unique_cluster_ids = sorted(list(set(data["cluster"].values.tolist())))

    # read in region signals
    if region_signal_file is not None:
        region_signals = pd.read_csv(region_signal_file, sep="\t", index_col=0)
        if is_ggr:
            region_signals = region_signals.drop("d05", axis=1)
    else:
        region_signals = None

    # read in rna signals
    if rna_signal_file is not None:
        rna_signals = pd.read_csv(rna_signal_file, sep="\t", index_col=0)
    else:
        rna_signals = None

    # filter gene set if filter file present
    if filter_gene_file is not None:
        filter_genes = pd.read_csv(
            filter_gene_file, sep="\t", index_col=0).index.values.tolist()
    else:
        filter_genes = []
        
    # go through each cluster
    for cluster_id in unique_cluster_ids:
        # set up bed file
        tmp_id_file = "{}/cluster-{}.txt.gz".format(tmp_dir, cluster_id)
        tmp_bed_file = "{}/cluster-{}.bed.gz".format(tmp_dir, cluster_id)
        ids_in_cluster = data[data["cluster"] == cluster_id]["id"]
        print "num regions in cluster:", ids_in_cluster.shape[0]
        ids_in_cluster.to_csv(
            tmp_id_file,
            sep="\t", compression="gzip",
            header=False, index=False)
        id_to_bed(tmp_id_file, tmp_bed_file, sort=True)

        # if filter by signal corr, figure out which genes match
        # trajectory (loosely) and use to filter gene set
        if filter_by_signal_corr:
            assert region_signals is not None
            assert rna_signals is not None
            cluster_mean_signal = region_signals[
                region_signals.index.isin(ids_in_cluster)].mean(axis=0)
            rna_corrs = rna_signals.corrwith(cluster_mean_signal, axis=1)
            filter_genes = rna_corrs[rna_corrs > 0].index.values.tolist()
        
        # bed file to gene set
        gene_set_file = "{}/cluster-{}.regions_to_genes.txt.gz".format(
            out_dir, cluster_id)
        regions_to_genes(
            tmp_bed_file,
            links_file,
            tss_file,
            gene_set_file,
            tmp_dir=tmp_dir,
            filter_genes=filter_genes,
            filter_by_score=filter_by_score,
            chromsizes=chromsizes,
            extend_len=extend_len)
        
        # run gprofiler
        gprofiler_dir = "{}/cluster-{}.enrichments".format(out_dir, cluster_id)
        run_shell_cmd("mkdir -p {}".format(gprofiler_dir))
        if run_enrichments:
            run_gprofiler(
                gene_set_file, background_gene_file,
                gprofiler_dir, ordered=True, header=True)

        # TODO: add in mean ATAC and mean RNA traj vals
        

    return


def build_confusion_matrix(
        atac_clusters_file,
        rna_clusters_file,
        traj_dir,
        prefix,
        normalize=True):
    """pull together confusion matrix
    """
    # read in atac data and get clusters
    atac = pd.read_csv(atac_clusters_file, sep="\t", header=0)
    atac_cluster_ids = sorted(
        list(set(atac["cluster"].values.tolist())))

    # read in rna data and get clusters
    rna = pd.read_csv(rna_clusters_file, sep="\t", header=0)
    all_rna_ids = rna["id"].values.tolist()
    num_total_genes = rna.shape[0]
    rna = rna.groupby("cluster")["id"].apply(list)
    rna_cluster_ids = rna.index.values.tolist()
    
    # go through clusters
    results = np.zeros(
        (len(atac_cluster_ids), len(rna_cluster_ids)))
    for atac_id in atac_cluster_ids:
        # pull regions to genes file and make sure only using genes that are in rna cluster file
        gene_file = "{}/cluster-{}.regions_to_genes.txt.gz".format(traj_dir, atac_id)
        genes = pd.read_csv(gene_file, sep="\t", index_col=0)
        genes = genes[genes.index.isin(all_rna_ids)]
        num_atac_regions = np.sum(atac["cluster"] == atac_id)
        #print "ATAC:", atac_id, num_atac_regions
        
        for rna_id in rna_cluster_ids:
            linked_genes = genes[genes.index.isin(rna.loc[rna_id])]
            num_possible_genes = len(rna.loc[rna_id])
            
            # calculate fisher exact test (counts may be low, so dont use chi square)
            atac_pos_rna_pos = linked_genes.shape[0]
            atac_pos_rna_neg = genes.shape[0] - atac_pos_rna_pos
            atac_neg_rna_pos = num_possible_genes - atac_pos_rna_pos
            atac_neg_rna_neg = num_total_genes - (
                atac_pos_rna_pos + atac_pos_rna_neg + atac_neg_rna_pos)
            table = [
                [atac_pos_rna_pos, atac_pos_rna_neg],
                [atac_neg_rna_pos, atac_neg_rna_neg]] # note doesn't matter row v column
            stat, pval = fisher_exact(table, alternative="greater") # stat is an odds ratio
            
            # normalized
            normalize = False
            if normalize:
                enrichment = linked_genes.shape[0] / (
                    float(num_atac_regions) * float(num_possible_genes))
            else:
                enrichment = linked_genes.shape[0] / float(num_possible_genes)
                
            #print "    RNA:", rna_id, linked_genes.shape[0], num_possible_genes, enrichment, -np.log10(pval), stat
            #print table
            results[atac_id-1,rna_id-1] = stat #-np.log10(pval) # stat
            
    # save out results
    results = pd.DataFrame(data=results)
    if False:
        results = results.div(results.sum(axis=0), axis=1) # normalize columns (gene count and linking expectations)
        results = results.div(results.sum(axis=1), axis=0) # normalize rows (prob across rows)
    else:
        results = results.div(results.sum(axis=1), axis=0) # normalize rows (prob across rows)
        results = results.div(results.sum(axis=0), axis=1) # normalize columns (gene count and linking expectations)
    
    out_file = "{}.mat.txt.gz".format(prefix)
    results.to_csv(out_file, sep="\t", compression="gzip")

    # and plot
    plot_file = "{}.pdf".format(prefix)
    plot_cmd = "plot.confusion_matrix.R {} {}".format(
        out_file, plot_file)
    print plot_cmd
    run_shell_cmd(plot_cmd)
    
    return


def build_correlation_matrix(
        atac_cluster_ids,
        atac_mat_file,
        rna_cluster_ids,
        rna_mat_file,
        out_file,
        is_ggr=True):
    """using the mean trajectory for ATAC and RNA, make a correlation
    matrix on those patterns
    """
    # read in signal data
    atac_data = pd.read_csv(atac_mat_file, sep="\t", index_col=0)
    if is_ggr:
        # for GGR, need to remove d05
        atac_data = atac_data.drop("d05", axis=1)

    # read in clusters and process
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
    run_shell_cmd(plot_cmd)
    
    return
