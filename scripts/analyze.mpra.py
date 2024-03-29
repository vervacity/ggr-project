#!/usr/bin/env python

import os
import re
import gzip

import numpy as np
import pandas as pd
import networkx as nx

from scipy.stats import wilcoxon
from scipy.stats import ranksums
from scipy.stats import ttest_1samp
from scipy.stats import ttest_ind


def trim_and_attach_umi(fastq1, fastq2, fastq_out):
    """set up dependent, here the UMI is the first 15 bp of read1, but we want the
    20 bp after 17 bps in read 2
    """
    with gzip.open(fastq1, "r") as fp1:
        with gzip.open(fastq2, "r") as fp2:
            with gzip.open(fastq_out, "w") as out:
                
                reads_seen = 0
                while True:

                    # progress
                    if reads_seen % 1000000 == 0:
                        print reads_seen

                    # get a read pair
                    try:
                        read_1 = []
                        for i in range(4):
                            read_1.append(fp1.readline().strip())

                        if read_1[0] == "":
                            break

                        read_2 = []
                        for i in range(4):
                            read_2.append(fp2.readline().strip())

                        if read_2[0] == "":
                            break

                    except:
                        break
                    
                    # find the umi (read 1)
                    umi = read_1[1][:15]

                    # trim lines
                    read_2[0] = read_2[0].split()[0] + ":" + umi
                    read_2[1] = read_2[1][17:]
                    read_2[3] = read_2[3][17:]

                    # write out
                    out.write("{}\n".format("\n".join(read_2)))

                    # continue
                    reads_seen += 1
                    
    return


def dedup_and_count(sam_file, dedup_file):
    """go through sorted sam file and get deduped counts
    """
    umis = set()
    barcode = ""
    count = 0
    with open(sam_file, "r") as fp:
        with gzip.open(dedup_file, "w") as out:

            for line in fp:

                # don't read header
                if line.startswith("@"):
                    continue

                # read in important details
                fields = line.split("\t")
                line_barcode = fields[2]
                umi = fields[0].split(":")[-1]
                
                # check if current barcode or not
                if line_barcode != barcode:
                    # save out latest
                    if barcode != "":
                        result = "{}\t{}\n".format(barcode, count)
                        out.write(result)
                    
                    # and reset variables
                    umis = set([umi])
                    barcode = line_barcode
                    count = 1
                else:
                    # check if umi present or not
                    if umi in umis:
                        continue
                    # add to count
                    umis.add(umi)
                    count += 1
    
    return


def fastq_to_counts(fastq_file, ref_fasta_file, prefix):
    """take a fastq file through to deduped counts
    """
    
    # align
    out_sai = "{}.sai".format(prefix)
    align = "bwa aln -t 24 {} {} > {}".format(ref_fasta_file, fastq_file, out_sai)
    if not os.path.isfile(out_sai):
        print align
        os.system(align)

    out_sam = "{}.sam".format(prefix)
    samse = "bwa samse {} {} {} > {}".format(ref_fasta_file, out_sai, fastq_file, out_sam)
    if not os.path.isfile(out_sam):
        print samse
        os.system(samse)

    # filter/sort
    filter_sam = "{}.filt.sam".format(prefix)
    filter_fn = "samtools view -h -F 4 {} | samtools sort -@24 -o {} -".format(
        out_sam, filter_sam)
    if not os.path.isfile(filter_sam):
        print filter_fn
        os.system(filter_fn)

    # get counts
    counts_file = "{}.dedup.counts.txt.gz".format(prefix)
    if not os.path.isfile(counts_file):
        dedup_and_count(filter_sam, counts_file)

    return None


def counts_to_mat(counts_files, out_prefix, metadata_file=None, norm_method="rlog"):
    """read in all files and merge to matrix. produces:
    - count matrix (no metadata)
    - redistributed count  matrix (no metadata)
    - normalized signal matrix (no metadata)
    - normalized signal matrix (metadata)
    """
    mat_data = None
    for key in sorted(counts_files.keys()):
        print key
        count_data = pd.read_csv(counts_files[key], sep="\t", header=None, index_col=0)
        count_data.columns = [key]
        if mat_data is None:
            mat_data = count_data.copy()
        else:
            mat_data = mat_data.merge(
                count_data, how="outer",
                left_index=True, right_index=True)

    # count matrix
    count_mat_file = "{}.counts.raw.mat.txt.gz".format(out_prefix)
    mat_data.to_csv(count_mat_file, sep="\t", compression="gzip", header=True, index=True)

    # normalize: divide by read total (column norm)
    mat_data_norm = mat_data.fillna(0)
    mat_data_norm = mat_data_norm / np.sum(mat_data_norm, axis=0)

    # normalize: multiply by plasmid fract
    plasmid_fracts = mat_data_norm["plasmid"]
    mat_data_norm = mat_data_norm.drop("plasmid", axis=1)
    mat_data_norm = mat_data_norm.multiply(plasmid_fracts.values, axis=0)
    
    # renormalize: make probability (sum to 1) then multiply by total read count
    mat_data_norm = mat_data_norm / np.sum(mat_data_norm, axis=0)
    mat_data_norm = mat_data_norm.multiply(
        np.sum(mat_data.drop("plasmid", axis=1), axis=0), axis=1)
    
    norm_mat_file = "{}.counts.norm.mat.txt.gz".format(out_prefix)
    mat_data_norm.astype(int).to_csv(norm_mat_file, sep="\t", compression="gzip", header=True, index=True)

    # log norm
    norm_signal_mat_file = "{}.signal.mat.txt.gz".format(out_prefix)
    if norm_method == "rlog":
        # also try rlog norm
        rlog_cmd = "~/git/ggr-project/R/normalize.rlog.R {} {}".format(
            norm_mat_file, norm_signal_mat_file)
        print rlog_cmd
        os.system(rlog_cmd)
        signal_mat = pd.read_csv(norm_signal_mat_file, sep="\t")
    else:
        # and make pseudo TPM (log2)
        signal_mat = mat_data_norm / np.sum(mat_data_norm, axis=0)
        signal_mat = signal_mat * 1000000
        signal_mat = np.log2(signal_mat)
        signal_mat[signal_mat < 0] = 0
        signal_mat.to_csv(
            norm_signal_mat_file, sep="\t", compression="gzip", header=True, index=True)
    
    # attach metadata
    if metadata_file is not None:
        metadata = pd.read_csv(metadata_file, header=0, sep="\t", index_col=0)
        signal_mat["unique_id"] = [
            re.sub("\\|.+", "", id_string) for id_string in signal_mat.index.values]
        mat_data_w_metadata = signal_mat.merge(metadata, on="unique_id")
        w_metadata_file = "{}.signal.w_metadata.mat.txt.gz".format(out_prefix)
        mat_data_w_metadata.to_csv(
            w_metadata_file, sep="\t",
            compression="gzip", header=True, index=False)
    
    return None


def _load_data_and_setup(data_file, drop_reps, filter_synergy):
    """load data and set up
    """
    # load data, make grammar col
    data = pd.read_csv(data_file, sep="\t", header=0)
    data["grammar"] = [re.sub("\\.\\d+$", "", name) for name in data["example_id"]]
    
    # if filter, then get just synergy set (remove distant and non-sig)
    data["synergy.max"] = data["synergy.scores.diff.sig.0"].apply(
        lambda x: np.max([int(val) for val in x.split(",")]))
    if filter_synergy:
        data = data[data["synergy.max"] != 0]
    
    # which headers to keep
    signal_headers = [name for name in data.columns if "_b" in name]
    for drop_rep in drop_reps:
        signal_headers = [name for name in signal_headers if drop_rep not in name]

    # get hypotheses
    grammars = list(set(data["grammar"].values.tolist()))
    grammars = sorted([name for name in grammars if "ggr" in name])

    return data, signal_headers, grammars


def _average_data(
        data, signal_headers, metadata_headers, signal_thresh,
        filter_synergy=False,
        days=["d0", "d3", "d6"]):
    """assumes you have signal data and corresponding metadata
    """
    # select out signal data and relevant metadata
    data = data[signal_headers+metadata_headers]
    if filter_synergy:
        data = data[data["synergy.max"] != 0]
    
    # apply signal thresh and remove rows of zeros
    data = data.set_index(metadata_headers)
    data[data < signal_thresh] = 0.0
    data = data[data.max(axis=1) > 0]
    data = data.replace(0, np.NaN) # ignore zeros, was not seen in replicate
    data = data.reset_index()
        
    # average the barcoded fragments in each REP first
    # this is because different reps may have different barcodes
    data = data.groupby(metadata_headers).mean()
    data = data.reset_index()
        
    # then average the reps
    data = data.set_index(metadata_headers)
    for day in days:
        day_headers = [header for header in signal_headers if day in header]
        data["{}_AVG".format(day)] = data[day_headers].mean(axis=1)
    data = data.drop(signal_headers, axis=1)
    data = data.reset_index()
    
    return data


def _normalize_for_trajectory_tests(
        grammar_data_avg, metadata_headers, method="d0_norm"):
    """normalize for trajectory plots/analysis
    """
    # mean center
    grammar_data_traj = grammar_data_avg.set_index(metadata_headers)
    grammar_data_traj[:] = np.subtract(
        grammar_data_traj.values,
        np.expand_dims(np.nanmean(grammar_data_traj.values, axis=1), axis=-1))
    grammar_data_traj = grammar_data_traj.reset_index()

    # melt
    grammar_data_traj = grammar_data_traj.melt(id_vars=metadata_headers)
    grammar_data_traj = grammar_data_traj.dropna(axis=0)
    grammar_data_traj["day"] = grammar_data_traj["variable"].apply(
        lambda x: re.sub("_.+", "", x))
    
    if method == "d0_norm":
        # then adjust to 0 median for each combo
        combos = ["-1", "0,0", "1,0", "0,1", "1,1"]
        for combo in combos:
            d0_median = grammar_data_traj[
                (grammar_data_traj["combos"] == combo) &
                (grammar_data_traj["day"] == "d0")]["value"].median()
            grammar_data_traj.loc[grammar_data_traj["combos"] == combo, "value"] = grammar_data_traj.loc[
                grammar_data_traj["combos"] == combo, "value"] - d0_median
    elif method == "d0_endog_norm":
        # normalize all to 0 median for endogenous
        d0_median = grammar_data_traj[
            (grammar_data_traj["combos"] == "0,0") &
            (grammar_data_traj["day"] == "d0")]["value"].median()
        grammar_data_traj["value"] = grammar_data_traj["value"] - d0_median
    elif method == "d0_double_mut_norm":
        # normalize all to 0 median for double mut
        d0_median = grammar_data_traj[
            (grammar_data_traj["combos"] == "1,1") &
            (grammar_data_traj["day"] == "d0")]["value"].median()
        grammar_data_traj["value"] = grammar_data_traj["value"] - d0_median
    else:
        raise ValueError, "method not implemented!"
        
    return grammar_data_traj


def _normalize_for_mutational_test(
        grammar_data_avg, metadata_headers,
        example_background_norm=True):
    """normalize for mutational plots/analysis
    normalize each example separately
    """
    examples = sorted(list(set(grammar_data_avg["example_id"].values.tolist())))
    grammar_data_norm = []
    for example in examples:
        # get example data
        example_data = grammar_data_avg[grammar_data_avg["example_id"] == example]
        
        # check that all necessary examples exist
        if example_data[example_data["combos"] == "1,1"].shape[0] != 1:
            continue

        # normalize relative to baseline (1,1)
        if example_background_norm:
            tmp_data = example_data.drop("example_id", axis=1)
            background_vals = tmp_data[tmp_data["combos"] == "1,1"].drop(
                ["combos", "synergy.max"], axis=1)
            example_data = example_data.set_index(metadata_headers)
            example_data[:] = np.subtract(example_data.values, background_vals.values)
            example_data = example_data.reset_index()
        else:
            # try mean center then subtract out null median val
            example_data_mean = example_data.set_index(metadata_headers)
            example_data_mean[:] = np.subtract(
                example_data_mean.values,
                np.expand_dims(np.nanmean(example_data_mean.values, axis=0), axis=0))
            example_data = example_data.reset_index()
            
        # append
        grammar_data_norm.append(example_data)

    if len(grammar_data_norm) == 0:
        print "no examples, can't run analysis"
        return None

    # concat
    grammar_data_norm = pd.concat(grammar_data_norm, axis=0)
    
    # and then adjust according to "1,1" mean
    if not example_background_norm:
        grammar_data_norm = grammar_data_norm.drop("index", axis=1)
        background_vals = grammar_data_norm[grammar_data_norm["combos"] == "1,1"].set_index(metadata_headers).mean(axis=0) # mean better
        grammar_data_norm = grammar_data_norm.set_index(metadata_headers)
        grammar_data_norm[:] = grammar_data_norm.values - np.expand_dims(background_vals.values, axis=0)
        grammar_data_norm = grammar_data_norm.reset_index()
    
    return grammar_data_norm


def _check_traj_endogenous_vs_double_mut(
        data_traj,
        pval_thresh=0.10,
        method="ranksums",
        days=["d0", "d3", "d6"]):
    """check for each time point, return true if different in at least one
    """
    passed = False
    for day in days:
        # pivot
        data_day = data_traj[data_traj["day"] == day]
        data_day = data_day.pivot(
            index="example_id", columns="combos")["value"].reset_index()

        # stats
        if method == "ranksums":
            endogenous = data_day["0,0"].dropna()
            double_mut = data_day["1,1"].dropna()
            stat, pval = ranksums(endogenous.values, double_mut.values)
        else:
            # only use those that have 1,1 and 0,0
            data_day = data_day[
                ~pd.isna(data_day["1,1"]) &
                ~pd.isna(data_day["0,0"])]
            endogenous = data_day["0,0"]
            double_mut = data_day["1,1"]
            diff = endogenous - double_mut.mean()
        
            # test
            stat, pval = wilcoxon(diff.values)
            
        # passed
        if pval < pval_thresh:
            passed = True
            
    return passed


def _check_traj_pattern(traj_medians, day_key):
    """make sure maximal expression at right point
    """
    # check if headed in right direction - check which is highest
    traj_match = np.max(traj_medians) == traj_medians[day_key]

    # if not, check again (d3 vs d6)
    if not traj_match:
        if day_key == "d3_AVG":
            check_match = np.max(traj_medians) == traj_medians["d6_AVG"]
            if check_match:
                day_key = "d6_AVG"
        elif day_key == "d6_AVG":
            check_match = np.max(traj_medians) == traj_medians["d3_AVG"]
            if check_match:
                day_key = "d3_AVG"
        else:
            check_match = False
    
    # compare
    traj_match = traj_match or check_match
            
    return traj_match, day_key


def _check_traj_sigdiff(traj_data, threshold=0.05):
    """check for statistical difference betweeen days, using wilcoxon
    """
    traj_data_clean = traj_data.set_index("example_id")
    num_columns = traj_data_clean.shape[1]

    # compare all timepoints to each other
    sig = False
    for i in range(num_columns):
        for j in range(num_columns):
            if i == j:
                continue
            
            timepoint_i = traj_data_clean.values[:,i]
            timepoint_i = timepoint_i[~np.isnan(timepoint_i)]
            timepoint_j = traj_data_clean.values[:,j]
            timepoint_j = timepoint_j[~np.isnan(timepoint_j)]
            
            # run test
            diff = timepoint_i - np.median(timepoint_j)
            stat, pval = wilcoxon(diff)
            if pval < threshold:
                sig = True
    
    return sig


def _check_traj_signull(traj_data, threshold=0.05):
    """check compared to null, whether significantly different
    """
    # pivot
    endog_norm = traj_data[traj_data["combos"]=="0,0"].drop(
        ["combos", "day"], axis=1)
    endog_norm = endog_norm.pivot(
        index="example_id", columns="variable")["value"] #.reset_index()

    null_norm = traj_data[traj_data["combos"]=="1,1"].drop(
        ["combos", "day"], axis=1)
    null_norm = null_norm.pivot(
        index="example_id", columns="variable")["value"] #.reset_index()
    
    # only check days 3,6
    sig = False
    for i in range(1,3):
        #print i
        endog_vals = endog_norm.values[:,i]
        endog_vals = endog_vals[~np.isnan(endog_vals)]
        null_vals = null_norm.values[:,i]
        null_vals = null_vals[~np.isnan(null_vals)]
        
        # test
        diff = null_vals - np.median(endog_vals)
        stat, pval = wilcoxon(diff)
        if pval < threshold:
            sig = True
            
    return sig


def _check_mut_interaction(
        expected, actual,
        pval_thresh=0.10,
        method="wilcox_unpaired"):
    """determine if interaction is synergy, additive, buffer
    """
    if method == "wilcox_paired":
        diff = actual - expected
        stat, pval = wilcoxon(diff)
        diff_sum = np.mean(diff)
    elif method == "wilcox_unpaired":
        diff = actual - expected.mean()
        diff = diff[~np.isnan(diff)]
        stat, pval = wilcoxon(diff)
        diff_sum = np.mean(diff)
    elif method == "ranksums":
        stat, pval = ranksums(actual, expected)
        diff = actual - expected.mean()
        diff = diff[~np.isnan(diff)]
        diff_sum = np.mean(diff)
    elif method == "ttest":
        diff = actual - expected.mean()
        diff = diff[~np.isnan(diff)]
        stat, pval = ttest_1samp(diff, 0, nan_policy="omit")
        diff_sum = np.mean(diff)
        
    # summary info
    if (pval < pval_thresh) and (diff_sum > 0):
        interaction = "SYNERGY"
    elif (pval < pval_thresh) and (diff_sum < 0):
        interaction = "BUFFER"
    else:
        interaction = "ADDITIVE"
        
    return interaction, pval


def _attach_null_result(results, interaction_msg):
    """attach null result
    """
    results["mut.medians"].append(np.array([0,0,0,0,0]))
    results["actual"].append(0)
    results["expected"].append(0)
    results["diff"].append(0)
    results["traj.sig"].append(0)
    results["pval"].append(1)
    results["interaction"].append(interaction_msg)

    return results


def analyze_combinations(
        data_file,
        grammar_dir,
        out_dir,
        signal_thresh=0.5,
        mut_pval_thresh=0.10,
        drop_reps=["b4"],
        filter_synergy=True,
        plot=False):
    """analyze results and do statistical tests
    """
    r_dir = "/users/dskim89/git/ggr-project/R"
    
    # load data and set up
    data, signal_headers, grammars = _load_data_and_setup(
        data_file, drop_reps, filter_synergy)
    metadata_headers = ["example_id", "combos", "synergy.max"]
    
    # results
    results = {
        "grammar": [],
        "day": [],
        "day_empirical": [],
        "pwm1": [],
        "pwm1_clean": [],
        "pwm2": [],
        "pwm2_clean": [],
        "traj.medians": [],
        "traj.sig": [],
        "pval": [],
        "mut.medians": [],
        "actual": [],
        "expected": [],
        "diff": [],
        "interaction": []}
    
    # go through hypotheses
    for grammar in grammars:
        print grammar
        grammar_prefix = "{}/{}".format(out_dir, grammar)
        
        # pull gml, get nodes, sort
        grammar_gml = "{}/{}.gml".format(grammar_dir, grammar)
        grammar_graph = nx.read_gml(grammar_gml)
        nodes = list(grammar_graph.nodes.data("pwmidx"))
        nodes.sort(key=lambda x: x[1])
        
        # 0,1 means the first pwm (index order) is still present
        pwm1 = nodes[0][0]
        pwm2 = nodes[1][0]

        # clean
        clean1 = re.sub("HCLUST-\\d+_", "", pwm1)
        pwm1_clean = re.sub(".UNK.0.A", "", clean1)
        clean2 = re.sub("HCLUST-\\d+_", "", pwm2)
        pwm2_clean = re.sub(".UNK.0.A", "", clean2)
        
        # get traj to day
        traj_to_day = {
            "0": "d0",
            "7": "d0",
            "8-10-11": "d0",
            "9": "d0",
            "12-13": "d3",
            "14": "d3",
            "1": "d3",
            "2": "d6",
            "3-4": "d6",
            "5": "d6"}
        traj_val = grammar.split(".grammar")[0].split("_x")[0].split("ggr.")[1].split("TRAJ_LABELS-")[1]
        try:
            day = traj_to_day[traj_val]
        except:
            print "ADD IN INTERSECTION {}".format(grammar)
            continue
        day_key = "{}_AVG".format(day)
        
        # extract relevant data
        grammar_data = data[data["grammar"] == grammar]

        # ignore triplets for now
        if "0,0" not in grammar_data["combos"].values.tolist():
            continue

        # start saving results
        results["grammar"].append(grammar)
        results["day"].append(day)
        results["pwm1"].append(pwm1)
        results["pwm1_clean"].append(pwm1_clean)
        results["pwm2"].append(pwm2)
        results["pwm2_clean"].append(pwm2_clean)
        
        # ----------------------
        # trajectory hypothesis
        # ----------------------
        # normalize to each combos' day 0 and plot
        # note that since the results are now relative to day 0,
        # it's important to normalize each combo to its day 0
        # and NOT normalize all to endog's day 0
        #grammar_data_avg = _average_data(
        #    grammar_data, signal_headers, metadata_headers, signal_thresh, filter_synergy=True)
        if "ggr.TRAJ_LABELS-1_x_TRAJ_LABELS-14_x_TRAJ_LABELS-2.42_x_21_x_58" in grammar:
            # tODO only filter synergy for klf/znf?
            grammar_data_avg = _average_data(
                grammar_data, signal_headers, metadata_headers, signal_thresh, filter_synergy=True)
        else:
            grammar_data_avg = _average_data(
                grammar_data, signal_headers, metadata_headers, signal_thresh, filter_synergy=filter_synergy)
        grammar_data_traj = _normalize_for_trajectory_tests(
            grammar_data_avg, metadata_headers, method="d0_norm")
        
        if plot:
            # drop 1 max and 1 min value (for extreme outliers)
            plot_data_file = "{}.traj.data.txt.gz".format(grammar_prefix)
            plot_data = grammar_data_traj.copy()
            # TODO - think about just adjusting plot area rather than removing this
            plot_data = plot_data[plot_data["value"] != plot_data["value"].max()]
            plot_data = plot_data[plot_data["value"] != plot_data["value"].min()]
            # do again (drop 2 maxes, 2 mins)
            plot_data = plot_data[plot_data["value"] != plot_data["value"].max()]
            plot_data = plot_data[plot_data["value"] != plot_data["value"].min()]

            # adjust combo names
            plot_data["pwms"] = "scr-CTL"
            plot_data.loc[plot_data["combos"] == "0,0", "pwms"] = "{0},{1}".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "0,1", "pwms"] = "{0},{1}-scr".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "1,0", "pwms"] = "{0}-scr,{1}".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "1,1", "pwms"] = "double scramble"
            plot_data.loc[plot_data["combos"] == "-1", "pwms"] = "scr-CTL"
            plot_data.to_csv(
                plot_data_file, sep="\t", compression="gzip", header=True, index=False)
            plot_cmd = "{}/plot.mpra.hypothesis.traj.R {} {}.traj".format(
                r_dir, plot_data_file, grammar_prefix)
            print plot_cmd
            os.system(plot_cmd)
            os.system("rm {}".format(plot_data_file))
            
        # save out medians for endogenous seq
        endog_norm = grammar_data_traj[grammar_data_traj["combos"]=="0,0"].drop(
            ["combos", "day"], axis=1)
        endog_norm = endog_norm.pivot(
            index="example_id", columns="variable")["value"].reset_index()
        traj_medians = endog_norm.median()
        results["traj.medians"].append(traj_medians)

        # check traj pattern
        passed_traj_sigdiff = _check_traj_sigdiff(endog_norm)
        # TODO check here
        passed_traj_pattern, day_key = _check_traj_pattern(traj_medians, day_key)
        if "ggr.TRAJ_LABELS-1_x_TRAJ_LABELS-14_x_TRAJ_LABELS-2.42_x_21_x_58" in grammar:
            day_key = "d6_AVG"
        results["day_empirical"].append(day_key)
        if not (passed_traj_pattern and passed_traj_sigdiff):
            results = _attach_null_result(results, "FAILED.TRAJ_PATTERN")
            continue

        # run test here, results saved out later
        #passed_traj_vs_double_mut = _check_traj_endogenous_vs_double_mut(
        #    grammar_data_traj, days=["d0", "d3", "d6"])

        # compare to null, is there a sig diff from null?
        # ^ test is above
        passed_traj_vs_double_mut = _check_traj_signull(grammar_data_traj)
        
        # compare to ATAC, H3K27ac
        if True:
            # plot the atac/h3k27ac
            keys = ["ATAC_SIGNALS.NORM", "H3K27ac_SIGNALS.NORM"]
            for key in keys:
                keep_keys = ["example_id"] + [key]
                ids = endog_norm["example_id"].values

                # parse epigenome data
                epigenome_data = data[data["example_id"].isin(ids)][keep_keys]
                epigenome_data = epigenome_data.drop_duplicates().set_index("example_id")
                epigenome_data = epigenome_data[key].str.split(",", expand=True).astype(float)
                if "ATAC" in key:
                    epigenome_data = epigenome_data[[0,6,12]]
                epigenome_data.columns = ["d0", "d3", "d6"]
                epigenome_medians = epigenome_data.median()

                # normalize
                epigenome_data = epigenome_data.sub(epigenome_data.mean(axis=1), axis=0)
                epigenome_data = epigenome_data.melt()
                epigenome_data["value"] = epigenome_data["value"] - epigenome_data[epigenome_data["variable"] == "d0"]["value"].median()
                
                # compare epigenome data to expression data
                plot_data_file = "{}.traj.epigenome_data.txt.gz".format(grammar_prefix)
                epigenome_data.to_csv(plot_data_file, sep="\t", compression="gzip", header=True, index=False)
                plot_cmd = "{}/plot.mpra.vs_epigenome.traj.R {} {}.{}".format(
                    r_dir, plot_data_file, grammar_prefix, key)
                print plot_cmd
                os.system(plot_cmd)
                os.system("rm {}".format(plot_data_file))
        
        # ----------------------
        # mutational hypothesis
        # ----------------------
        # normalize and only keep those that have 1,1 (baseline val)
        grammar_data_avg = _average_data(
            grammar_data, signal_headers, metadata_headers, signal_thresh)
        grammar_data_mut = _normalize_for_mutational_test(
            grammar_data_avg, metadata_headers)
        
        # don't continue if no examples
        if grammar_data_mut is None:
            results = _attach_null_result(results, "FAILED.CANT_TEST")
            continue
        
        # plot
        if plot:
            # remove 1 max and 1 min val (extreme outliers)
            plot_data_file = "{}.mut.agg_data.txt.gz".format(grammar_prefix)
            plot_data = grammar_data_mut.copy()
            day_keys = ["d0_AVG", "d3_AVG", "d6_AVG"]
            max_val = plot_data[day_keys].max().max()
            min_val = plot_data[day_keys].min().min()
            for day in day_keys: # NOTE: day_key is an important variable 
                plot_data.loc[plot_data[day] == max_val, day] = np.nan
                plot_data.loc[plot_data[day] == min_val, day] = np.nan
            # adjust combo names
            plot_data["pwms"] = "scr-CTL"
            plot_data.loc[plot_data["combos"] == "0,0", "pwms"] = "{0},{1}".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "0,1", "pwms"] = "{0},scr".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "1,0", "pwms"] = "scr,{1}".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "1,1", "pwms"] = "double scramble".format(pwm1_clean,pwm2_clean)
            plot_data.loc[plot_data["combos"] == "-1", "pwms"] = "scr-CTL"
            plot_data.to_csv(
                plot_data_file, sep="\t", compression="gzip", header=True, index=False)
            plot_cmd = "{0}/plot.mpra.hypothesis.mut.R {1} {2} {3}.mut.{2}".format(
                r_dir, plot_data_file, day_key, grammar_prefix)
            print plot_cmd
            os.system(plot_cmd)
            os.system("rm {}".format(plot_data_file))

        # testing - pivot table
        check_data = grammar_data_mut.set_index(
            metadata_headers)[[day_key]].reset_index()
        check_data = check_data.pivot(index="example_id", columns="combos")[day_key].reset_index()
        check_data = check_data.set_index("example_id")
        
        # don't continue if endogenous is not positive
        passed_positive_over_double_mut = check_data["0,0"].mean() >= 0 # saves more than median
        #passed_positive_over_double_mut = (
        #    check_data["0,0"].mean() >= 0) and (check_data["0,0"].median() >= 0)
        if not passed_positive_over_double_mut:
            print "background higher than endogenous, not well controlled, continue"
            results = _attach_null_result(results, "FAILED.NEGATIVE_ENDOG")
            continue

        # mut medians
        mut_medians = check_data.median()
        results["mut.medians"].append(mut_medians.values)

        # determine if expected and actual differ
        expected = check_data["0,1"] + check_data["1,0"]
        actual = check_data["0,0"]
        interaction, pval = _check_mut_interaction(
            expected, actual,
            pval_thresh=mut_pval_thresh)

        # save out
        results["actual"].append(actual.mean())
        results["expected"].append(expected.mean())
        results["diff"].append((actual - expected).mean())
        results["interaction"].append(interaction)

        # separately consider traj sig        
        if passed_traj_vs_double_mut:
            results["traj.sig"].append(1)
        else:
            results["traj.sig"].append(0)
        results["pval"].append(pval)
        
    # create summary
    results_adj = []
    days = ["d0", "d3", "d6"]
    for key in sorted(results.keys()):
        print key
        data = np.stack(results[key], axis=0)
        print key, data.shape
        if "traj.medians" in key:
            header = ["{}.{}".format(key, day) for day in days]
            data = pd.DataFrame(data=data, columns=header)
        elif "mut.medians" in key:
            header = ["-1", "0,0", "0,1", "1,0", "1,1"]
            data = pd.DataFrame(data=data, columns=header)
        else:
            data = pd.DataFrame(data=data, columns=[key])
        results_adj.append(data)
    results = pd.concat(results_adj, axis=1)

    # some extra calcs
    traj_headers = ["traj.medians.{}".format(day) for day in days]
    results["max_diff"] = results[traj_headers].abs().max(axis=1)
    results = results.sort_values(["interaction", "max_diff"])
    
    # save out
    summary_file = "{}/summary.txt.gz".format(out_dir)
    results.to_csv(summary_file, sep="\t", compression="gzip", index=False, header=True)
    
    return None


def _plot_epigenome(data_file, plot_prefix, keep_rules=[], key="ATAC_SIGNALS.NORM"):
    """pull in normalized file and plot
    """
    # load
    data = pd.read_csv(data_file, sep="\t")
    data["rule"] = [re.sub("\\.\\d+$", "", name) for name in data["example_id"]]
    
    # which headers to keep
    signal_headers = [name for name in data.columns if "_b" in name]
    metadata_headers = [key, "example_metadata", "combos", "rule"]

    # filt
    data_filt = data[signal_headers + metadata_headers]
    data_filt = data_filt[data_filt[key] != "0"]
    data_filt = data_filt[data_filt["combos"] == "0,0"]
    data_filt = data_filt.drop("combos", axis=1)
    data_filt = data_filt[data_filt[signal_headers].mean(axis=1) != 0]
    keep_rules = []
    if len(keep_rules) > 0:
        data_filt = data_filt[data_filt["rule"].isin(keep_rules)]
    
    # get epigenome data
    epigenome_data = data_filt[key].str.split(",", expand=True).astype(float)
    if "ATAC" in key:
        epigenome_data = epigenome_data[[0,6,12]]
    if "logit" in key:
        epigenome_data = epigenome_data[[0,6,12]]
    
    # plot per timepoint
    days = ["d0", "d3", "d6"]
    for day_idx in range(len(days)):
        day = days[day_idx]
        # get mpra data
        mpra_day_headers = [name for name in data_filt.columns if day in name]
        mpra_vals = data_filt[mpra_day_headers].median(axis=1) # median?

        # get epigenome vals
        epigenome_vals = epigenome_data.values[:,day_idx]

        # make df
        plot_data = pd.DataFrame({
            "example_ids": data_filt["example_metadata"].values,
            "mpra": mpra_vals.values,
            "epigenome": epigenome_vals
        })

        # reduce
        plot_data = plot_data.groupby("example_ids").median() # mean?

        # save out and plot
        plot_data_file = "{}.{}.data.txt".format(plot_prefix, day)
        plot_file = "{}.{}.pdf".format(plot_prefix, day)
        plot_data.to_csv(plot_data_file, sep="\t", header=True)
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.mpra.compare.R {} {}".format(
            plot_data_file, plot_file)
        print plot_cmd
        os.system(plot_cmd)
        
    return


def _plot_epigenome_both(atac_prefix, h3k27ac_prefix, plot_prefix):
    """try a sum of both
    """
    days = ["d0", "d3", "d6"]
    for day_idx in range(len(days)):
        day = days[day_idx]

        atac_file = "{}.{}.data.txt".format(atac_prefix, day)
        h3k27ac_file = "{}.{}.data.txt".format(h3k27ac_prefix, day)

        # read in
        atac_data = pd.read_csv(atac_file, sep="\t").set_index("example_ids")
        atac_data.columns = ["ATAC", "mpra"]
        h3k27ac_data = pd.read_csv(h3k27ac_file, sep="\t").set_index("example_ids")
        h3k27ac_data.columns = ["H3K27ac", "mpra"]
        h3k27ac_data = h3k27ac_data.drop("mpra", axis=1)
        
        # merge
        data = atac_data.merge(
            h3k27ac_data, how="outer", left_index=True, right_index=True)
        data = data.fillna(0)
        data["epigenome"] = data["ATAC"] + data["H3K27ac"]
        #print data
        #quit()

        # save out and plot
        plot_data_file = "{}.{}.data.txt".format(plot_prefix, day)
        plot_file = "{}.{}.pdf".format(plot_prefix, day)
        data.to_csv(plot_data_file, sep="\t", header=True)
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.mpra.compare.R {} {}".format(
            plot_data_file, plot_file)
        print plot_cmd
        os.system(plot_cmd)

    return


def main():
    """analyze mpra data
    """
    # inputs
    PLASMID_DIR = "./raw/plasmid"
    DATA_DIR = "./raw"
    OUT_DIR = "./results"
    ref_fasta_file = "./raw/annotations/sequences.fa"
    metadata_file = "./raw/annotations/mpra.seqs.run-0.txt"
    grammar_dir = "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/dmim.shuffle.OLD/grammars.annotated.manual_filt.merged.final"
    
    os.system("mkdir -p {}".format(OUT_DIR))
    
    # make sure to load bwa 0.7.10

    # generate index if needed
    if not os.path.isfile("{}.bwt".format(ref_fasta_file)):
        bwa_index_cmd = "bwa index {}".format(ref_fasta_file)
        print bwa_index_cmd
        os.system(bwa_index_cmd)
        
    # ===========
    # GET DNA PLASMID COUNTS
    # ===========
    
    # analyze the plasmid library (to get normalization fractions)
    # need to look at r1 and r2 to get the UMIs to remove PCR duplicates
    plasmid_fastq1 = "{}/DK-C2-plasmidlib_S1_L001_R1_001.fastq.gz".format(PLASMID_DIR)
    plasmid_fastq2 = "{}/DK-C2-plasmidlib_S1_L001_R2_001.fastq.gz".format(PLASMID_DIR)
    plasmid_trim_fastq = "{}/plasmid.trim.fastq.gz".format(PLASMID_DIR)

    # trim and get umis
    if not os.path.isfile(plasmid_trim_fastq):
        trim_and_attach_umi(plasmid_fastq1, plasmid_fastq2, plasmid_trim_fastq)

    # get counts
    plasmid_counts = "{}/plasmid.dedup.counts.txt.gz".format(OUT_DIR)
    if not os.path.isfile(plasmid_counts):
        fastq_to_counts(plasmid_trim_fastq, ref_fasta_file, "{}/plasmid".format(OUT_DIR))
        
    # ===========
    # GET EXPRESSION COUNTS
    # ===========
    
    # args
    mixed_fastq_file = "{}/non_demux/pooled_R1_001.fastq.gz".format(DATA_DIR)
    demux_fastqs = {
        "d0_b1": "{}/demux/raw_data/S1/S1_DK-1-D0_S1_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d3_b1": "{}/demux/raw_data/S2/S2_DK-1-D3_S2_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d6_b1": "{}/demux/raw_data/S3/S3_DK-1-D6_S3_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),

        "d0_b2": "{}/demux/raw_data/S4/S4_DK-2-D0_S4_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d3_b2": "{}/demux/raw_data/S5/S5_DK-2-D3_S5_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d6_b2": "{}/demux/raw_data/S6/S6_DK-2-D6_S6_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),

        "d0_b3": "{}/demux/raw_data/S7/S7_DK-3-D0_S7_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d3_b3": "{}/demux/raw_data/S8/S8_DK-3-D3_S8_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d6_b3": "{}/demux/raw_data/S9/S9_DK-3-D6_S9_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),

        "d0_b4": "{}/demux/raw_data/S10/S10_DK-4-D0_S10_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d3_b4": "{}/demux/raw_data/S11/S11_DK-4-D3_S11_HGJLVDRXX_L1.fq.gz".format(DATA_DIR),
        "d6_b4": "{}/demux/raw_data/S12/S12_DK-4-D6_S12_HGJLVDRXX_L1.fq.gz".format(DATA_DIR)}
    
    # chop the fastq file AND demultiplex AND attach UMIs to ID
    print "demultiplex, trim, attach UMIs to ID"
    trimmed_fastqs = {}
    for key in sorted(demux_fastqs.keys()):
        trimmed_fastq = "{}/{}.trim.fastq.gz".format(OUT_DIR, key)
        trimmed_fastqs[key] = trimmed_fastq
        
    reads_seen = 0
    skipped = 0
    if not os.path.isfile(trimmed_fastq):

        # load in demultiplexed IDs
        print "loading demultiplexed IDs"
        read_ids = {}
        for key in sorted(demux_fastqs.keys()):
            fastq_file = demux_fastqs[key]
            print fastq_file
            line_num = 0
            with gzip.open(fastq_file, "r") as fp:
                for line in fp:
                    if line_num % 4000000 == 0:
                        print "reads:", line_num / 4
                    if line_num % 4 == 0:
                        read_id = line.strip().split()[0]
                        read_ids[read_id] = key
                    line_num += 1

        # read in non demultiplexed file and split out
        with gzip.open(mixed_fastq_file, "r") as fp:
            while True:
                if reads_seen % 1000000 == 0:
                    print reads_seen
                
                # get a read (4 lines)
                read = []
                for i in range(4):
                    read.append(fp.readline().strip())

                    if read[0] == "":
                        break
                    
                if read[0] == "":
                    break
                    
                # keep id
                read_id = read[0].split()[0]
                
                # attach UMI to read name
                new_id = read[0].split()[0] + ":" + read[0][-12:]
                read[0] = new_id
                
                # trim lines
                read[1] = read[1].strip()[0:20]
                read[3] = read[3].strip()[0:20]

                try:
                    fastq_key = read_ids[read_id]
                    with gzip.open(trimmed_fastqs[fastq_key], "a") as out:
                        out.write("{}\n".format("\n".join(read)))
                except:
                    skipped += 1
                    pass

                reads_seen += 1

        print "reads seen", reads_seen
        print "skipped", skipped

    # get counts
    print "fastq -> counts"
    counts_files = {"plasmid": plasmid_counts}
    for key in sorted(trimmed_fastqs.keys()):
        out_prefix = "{}/{}".format(OUT_DIR, key)
        fastq_to_counts(trimmed_fastqs[key], ref_fasta_file, out_prefix)
        counts_files[key] = "{}.dedup.counts.txt.gz".format(out_prefix)

    # and then merge all counts into a matrix file
    out_prefix = "{}/ggr".format(OUT_DIR)
    count_mat_file = "{}.counts.raw.mat.txt.gz".format(out_prefix)
    signal_mat_file = "{}.signal.mat.txt.gz".format(out_prefix)
    if not os.path.isfile(signal_mat_file):
        counts_to_mat(counts_files, out_prefix, metadata_file=metadata_file)
        
    # run QC
    qc_cmd = "Rscript ~/git/ggr-project/R/qc.mpra.R {} {}".format(
        count_mat_file, signal_mat_file)
    print qc_cmd
    #os.system(qc_cmd)
    
    # signal file w metadata
    signal_mat_file = "{}.signal.w_metadata.mat.txt.gz".format(out_prefix)

    # analyze
    results_dir = "{}/combinatorial_rules".format(OUT_DIR)
    os.system("mkdir -p {}".format(results_dir))
    if False:
        analyze_combinations(
            signal_mat_file,
            grammar_dir,
            results_dir,
            signal_thresh=0,
            drop_reps=[],
            filter_synergy=False,
            plot=True)
    
    # plot expected vs actual
    summary_file = "{}/summary.txt.gz".format(results_dir)
    summary_plot_file = "{}/summary.expected_v_actual.pdf".format(results_dir)
    plot_cmd = "Rscript ~/git/ggr-project/R/plot.mpra.summary.R {} {}".format(
        summary_file, summary_plot_file)
    print plot_cmd
    #os.system(plot_cmd)

    # plot rule map
    map_plot_file = "{}/summary.rule_map.pdf".format(results_dir)
    plot_cmd = "Rscript ~/git/ggr-project/R/plot.mpra.rule_map.R {} {}".format(
        summary_file, map_plot_file)
    print plot_cmd
    #os.system(plot_cmd)

    if True:
        results_data = pd.read_csv(summary_file, sep="\t")
        results_data = results_data[results_data["actual"] != 0]
        keep_rules = results_data["grammar"].values.tolist()

        # plot mpra vs ATAC
        atac_prefix = "{}/compare.mpra_vs_atac".format(results_dir)
        _plot_epigenome(signal_mat_file, atac_prefix,
                        keep_rules=keep_rules, key="ATAC_SIGNALS.NORM")

        # plot mpra vs H3K27ac
        h3k27ac_prefix = "{}/compare.mpra_vs_h3k27ac".format(results_dir)
        _plot_epigenome(signal_mat_file, h3k27ac_prefix,
                        keep_rules=keep_rules, key="H3K27ac_SIGNALS.NORM")

        # plot mpra vs logit
        logit_prefix = "{}/compare.mpra_vs_logit".format(results_dir)
        _plot_epigenome(signal_mat_file, logit_prefix,
                        keep_rules=keep_rules, key="logits.norm")

        
        # try plots considering both
        #plot_prefix = "{}/compare.mpra_v_epigenome".format(results_dir)
        #_plot_epigenome_both(atac_prefix, h3k27ac_prefix, plot_prefix)
    
    return

main()
