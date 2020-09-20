
import glob

import numpy as np
import pandas as pd

from scipy.stats import pearsonr
from scipy.stats import spearmanr


def main():
    """get stats from PAS-seq
    
    - num reads per file
    - gene quant level spearman correlations

    """
    
    # files
    DATA_DIR = "/mnt/lab_data/kundaje/projects/skin/data/bds"
    ALIGN_DIR = "{}/processed.rna.2016-07-25.passeq_alignments".format(DATA_DIR)
    QUANT_DIR = "{}/processed.rna.2016-07-30.passeq_quants".format(DATA_DIR)
    rsem_mat_file = "{}/matrices/primary_keratinocyte.GGR.Stanford_Khavari.PAS-seq.Quants_RSEM.counts".format(QUANT_DIR)
    rsem_norm_mat_file = ""
    
    # load rsem mat file
    rsem_mat = pd.read_csv(rsem_mat_file, sep="\t", index_col=0)
    
    # params
    days = np.arange(0, 6.5, 0.5)
    days = ["{0:3.1f}".format(day) for day in days]
    days = ["d{}".format(day).replace(".", "") for day in days]
    reps = ["b1", "b2"]

    # results 
    results = {}
    results["timepoint"] = []
    results["replicate"] = []
    results["num_input_reads"] = []
    results["num_unique_mappers"] = []
    results["num_multi_mappers"] = []
    results["rsem_total_counts"] = []
    results["rep_corr_pearsonr"] = []
    results["rep_corr_spearmanr"] = []
    
    # collect for each day
    for day in days:
        for rep in reps:
            
            # timepoint, rep
            results["timepoint"].append(day)
            results["replicate"].append(rep)

            # align log
            rep_align_log = glob.glob(
                "{}/*{}*{}*/*Log.final.out".format(ALIGN_DIR, day, rep))[0]
            with open(rep_align_log, "r") as fp:
                for line in fp:
                    if "Number of input reads" in line:
                        num_input_reads = line.split("|")[1].strip()
                        results["num_input_reads"].append(num_input_reads)
                    elif "Uniquely mapped reads number" in line:
                        num_unique_mappers = line.split("|")[1].strip()
                        results["num_unique_mappers"].append(num_unique_mappers)
                    elif "Number of reads mapped to multiple loci" in line:
                        num_multi_mappers = line.split("|")[1].strip()
                        results["num_multi_mappers"].append(num_multi_mappers)

            # rsem counts - should be about same as total alignments, each read gets distributed
            rsem_header = "{}_{}".format(day, rep)
            results["rsem_total_counts"].append(rsem_mat[rsem_header].sum().astype(int))

            # rsem corr
            b1_header = "{}_b1".format(day)
            b2_header = "{}_b2".format(day)
            pearson_r = pearsonr(rsem_mat[b1_header].values, rsem_mat[b2_header].values)[0]
            results["rep_corr_pearsonr"].append(pearson_r)
            spearman_r = spearmanr(rsem_mat[b1_header].values, rsem_mat[b2_header].values)[0]
            results["rep_corr_spearmanr"].append(spearman_r)
            
    # dataframe
    results = pd.DataFrame(results)
    ordered_headers = [
        "timepoint",
        "replicate",
        "num_input_reads",
        "num_unique_mappers",
        "num_multi_mappers",
        "rsem_total_counts",
        "rep_corr_pearsonr",
        "rep_corr_spearmanr"]
    results = results[ordered_headers]
    out_file = "ggr.PAS-seq.QC.summary.txt"
    results.to_csv(out_file, sep="\t", header=True, index=False)
    
    
    return

main()
