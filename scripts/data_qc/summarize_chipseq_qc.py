
import gzip
import glob

import numpy as np
import pandas as pd

from scipy.stats import pearsonr
from scipy.stats import spearmanr


def get_num_lines_gz(filename):
    num_lines = 0
    with gzip.open(filename, "r") as fp:
        for line in fp:
            num_lines += 1
    
    return num_lines


def main():
    """get stats from PAS-seq
    
    - num reads per file
    - gene quant level spearman correlations

    """
    
    # files
    DATA_DIR = "/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones"
    
    # params
    marks = ["H3K27ac", "H3K4me1", "H3K27me3", "CTCF"]
    days = np.arange(0, 7, 3)
    days = ["d{}".format(day).replace(".", "") for day in days]
    reps = ["1", "2"]
    
    # results 
    results = {}
    results["mark_or_tf"] = [] 
    results["timepoint"] = []
    results["replicate"] = []
    #results["num_input_reads"] = []
    results["num_nodup_reads"] = []
    results["NRF"] = []
    results["PBC1"] = []
    results["PBC2"] = []
    results["num_macs2_peaks"] = []
    results["num_overlap_peaks"] = []
    results["num_idr_peaks"] = []

    for mark in marks:
        print mark
        for day in days:
            for rep in reps:
            
                # timepoint, rep
                results["mark_or_tf"].append(mark)
                results["timepoint"].append(day)
                results["replicate"].append(rep)

                # nodup reads
                nodup_log = glob.glob(
                    "{}/*{}*{}*/qc/rep{}/*nodup.flagstat.qc".format(
                        DATA_DIR, day, mark, rep))[0]
                with open(nodup_log, "r") as fp:
                    for line in fp:
                        if "in total" in line:
                            num_nodup_reads = line.split("+")[0].strip()
                            results["num_nodup_reads"].append(num_nodup_reads)

                # NRF/PBC1/PBC2
                lib_log = glob.glob(
                    "{}/*{}*{}*/qc/rep{}/*nodup.pbc.qc".format(
                        DATA_DIR, day, mark, rep))[0]
                with open(lib_log, "r") as fp:
                    # cols 5,6,7 is NRF/PBC1/PBC2
                    for line in fp:
                        fields = line.strip().split()
                    results["NRF"].append(fields[4])
                    results["PBC1"].append(fields[5])
                    results["PBC2"].append(fields[6])

                # peak files
                macs2_peaks = glob.glob(
                    "{}/*{}*{}*/peak/macs2/rep{}/*narrowPeak.gz".format(
                        DATA_DIR, day, mark, rep))[0]
                num_macs2 = get_num_lines_gz(macs2_peaks)
                results["num_macs2_peaks"].append(num_macs2)

                if "CTCF" in mark:
                    idr_peaks = glob.glob(
                        "{}/*{}*{}*/peak/idr/true_reps/rep1-rep2/*filt.narrowPeak.gz".format(
                            DATA_DIR, day, mark))[0]
                    num_idr = get_num_lines_gz(idr_peaks)
                    results["num_idr_peaks"].append(num_idr)
                    results["num_overlap_peaks"].append("NA")
                else:
                    results["num_idr_peaks"].append("NA")
                    overlap_peaks = glob.glob(
                        "{}/*{}*{}*/peak/macs2/overlap/*filt.narrowPeak.gz".format(
                            DATA_DIR, day, mark, rep))[0]
                    num_overlap = get_num_lines_gz(overlap_peaks)
                    results["num_overlap_peaks"].append(num_overlap)
                
            
    # dataframe
    results = pd.DataFrame(results)
    ordered_headers = [
        "mark_or_tf",
        "timepoint",
        "replicate",
        #"num_input_reads",
        "num_nodup_reads",
        "NRF",
        "PBC1",
        "PBC2",
        "num_macs2_peaks",
        "num_overlap_peaks",
        "num_idr_peaks"]
    results = results[ordered_headers]
    out_file = "ggr.ChIP-seq.QC.summary.txt"
    results.to_csv(out_file, sep="\t", header=True, index=False)
    
    
    return

main()
