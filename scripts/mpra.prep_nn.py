
import os

import numpy as np
import pandas as pd

from numpy.random import RandomState

_BASES = ["A", "C", "G", "T"]

def generate_random_seq(rand_seed, length):
    """generate random seq
    """
    rand_state = RandomState(rand_seed)
    rand_seed += 1
    random_seq = rand_state.choice(
        _BASES, size=length)
    random_seq = "".join(random_seq)
        
    return random_seq, rand_seed


def main():
    """set up files for running predictions through tronn
    """
    DATA_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2019-10-22.results"
    OUT_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/validation/mpra.2020-06-29.fit_nn"

    # copy over the results file
    results_file = "{}/results/ggr.signal.w_metadata.mat.txt.gz".format(DATA_DIR)
    os.system("rsync -avz --progress {} {}/".format(results_file, OUT_DIR))

    # reduce the results file
    results = pd.read_csv(results_file, sep="\t")
    vals_headers = [header for header in results.columns if "_b" in header] + ["example_combo_id"]
    results_vals = results[vals_headers].set_index("example_combo_id")
    results_vals[results_vals == 0] = np.nan
    
    # average barcodes first
    results_vals_avg = results_vals.reset_index().groupby("example_combo_id").mean()

    # then average reps
    days = ["d0", "d3", "d6"]
    rep_avgs = {}
    for day in days:
        day_headers = [header for header in results_vals_avg.columns if day in header]
        day_results = results_vals_avg[day_headers]
        rep_avgs[day] = day_results.mean(axis=1).values
    rep_avgs = pd.DataFrame(rep_avgs)
    rep_avgs.index = results_vals_avg.index

    # clean up
    rep_avgs = rep_avgs[~np.all(np.isnan(rep_avgs), axis=1)]
    rep_avgs = rep_avgs.fillna(0)

    # now need to attach relevant metadata
    metadata_headers = [
        "sequence.nn",
        #"sequence.mpra",
        "ATAC_SIGNALS.NORM",
        "H3K27ac_SIGNALS.NORM",
        "example_metadata",
        "example_combo_id"]
    results_metadata = results[metadata_headers]
    results_metadata = results_metadata.drop_duplicates()
    results_metadata = results_metadata.set_index("example_combo_id")

    # merge
    rep_avgs = rep_avgs.merge(
        results_metadata, how="left", left_index=True, right_index=True).reset_index()

    # adjust combo id
    rep_avgs["example_combo_id"] = rep_avgs["example_combo_id"].str.replace("-", "_")
    
    # save this file out
    summ_file = "{}/ggr.mpra_signal.summarized.txt.gz".format(OUT_DIR)
    rep_avgs.to_csv(summ_file, sep="\t", header=True, index=False, compression="gzip")

    # use this file to set up new fasta
    final_len = 1000
    mpra_seq_len = 160
    left_extend = (1000 - mpra_seq_len) / 2
    left_random_seq, _ = generate_random_seq(100, left_extend) # seed 0
    right_extend = (1000 - mpra_seq_len) / 2
    right_random_seq, _ = generate_random_seq(50, right_extend) # seed 1

    new_fasta_file = "{}/sequences_adjusted.fa".format(OUT_DIR)
    with open(new_fasta_file, "w") as out:
        for line_idx in range(rep_avgs.shape[0]):
            example = rep_avgs.iloc[line_idx]
            
            # header
            out.write(">{}\n".format(example["example_combo_id"]))

            # sequence
            sequence = str(example["sequence.nn"])
            if len(sequence) != mpra_seq_len:
                print "problem", example["example_combo_id"], sequence
                print example
                quit()
            
            sequence = left_random_seq + sequence + right_random_seq + "\n"
            out.write(sequence)

    return

main()
