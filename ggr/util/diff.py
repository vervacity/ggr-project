"""Functions to help with downstream differential analyses
"""

import os
import pandas as pd

#debug tool
import numpy as np


def join_diff_region_lists_to_mat(
        master_list,
        sig_up_down_pairs,
        out_file,
        diff_name_col=3):
    """Takes pairs of (sig up, sig down) files (zip them ordered together)
    and saves them out to out_file. includes non changed ids from master list
    """
    assert ".gz" in out_file

    master = pd.read_table(master_list, header=None)
    master.columns = ["chrom", "start", "stop"]
    master["id"] = master["chrom"] + ":" + master["start"].map(str) + "-"  + master["stop"].map(str)
    del master["chrom"]
    del master["start"]
    del master["stop"]

    for sig_up_file, sig_down_file in sig_up_down_pairs:

        # sig up file
        up_col_name = os.path.basename(sig_up_file).split(".")[diff_name_col]
        sig_up = pd.read_table(sig_up_file, header=None)
        sig_up.columns = ["id"]
        sig_up[up_col_name] = [1 for i in range(sig_up.shape[0])]

        # merge in and replace NaN with zeros
        master = master.merge(sig_up, how="left", on="id")
        master = master.fillna(0)

        # sig down file
        down_col_name = os.path.basename(sig_down_file).split(".")[diff_name_col]
        sig_down = pd.read_table(sig_down_file, header=None)
        sig_down.columns = ["id"]
        sig_down[down_col_name] = [-1 for i in range(sig_down.shape[0])]

        # merge in and replace NaN with zeros
        master = master.merge(sig_down, how="left", on="id")
        master = master.fillna(0)

        # now pair them up
        final_name = up_col_name.split("_sig")[0]
        master[final_name] = master[up_col_name] + master[down_col_name]
        del master[up_col_name]
        del master[down_col_name]

    master.index = master["id"]
    del master["id"]

    # set up cluster nums
    master["cluster"] = [0 for i in range(master.shape[0])]
    num_diff_cols = master.shape[1] - 1
    for i in range(num_diff_cols):
        master["cluster"] = master["cluster"] + (3**i) * (-1 * master.iloc[:,num_diff_cols - i - 1] + 1)
    master["cluster"] = master["cluster"].astype(int)

    # save
    master.to_csv(out_file, sep='\t', compression="gzip")
    
    return
