"""Contains functions for filtering
"""

import math
import pandas as pd


def remove_media_timepoints(in_mat_file, args, out_mat_file):
    """Takes an input matrix and removes specific timepoints
    """
    media_timepoints = args.misc["media_timepoints"]
    data = pd.read_table(in_mat_file)
    keep_columns = [column for column in data.columns
                    if column.split('_')[0] not in media_timepoints]
    data_filt = data.filter(items=keep_columns)
    data_filt.to_csv(out_mat_file, sep='\t', compression="gzip")

    return None


def filter_for_ids(mat_file, keep_ids_file, gz_out_file, opposite=False, counts=False):
    """Given a list of ids, filter the  matrix file to keep
    only those ids. First column must be ids.
    """
    data_df = pd.read_csv(mat_file, sep='\t', header=0, index_col=0)
    keep_ids = pd.read_csv(keep_ids_file, header=None)
    if opposite:
    	keep_data = data_df.loc[~data_df.index.isin(keep_ids[0])]
    else:
    	keep_data = data_df.loc[data_df.index.isin(keep_ids[0])]
    if counts:
        keep_data = keep_data.applymap(int)
    keep_data.to_csv(gz_out_file, sep='\t', compression='gzip')

    return None


def get_ordered_subsample(in_file, out_file, out_nrow=2000):
    """Given an input text file, grab an ordered sample
    """
    num_lines = 0
    with open(in_file, "r") as fp:
        for line in fp:
            num_lines += 1

    skip = math.ceil(float(num_lines) / out_nrow)

    num_lines = 0
    with open(out_file, "w") as out:
        with open(in_file, "r") as fp:
            for line in fp:
                if num_lines % skip == 0:
                    out.write(line)
                num_lines += 1
    
    return None
