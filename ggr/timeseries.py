"""Contains time series code
"""

import os
import pandas as pd

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel

def get_consistent_dpgp_trajectories(
        rep1_timeseries_file,
        rep2_timeseries_file,
        pooled_timeseries_file,
        out_dir,
        prefix):
    """Given count matrices for replicates, gets 
    consistent trajectories
    """
    assert ".gz" in rep1_timeseries_file
    assert ".gz" in rep2_timeseries_file
    assert ".gz" in pooled_timeseries_file

    # setup names
    timeseries_zipped_files = [
        rep1_timeseries_file,
        rep2_timeseries_file,
        pooled_timeseries_file]

    timeseries_names = [
        "rep1",
        "rep2",
        "pooled"]

    # unzip files as needed
    timeseries_files = []
    for timeseries_file in timeseries_zipped_files:
        timeseries_unzipped_file = "{}/{}.tmp".format(
            out_dir,
            os.path.basename(timeseries_file).split('.mat')[0])
        unzip_mat = ("zcat {0} > {1}").format(
            timeseries_file, timeseries_unzipped_file)
        print unzip_mat
        os.system(unzip_mat)
        timeseries_files.append(timeseries_unzipped_file)
    
    # set up queue and add
    dpgp_queue = setup_multiprocessing_queue()
    cluster_files = []
    for i in range(len(timeseries_files)):
        rep_dir = "{}/{}".format(out_dir, timeseries_names[i])

        cluster_file = "{0}/{1}.{2}_optimal_clustering.txt".format(
            rep_dir, prefix, timeseries_names[i])
        cluster_files.append[cluster_file]
        
        if not os.path.isdir(rep_dir):
            os.system("mkdir -p {}".format(rep_dir))
            cluster = (
                "DP_GP_cluster.py "
                "-i {} "
                "-o {} "
                "-p png --plot").format(
                    timeseries_files[i],
                    "{0}/{1}.{2}".format(
                        rep_dir,
                        prefix,
                        timeseries_names[i]))

            dpgp_queue.put([os.system, [cluster]])
            
    # run queue
    run_in_parallel(dpgp_queue, parallel=3)

    # delete the unzipped tmp files
    os.system("rm {}/*.tmp".format(out_dir))
    
    # now get cluster similarities to each other
    # for each cluster file, extract average profile of cluster

    # compute correlation for each to each
    # and plot



    # then join all cluster nums (rep1, rep2, pooled) into file
    # then for each row, compare rep1/rep2 clusters to pool. if they overlap,
    # then keep. if they don't, label as inconclusive, -1


    # feed this out into a clustering file that represents consistent clusters


    # TODO return the filename of the optimal clustering
    
    return 
