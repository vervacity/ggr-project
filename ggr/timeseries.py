"""Contains time series code
"""

import os
import pandas as pd

from scipy.stats import zscore
from scipy.stats import spearmanr

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel


def merge_cluster_files(cluster_files, out_file):
    """Given a list of cluster files, merge
    """
    clusters = pd.read_table(cluster_files[0], sep='\t')
    col_name = os.basename(cluster_files[0]).split("_optimal")[0]
    clusters.index = [col_name, "region"]
    
    for cluster_file in cluster_files[1:]:
        clusters_tmp = pd.read_table(cluster_file, sep='\t')
        col_name = os.basename(cluster_file).split("_optimal")[0]
        clusters_tmp.index = [col_name, "region"]

        clusters = clusters.merge(clusters_tmp, on="region")

    clusters.to_csv(out_file, sep='\t', header=False, index=False)

    return


def get_cluster_means(cluster_file, timeseries_file):
    """From cluster file and timeseries file, get 
    means of the clusters and return as numpy array
    """
    # open up cluster file
    clusters_df = pd.read_table(clusters_file, sep='\t')
    clusters_df.columns = ["cluster", "region"]
    cluster_num = list(set(clusters_df["cluster"].tolist()))

    # open up timeseries file
    timeseries_df = pd.read_table(timeseries_file, sep='\t')
    timeseries_df["regions"] = timeseries_df.index

    # merge the two
    merged_data = timeseries_df.merge(clusters_df, on="regions")

    # set up numpy out matrix
    cluster_means = np.zeros((cluster_num, timeseries_df.shape[1]))

    # get cluster means
    for cluster_idx in range(pooled_cluster_num):
        cluster_num = cluster_idx + 1

        cluster_data = pd.DataFrame(merged_data[merged_data["cluster"] == cluster_num])
        del cluster_data["cluster"]
        
        cluster_means[cluster_idx,:] = cluster_data.mean(axis=1)

    # if zscore, do so on rows
    cluster_means_z = zscore(cluster_means, axis=1)

    import ipdb
    ipdb.set_trace()
    
        
    return cluster_means_z


def get_corr_mat(mat_a, mat_b):
    """given two matrices, calculate correlations between
    the rows of the two
    """
    corr_mat = np.zeros((mat_a.shape[0], mat_b.shape[0]))

    for mat_a_idx in range(mat_a.shape[0]):
        for mat_b_idx in range(mat_b.shape[0]):
            corr_mat[mat_a_idx, mat_b_idx] = spearmanr(
                mat_a[mat_a_idx,:], mat_b[mat_b_idx,:])

    return corr_mat


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
        cluster_files.append(cluster_file)
        
        if not os.path.isdir(rep_dir):
            os.system("mkdir -p {}".format(rep_dir))
            cluster = (
                "DP_GP_cluster.py -i {} -o {} -p png --plot").format(
                    timeseries_files[i],
                    "{0}/{1}.{2}".format(
                        rep_dir, prefix,
                        timeseries_names[i]))
            print cluster
            dpgp_queue.put([os.system, [cluster]])
            
    # run queue
    run_in_parallel(dpgp_queue, parallel=3)

    # delete the unzipped tmp files
    os.system("rm {}/*.tmp".format(out_dir))

    # TODO save out all cluster nums into one file
    cluster_mat_file = "{}/{}.clusterings.txt".format(out_dir, prefix)
    if not os.path.isfile(cluster_mat_file):
        merge_cluster_files(cluster_files, cluster_mat_file)
    
    # now get cluster similarities to each other
    # for each cluster file, extract average profile of cluster

    # Set up cluster means
    cluster_files = []
    cluster_means = []
    for rep_idx in range(len(timeseries_files)):
        name_string = timeseries_names[rep_idx]
        cluster_files.append("{}/{}/{}.{}_optimal_clustering.txt".format(
            out_dir, name_string, prefix, name_string))
        cluster_means.append(get_cluster_means(
            cluster_files[name_string],
            timeseries_files[rep_idx]))
    
    # compute correlations, rep1-pooled and rep2-pooled
    true_reps = ["rep1", "rep2"]
    cluster_corrs = []
    cluster_max_matches = []
    for rep_idx in range(len(true_reps)):
        cluster_corrs.append(get_corr_mat(
            cluster_means[rep_idx], cluster_means[2]))
        cluster_max_matches.append(np.argmax(cluster_corrs[rep_idx], axis=1))
        
    # TODO: plot similarities out?
    
    
    # now read out to clustering file
    cluster_repr_mat_file = "{}/{}.clusterings.reproducible.txt".format(out_dir, prefix)
    with open(cluster_repr_mat_file, 'w') as out:
        with open(cluster_mat_file, 'r') as fp:
            for line in fp:
                fields = line.strip().split('\t')
                region = fields[0]
                rep1_cluster = int(fields[1]) + 1
                rep2_cluster = int(fields[2]) + 1
                pooled_cluster = int(fields[3]) + 1
                
                rep1_pooled_cluster = cluster_max_matches[0][rep1_cluster]
                rep2_pooled_cluster = cluster_max_matches[0][rep2_cluster]

                if rep1_pooled_cluster == rep2_pooled_cluster:
                    assert pooled_cluster == rep1_pooled_cluster
                    out.write(line)
                else:
                    # inconsistent
                    out.write("{}\t{}\t{}\t{}\n".format(region, rep1_cluster, rep2_cluster, "-1"))
    
    return cluster_repr_mat_file
