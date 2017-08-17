"""Contains time series code
"""

import os
import gzip
import numpy as np
import pandas as pd

from scipy.stats import zscore
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import multivariate_normal
from scipy.stats import chi2

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel


def merge_cluster_files(cluster_files, out_file, cluster_name_pos=2):
    """Given a list of cluster files, merge
    """
    # merge all clusters
    clusters = None
    for cluster_file in cluster_files:
        clusters_tmp = pd.read_table(cluster_file, sep='\t', index_col=0, header=None)
        col_name = os.path.basename(cluster_file).split(".")[cluster_name_pos]
        clusters_tmp[col_name] = [1 for i in range(clusters_tmp.shape[0])]
        clusters_tmp["region"] = clusters_tmp.index
        clusters_tmp = clusters_tmp[["region", col_name]]    
        
        if clusters is None:
            clusters = clusters_tmp
        else:
            clusters = clusters.merge(clusters_tmp, how="outer", on="region")

    clusters = clusters.fillna(0)
    clusters.index = clusters["region"]
    del clusters["region"]
    
    # and set up master cluster id (for sorting)
    clusters["atac_cluster"] = [0 for i in range(clusters.shape[0])]
    num_clusters = clusters.shape[1] - 1
    for i in range(num_clusters):
        clusters["atac_cluster"] = clusters["atac_cluster"] + (10**i) * (clusters.iloc[:,i])
    clusters["atac_cluster"] = clusters["atac_cluster"].astype(int)

    # and sort by ID column
    clusters = clusters.sort_values("atac_cluster")
    clusters.to_csv(out_file, sep='\t', compression="gzip")

    return


def get_cluster_means(cluster_file, timeseries_file):
    """From cluster file and timeseries file, get 
    means of the clusters and return as numpy array
    """
    # open up cluster file
    clusters_df = pd.read_table(cluster_file, sep='\t')
    clusters_df.columns = ["cluster", "regions"]
    cluster_nums = list(set(clusters_df["cluster"].tolist()))
    
    # open up timeseries file
    timeseries_df = pd.read_table(timeseries_file, sep='\t', index_col=0)
    timeseries_df["regions"] = timeseries_df.index
    
    # merge the two
    merged_data = timeseries_df.merge(clusters_df, on="regions")
    del merged_data["regions"]
    
    # set up numpy out matrix
    cluster_means = np.zeros((len(cluster_nums), merged_data.shape[1]-1))

    # get cluster means
    for cluster_idx in range(len(cluster_nums)):
        cluster_num = cluster_idx + 1

        cluster_data = pd.DataFrame(merged_data[merged_data["cluster"] == cluster_num])
        del cluster_data["cluster"]
        
        cluster_means[cluster_idx,:] = cluster_data.mean(axis=0)

    # if zscore, do so on rows
    cluster_means_z = zscore(cluster_means, axis=1)
        
    return cluster_means_z


def get_cluster_sufficient_stats(cluster_file, timeseries_file):
    """From cluster file and timeseries file, get 
    means and covariances of the clusters and return as numpy array
    """
    # open up cluster file
    clusters_df = pd.read_table(cluster_file, sep='\t')
    clusters_df.columns = ["cluster", "regions"]
    cluster_nums = list(set(clusters_df["cluster"].tolist()))
    
    # open up timeseries file and zscore
    timeseries_df = pd.read_table(timeseries_file, sep='\t', index_col=0)
    timeseries_df = timeseries_df.apply(zscore, axis=1)
    timeseries_df["regions"] = timeseries_df.index
    
    # merge the two
    merged_data = timeseries_df.merge(clusters_df, on="regions")
    del merged_data["regions"]
    
    # set up numpy out matrix
    cluster_means = []
    cluster_covariances = []
    cluster_sizes = []
    cluster_names = []
    
    # per cluster, get info
    for cluster_idx in range(len(cluster_nums)):
        cluster_num = cluster_idx + 1

        cluster_data = pd.DataFrame(merged_data[merged_data["cluster"] == cluster_num])
        del cluster_data["cluster"]

        # mean and covariance
        cluster_means.append(
            cluster_data.mean(axis=0).as_matrix())
        
        cluster_covariances.append(
            cluster_data.cov().as_matrix()) # note pandas already normalizes by N-1

        # other useful info
        cluster_sizes.append(cluster_data.shape[0])
        cluster_names.append(cluster_num)
        
    return cluster_means, cluster_covariances, cluster_sizes, cluster_names


def filter_null_and_small_clusters(
        cluster_means,
        cluster_covariances,
        cluster_sizes,
        cluster_names,
        ci=0.999, # sort of multiple hypothesis correction
        size_cutoff=1000):
    """Quick filters to remove trajectory groups whose multivariate Gaussian
    does NOT reject the null (ie, vector of zeros falls in the confidence
    interval) as well as small clusters
    """
    pdf_cutoff = chi2.pdf(1-ci, cluster_means[0].shape[0])
    indices_to_delete = []
    for cluster_idx in range(len(cluster_means)):

        if cluster_sizes[cluster_idx] <= size_cutoff:
            indices_to_delete.append(cluster_idx)
            continue

        pdf_val = multivariate_normal.pdf(
            np.array([0 for i in range(cluster_means[0].shape[0])]),
            mean=cluster_means[cluster_idx],
            cov=cluster_covariances[cluster_idx],
            allow_singular=True)
        
        if pdf_val > pdf_cutoff:
            indices_to_delete.append(cluster_idx)
            continue

    print cluster_sizes
    print indices_to_delete
    
    # fix clusters
    for index in sorted(indices_to_delete, reverse=True):
        del cluster_means[index]
        del cluster_covariances[index]
        del cluster_sizes[index]
        del cluster_names[index]
    
    return cluster_means, cluster_covariances, cluster_sizes, cluster_names



def get_consistent_soft_clusters_by_region(
        cluster_means,
        cluster_covariances,
        cluster_names,
        timepoint_vectors,
        pdf_cutoff,
        corr_cutoff=0.05):
    """Given a list of timepoint vectors, compare to all clusters
    then perform an intersect to get consistent clusters
    utilizing a multivariate normal distribution
    """
    soft_cluster_sets = []
    
    for timepoint_vector in timepoint_vectors:
        cluster_set = []
        
        for cluster_idx in range(len(cluster_means)):
            cluster_num = cluster_names[cluster_idx]

            # determine if in confidence interval
            pdf_val = multivariate_normal.pdf(
                timepoint_vector,
                mean=cluster_means[cluster_idx],
                cov=cluster_covariances[cluster_idx],
                allow_singular=True) # this has to do with scaling? stack overflow 35273908

            if pdf_val < pdf_cutoff: # if not in confidence interval, continue
                continue

            # pearson: check relative differences
            pearson_corr, pearson_pval = pearsonr(
                timepoint_vector, cluster_means[cluster_idx])
            
            if pearson_pval > corr_cutoff: # if pval is not significant, continue
                continue
            
            # spearman: check rank order of timepoints
            spearman_corr, spearman_pval = spearmanr(
                timepoint_vector, cluster_means[cluster_idx])

            if spearman_pval > corr_cutoff: # if pval is not significant, continue
                continue

            # if all conditions are met, add to cluster set
            cluster_set.append(cluster_num)

        # append set to all cluster sets
        soft_cluster_sets.append(set(cluster_set))
        
    # intersect
    consistent_clusters = list(set.intersection(*soft_cluster_sets))

    #import ipdb
    #ipdb.set_trace()
    
    return consistent_clusters


def get_consistent_soft_clusters(
        rep1_timeseries_file,
        rep2_timeseries_file,
        pooled_timeseries_file,
        out_dir,
        prefix,
        cluster_means,
        cluster_covariances,
        cluster_names,
        ci=0.95,
        corr_cutoff=0.05):
    """Given rep1/rep2/pooled timeseries files, go through
    regions and check for consistency after soft clustering
    """
    # get pdf val for confidence interval
    # note that the multivariate normal pdf is distributed
    # as chi2(alpha) where alpha is the chosen significance
    pdf_cutoff = chi2.pdf(1-ci, cluster_means[0].shape[0])

    # open all files to stream regions
    region_idx = 0
    with gzip.open(rep1_timeseries_file, 'r') as rep1:
        with gzip.open(rep2_timeseries_file, 'r') as rep2:
            with gzip.open(pooled_timeseries_file, 'r') as pooled:

                while True:
                    
                    # read lines
                    rep1_line = rep1.readline()
                    rep2_line = rep2.readline()
                    pooled_line = pooled.readline()
                
                    # break if end of file
                    if rep1_line == "":
                        break
                
                    # ignore header line
                    if rep1_line.strip().startswith("d"):
                        continue

                    # separate to fields
                    rep1_fields = rep1_line.strip().split("\t")
                    rep2_fields = rep2_line.strip().split("\t")
                    pooled_fields = pooled_line.strip().split("\t")

                    # split out regions and timepoints
                    rep1_region, rep1_timepoints = rep1_fields[0], zscore(np.array(map(float, rep1_fields[1:])))
                    rep2_region, rep2_timepoints = rep2_fields[0], zscore(np.array(map(float, rep2_fields[1:])))
                    pooled_region, pooled_timepoints = pooled_fields[0], zscore(np.array(map(float, pooled_fields[1:])))

                    assert rep1_region == rep2_region
                    assert rep1_region == pooled_region
                    assert rep2_region == pooled_region

                    # now for each cluster, check consistency
                    consistent_clusters = get_consistent_soft_clusters_by_region(
                        cluster_means,
                        cluster_covariances,
                        cluster_names,
                        [rep1_timepoints, rep2_timepoints, pooled_timepoints],
                        pdf_cutoff,
                        corr_cutoff=corr_cutoff)
                    
                    for cluster in consistent_clusters:
                        
                        # write out to file
                        cluster_file = "{}/soft/{}.cluster_{}.soft.txt.gz".format(out_dir, prefix, cluster)
                        with gzip.open(cluster_file, 'a') as out:
                            out.write("{}\t{}\n".format(
                                pooled_region,
                                "\t".join(map(str, pooled_timepoints.tolist()))))

                    region_idx += 1
                    if region_idx % 1000 == 0:
                        print region_idx

    return


def get_corr_mat(mat_a, mat_b):
    """given two matrices, calculate correlations between
    the rows of the two
    """
    corr_mat = np.zeros((mat_a.shape[0], mat_b.shape[0]))

    for mat_a_idx in range(mat_a.shape[0]):
        for mat_b_idx in range(mat_b.shape[0]):
            corr_mat[mat_a_idx, mat_b_idx] = spearmanr(
                mat_a[mat_a_idx,:], mat_b[mat_b_idx,:])[0]

    return corr_mat


def get_consistent_dpgp_trajectories(
        rep1_timeseries_file,
        rep2_timeseries_file,
        pooled_timeseries_file,
        out_dir,
        prefix):
    """Given count matrices for replicates, gets 
    consistent trajectories

    Do this by calling trajectories using the pooled data
    (to increase "read depth") this leads to more fine grained
    trajectories. Then, for each region, get the max similarity
    to a trajectory (mean val). If consistent, keep. Else, 
    throw away.


    """
    assert ".gz" in rep1_timeseries_file
    assert ".gz" in rep2_timeseries_file
    assert ".gz" in pooled_timeseries_file

    # unzip files as needed
    pooled_unzipped_file = "{}/{}.tmp".format(
        out_dir,
        os.path.basename(pooled_timeseries_file).split(".mat")[0])
    unzip_mat = ("zcat {0} > {1}").format(
        pooled_timeseries_file, pooled_unzipped_file)
    print unzip_mat
    os.system(unzip_mat)

    # run DP GP on pooled data
    pooled_dir = "{}/pooled".format(out_dir)
    pooled_cluster_file = "{0}/{1}.pooled_optimal_clustering.txt".format(
        pooled_dir, prefix)
    if not os.path.isdir(pooled_dir):
        os.system("mkdir -p {}".format(pooled_dir))
        cluster = (
            "DP_GP_cluster.py -i {} -o {} -p png --plot").format(
                pooled_unzipped_file,
                "{0}/{1}.pooled".format(
                    pooled_dir, prefix))
        print cluster
        os.system(cluster)

    # delete the unzipped tmp files
    os.system("rm {}/*.tmp".format(out_dir))

    # stage 2: soft clustering
    soft_cluster_dir = "{}/soft".format(out_dir)
    if not os.path.isdir(soft_cluster_dir):
        os.system("mkdir -p {}".format(soft_cluster_dir))
    
        # for each cluster (from POOLED clusters), get the mean and covariance matrix (assuming multivariate normal)
        cluster_means, cluster_covariances, cluster_sizes, cluster_names = get_cluster_sufficient_stats(
            pooled_cluster_file, pooled_timeseries_file)

        # filter clusters: remove low membership (<1k) clusters and those that don't reject the null
        # note that "low membership" was empirical - strong split between <500 and >1400.
        # rejecting null - confidence interval of .999 (this is stronger as an FDR correction)
        cluster_means, cluster_covariances, cluster_sizes, cluster_names = filter_null_and_small_clusters(
            cluster_means,
            cluster_covariances,
            cluster_sizes,
            cluster_names)
        
        # now extract consistent soft clusters
        # ie, for each region, check each cluster to see if it belongs after checking
        # confidence interval (multivariate normal), pearson and spearmans.
        get_consistent_soft_clusters(
            rep1_timeseries_file,
            rep2_timeseries_file,
            pooled_timeseries_file,
            out_dir,
            prefix,
            cluster_means,
            cluster_covariances,
            cluster_names)
                            
    # sanity check - replot the soft clusters
    plot_dir = "{}/soft/plot".format(out_dir)
    if not os.path.isdir(plot_dir):
        os.system("mkdir -p {}".format(plot_dir))
        
        # use R to plot 
        plot_soft_clusters = ("plot_soft_trajectories.R {0}/soft/plot {0}/soft/{1}").format(
            out_dir, "*.gz")
        print plot_soft_clusters
        os.system(plot_soft_clusters)

    # TODO convert into BED files with cluster numbers?
        

    return
