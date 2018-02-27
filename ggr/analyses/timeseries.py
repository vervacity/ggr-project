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

from ggr.util.utils import run_shell_cmd


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
    # TODO make sure to remove examples that are not among the clusters
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


def filter_null_and_small_clusters_old(
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


def filter_null_and_small_clusters(
        cluster_file,
        mat_file,
        out_cluster_file,
        ci=0.999, # sort of multiple hypothesis correction
        size_cutoff=1000):
    """Quick filters to remove trajectory groups whose multivariate Gaussian
    does NOT reject the null (ie, vector of zeros falls in the confidence
    interval) as well as small clusters
    
    Args:
      mat_file: file with the timepoint data, with 0-column with ids
      cluster_file: file with clusters, 0-column is clusters, 1-column is ids
    """
    # first get sufficient statistics for the clusters
    cluster_means, cluster_covariances, cluster_sizes, cluster_names = get_cluster_sufficient_stats(
            cluster_file, mat_file)

    # set up size cutoff
    size_cutoff = size_cutoff * np.sum(cluster_sizes)
    
    # calculate what the chi2 cutoffs are based on the means and confidence intervals
    pdf_cutoff = chi2.pdf(1-ci, cluster_means[0].shape[0])
    indices_to_delete = []
    for cluster_idx in range(len(cluster_means)):

        # remove clusters by cluster size
        if cluster_sizes[cluster_idx] <= size_cutoff:
            indices_to_delete.append(cluster_idx)
            continue

        # calculate the multivariate normal (GP) and determine if cluster
        # does not reject the null.
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
    cluster_names_to_delete = (np.array(indices_to_delete) + 1).tolist()
    
    # return cluster file
    cluster_list = pd.read_table(cluster_file)
    cluster_list = cluster_list[~cluster_list["cluster"].isin(cluster_names_to_delete)]

    # and now renumber
    clusters_remaining = cluster_list["cluster"].unique().tolist()
    renumbering = dict(zip(clusters_remaining, range(1, len(clusters_remaining)+1)))
    print renumbering
    cluster_list["cluster"].replace(renumbering, inplace=True)
    cluster_list.columns = ["cluster", "id"]
    cluster_list = cluster_list.sort_values("cluster")

    # save out
    cluster_list.to_csv(out_cluster_file, sep="\t", index=False)
    
    return None



def get_consistent_soft_clusters_by_region(
        cluster_means,
        cluster_covariances,
        cluster_names,
        timepoint_vectors,
        pdf_cutoff,
        corr_cutoff,
        epsilon):
    """Given a list of timepoint vectors, compare to all clusters
    then perform an intersect to get consistent clusters
    utilizing a multivariate normal distribution
    """
    soft_cluster_sets = []
    
    for timepoint_vector in timepoint_vectors:
        cluster_set = []
        
        for cluster_idx in range(len(cluster_means)):

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
            cluster_set.append(cluster_idx)

        # append set to all cluster sets
        soft_cluster_sets.append(set(cluster_set))
        
    # intersect
    consistent_cluster_indices = list(set.intersection(*soft_cluster_sets))
    consistent_clusters = [cluster_names[cluster_idx]
                           for cluster_idx in consistent_cluster_indices]
    
    # With the consistent clusters, figure out best fit (to have a hard cluster assignment)
    best_pdf_val = 0
    hard_cluster = None
    for cluster_idx in consistent_cluster_indices:
        # get pdf val
        pdf_val = multivariate_normal.pdf(
            timepoint_vectors[2],
            mean=cluster_means[cluster_idx],
            cov=cluster_covariances[cluster_idx],
            allow_singular=True)
        if pdf_val > best_pdf_val:
            best_pdf_val = pdf_val
            hard_cluster = cluster_names[cluster_idx]

    # TODO(dk) now filter consistent clusters - if they are within {error} of pdf val, then keep
    # this would mean the region is fairly close to both
    filtered_consistent_clusters = []
    for cluster_idx in consistent_cluster_indices:
        # get pdf val
        pdf_val = multivariate_normal.pdf(
            timepoint_vectors[2],
            mean=cluster_means[cluster_idx],
            cov=cluster_covariances[cluster_idx],
            allow_singular=True)
        if (best_pdf_val - pdf_val) < (epsilon * best_pdf_val):
            filtered_consistent_clusters.append(cluster_names[cluster_idx])
    
    #return consistent_clusters, hard_cluster
    return filtered_consistent_clusters, hard_cluster


def get_consistent_soft_clusters_old(
        rep1_timeseries_file,
        rep2_timeseries_file,
        pooled_timeseries_file,
        out_dir,
        prefix,
        cluster_means,
        cluster_covariances,
        cluster_names,
        ci=0.95,
        corr_cutoff=0.05,
        epsilon=0.10):
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

                    # now for each cluster, check consistency and get soft/hard clusters
                    consistent_clusters, hard_cluster = get_consistent_soft_clusters_by_region(
                        cluster_means,
                        cluster_covariances,
                        cluster_names,
                        [rep1_timepoints, rep2_timepoints, pooled_timepoints],
                        pdf_cutoff,
                        corr_cutoff,
                        epsilon)

                    # write to cluster file
                    
                    for cluster in consistent_clusters:
                        
                        # write out to file
                        cluster_file = "{}/soft/{}.cluster_{}.soft.txt.gz".format(out_dir, prefix, cluster)
                        with gzip.open(cluster_file, 'a') as out:
                            out.write("{}\t{}\n".format(
                                pooled_region,
                                "\t".join(map(str, pooled_timepoints.tolist()))))
                            
                    # save out hard clusters
                    if hard_cluster is not None:
                        hard_cluster_file = "{}/hard/{}.clusters.hard.all.txt.gz".format(out_dir, prefix)
                        with gzip.open(hard_cluster_file, "a") as out:
                            out.write("{}\t{}\t{}\n".format(
                                pooled_region,
                                "\t".join(map(str, pooled_timepoints.tolist())),
                                hard_cluster))

                        # Write out to individual hard cluster file
                        cluster_file = "{}/hard/{}.cluster_{}.hard.txt.gz".format(out_dir, prefix, hard_cluster)
                        with gzip.open(cluster_file, 'a') as out:
                            out.write("{}\t{}\n".format(
                                pooled_region,
                                "\t".join(map(str, pooled_timepoints.tolist()))))
                        
                    region_idx += 1
                    if region_idx % 1000 == 0:
                        print region_idx

    return


def get_reproducible_clusters(
        clusters_file,
        pooled_mat_file,
        rep1_timeseries_file,
        rep2_timeseries_file,
        out_soft_clusters_file,
        out_hard_clusters_file,
        out_dir,
        prefix,
        ci=0.95,
        corr_cutoff=0.05,
        epsilon=0.10):
    """Given rep1/rep2/pooled timeseries files, go through
    regions and check for consistency after soft clustering
    """
    run_shell_cmd("mkdir -p {0}/soft {0}/hard".format(out_dir))
    cluster_means, cluster_covariances, cluster_sizes, cluster_names = get_cluster_sufficient_stats(
        clusters_file, pooled_mat_file)

    # get pdf val for confidence interval
    # note that the multivariate normal pdf is distributed
    # as chi2(alpha) where alpha is the chosen significance
    pdf_cutoff = chi2.pdf(1-ci, cluster_means[0].shape[0])

    with open(out_hard_clusters_file, "a") as out:
        out.write("cluster\tid\n")
    with open(out_soft_clusters_file, "a") as out:
        out.write("cluster\tid\n")
    
    # open all files to stream regions
    region_idx = 0
    with gzip.open(rep1_timeseries_file, 'r') as rep1:
        with gzip.open(rep2_timeseries_file, 'r') as rep2:
            with gzip.open(pooled_mat_file, 'r') as pooled:

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
                    rep1_region, rep1_timepoints = rep1_fields[0], zscore(
                        np.array(map(float, rep1_fields[1:])))
                    rep2_region, rep2_timepoints = rep2_fields[0], zscore(
                        np.array(map(float, rep2_fields[1:])))
                    pooled_region, pooled_timepoints = pooled_fields[0], zscore(
                        np.array(map(float, pooled_fields[1:])))

                    assert rep1_region == rep2_region
                    assert rep1_region == pooled_region
                    assert rep2_region == pooled_region

                    # now for each cluster, check consistency and get soft/hard clusters
                    # TODO: rename this function - assign_region_to_cluster
                    consistent_clusters, hard_cluster = get_consistent_soft_clusters_by_region(
                        cluster_means,
                        cluster_covariances,
                        cluster_names,
                        [rep1_timepoints, rep2_timepoints, pooled_timepoints],
                        pdf_cutoff,
                        corr_cutoff,
                        epsilon)

                    # TODO: save out to cluster files instead of splitting
                    if len(consistent_clusters) > 0:
                        with open(out_soft_clusters_file, "a") as out:
                            for cluster in consistent_clusters:
                                out.write("{}\t{}\n".format(cluster, pooled_region))

                    if hard_cluster is not None:
                        with open(out_hard_clusters_file, "a") as out:
                            out.write("{}\t{}\n".format(hard_cluster, pooled_region))

                    region_idx += 1
                    if region_idx % 1000 == 0:
                        print region_idx

    return None



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


def run_dpgp(mat_file, prefix, out_dir, tmp_dir, subsample=False, subsample_num=5000):
    """Run DP-GP while managing files correctly
    """
    run_shell_cmd("mkdir -p {}".format(out_dir))
    run_shell_cmd("mkdir -p {}".format(tmp_dir))

    # unzip file if needed
    input_mat_file = "{}/{}.dpgp.tmp".format(
        tmp_dir, os.path.basename(mat_file).split(".txt")[0])
    if mat_file.endswith(".gz"):
        unzip_mat = "zcat {0} > {1}".format(
            mat_file, input_mat_file)
        run_shell_cmd(unzip_mat)
    else:
        unzip_mat = "cat {0} > {1}".format(
            mat_file, input_mat_file)
        run_shell_cmd(unzip_mat)

    # if subsample
    if subsample:
        subsampled_file = "{}.subsampled.tmp".format(input_mat_file.split(".tmp")[0])
        # keep header
        echo_header = "cat {0} | awk 'NR < 2' > {1}".format(
            input_mat_file, subsampled_file)
        os.system(echo_header)
        # subsample
        run_subsample = "cat {0} | awk 'NR > 1' | shuf -n {1} >> {2}".format(
            input_mat_file, subsample_num, subsampled_file)
        os.system(run_subsample)
        os.system("rm {}".format(input_mat_file))
        input_mat_file = subsampled_file
        
    # TODO change header if need be?
    cluster = (
        "DP_GP_cluster.py -i {} -o {} --fast -p pdf --plot").format(
            input_mat_file,
            "{0}/{1}".format(out_dir, prefix))
    run_shell_cmd(cluster)

    # delete the unzipped tmp files
    run_shell_cmd("rm {}/*.dpgp.tmp".format(tmp_dir))

    # outfile name
    optimal_clusters = "{}/{}_optimal_clustering.txt".format(out_dir, prefix)

    return optimal_clusters


# TODO this should be a workflow, not an analysis
def get_consistent_dpgp_trajectories(
        matrices,
        out_dir,
        prefix,
        raw_cluster_min_size=1000,
        raw_cluster_reject_null_ci_interval=0.999,
        rep_to_cluster_ci_interval=0.95,        
        rep_to_cluster_corr_cutoff=0.05,
        epsilon=0.10):
    """Given count matrices for replicates, gets 
    consistent trajectories

    Do this by calling trajectories using the pooled data
    (to increase "read depth") this leads to more fine grained
    trajectories. Then, for each region, get the max similarity
    to a trajectory (mean val). If consistent, keep. Else, 
    throw away.


    """
    run_shell_cmd("mkdir -p {}".format(out_dir))
    
    rep1_timeseries_file = matrices[0]
    rep2_timeseries_file = matrices[1]
    pooled_timeseries_file = matrices[2]
    
    assert ".gz" in rep1_timeseries_file
    assert ".gz" in rep2_timeseries_file
    assert ".gz" in pooled_timeseries_file

    # run DP GP on pooled data
    pooled_dir = "{}/pooled".format(out_dir)
    pooled_cluster_file = "{0}/{1}.pooled_optimal_clustering.txt".format(
        pooled_dir, prefix)
    if not os.path.isdir(pooled_dir):

        # unzip files as needed
        pooled_unzipped_file = "{}/{}.tmp".format(
            out_dir,
            os.path.basename(pooled_timeseries_file).split(".mat")[0])
        unzip_mat = ("zcat {0} > {1}").format(
            pooled_timeseries_file, pooled_unzipped_file)
        run_shell_cmd(unzip_mat)

        # TODO(dk) convert header as needed?
        

        # cluster
        run_shell_cmd("mkdir -p {}".format(pooled_dir))
        cluster = (
            "DP_GP_cluster.py -i {} -o {} -p pdf --plot --fast").format(
                pooled_unzipped_file,
                "{0}/{1}.pooled".format(
                    pooled_dir, prefix))
        run_shell_cmd(cluster)

        # delete the unzipped tmp files
        run_shell_cmd("rm {}/*.tmp".format(out_dir))

    # stage 2: soft clustering
    # TODO(dk) separate this out from DPGP consistency?
    soft_cluster_dir = "{}/soft".format(out_dir)
    hard_cluster_dir = "{}/hard".format(out_dir)
    if not os.path.isdir(soft_cluster_dir):
        run_shell_cmd("mkdir -p {}".format(soft_cluster_dir))
        run_shell_cmd("mkdir -p {}".format(hard_cluster_dir))
    
        # for each cluster (from POOLED clusters), get the mean and covariance matrix (assuming multivariate normal)
        cluster_means, cluster_covariances, cluster_sizes, cluster_names = get_cluster_sufficient_stats(
            pooled_cluster_file, pooled_timeseries_file)

        # filter clusters: remove low membership (<1k) clusters and those that don't reject the null
        # note that "low membership" was empirical - strong split between <500 and >1400.
        # rejecting null - confidence interval of .999 (this is stronger as an FDR correction)
        # TODO(dk) pull out params here
        cluster_means, cluster_covariances, cluster_sizes, cluster_names = filter_null_and_small_clusters(
            cluster_means,
            cluster_covariances,
            cluster_sizes,
            cluster_names,
            ci=raw_cluster_reject_null_ci_interval,
            size_cutoff=raw_cluster_min_size)
        
        # now extract consistent soft clusters
        # ie, for each region, check each cluster to see if it belongs after checking
        # confidence interval (multivariate normal), pearson and spearmans.
        # TODO(dk) pull out params here
        get_consistent_soft_clusters(
            rep1_timeseries_file,
            rep2_timeseries_file,
            pooled_timeseries_file,
            out_dir,
            prefix,
            cluster_means,
            cluster_covariances,
            cluster_names,
            ci=rep_to_cluster_ci_interval,
            corr_cutoff=rep_to_cluster_corr_cutoff,
            epsilon=epsilon)
    
    # TODO separate this out as a viz analysis
    # sanity check - replot the clusters
    plot_dir = "{}/soft/plot".format(out_dir)
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))
        # use R to plot 
        plot_soft_clusters = ("plot_soft_trajectories.R {0}/soft/plot {0}/soft/{1}").format(
            out_dir, "*soft*.gz")
        print plot_soft_clusters
        run_shell_cmd(plot_soft_clusters)

    plot_dir = "{}/hard/plot".format(out_dir)
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))
        # use R to plot 
        plot_hard_clusters = ("plot_soft_trajectories.R {0}/hard/plot {0}/hard/{1}").format(
            out_dir, "*cluster_*hard*.gz")
        print plot_hard_clusters
        run_shell_cmd(plot_hard_clusters)

    # TODO(dk): return the soft and hard clusters as lists
        
        
    return


def split_mat_to_clusters(cluster_file, cluster_mat, out_dir, prefix):
    """helper function to split a mat (given a cluster file)
    into cluster mat files and also lists
    """
    # use pandas
    clusters = pd.read_table(cluster_file, header=True)
    data = pd.read_table(cluster_mat, header=True)

    import ipdb
    ipdb.set_trace()

    return None


def plot_clusters(
        cluster_file,
        cluster_subsample_file,
        cluster_mat,
        out_dir,
        prefix):
    """plots clusters given in the cluster file using the
    data in cluster mat
    """
    # assertions
    assert os.path.isdir(out_dir)
    
    # heatmap plot
    r_plot_heatmap = (
        "viz.plot_timeseries_heatmap.R "
        "{0} {1} {2} {3}").format(
            cluster_subsample_file, cluster_mat, out_dir, prefix)
    run_shell_cmd(r_plot_heatmap)

    # individual clusters
    r_plot_clusters = (
        "viz.plot_timeseries_clusters.R "
        "{0} {1} {2} {3}").format(
            cluster_file, cluster_mat, out_dir, prefix)
    run_shell_cmd(r_plot_clusters)

    # print "plot_clusters function (analyses/timeseries.py)"
    
    return None


def flip_positions(position_list, start, middle, stop):
    """for reorder clusers below: flip from start:stop around middle.
    """
    tmp = position_list[start:stop]
    position_list[start:middle] = tmp[middle:stop]
    position_list[middle:stop] = tmp[start:middle]
    
    return position_list


def get_ordered_tree_nodes(tree):
    """Recursively go through tree to collect nodes
    This will be in the order of leaveslist
    """
    nodes = []
    # check left
    if tree.left.count == 1:
        nodes.append(tree.left.id)
    else:
        nodes += get_ordered_tree_nodes(tree.left)
        
    # check right
    if tree.right.count == 1:
        nodes.append(tree.right.id)
    else:
        nodes += get_ordered_tree_nodes(tree.right)
        
    # sum up and return
    return nodes


def reorder_tree(tree, cluster_means):
    """Recursive tool to reorder tree
    """
    # go left
    if tree.left.count == 1:
        tree.left = tree.left
        left_nodes = [tree.left.id]
    else:
        # adjust the tree
        tree.left = reorder_tree(tree.left, cluster_means)
        left_nodes = get_ordered_tree_nodes(tree.left)
        
    # go right
    if tree.right.count == 1:
        tree.right = tree.right
        right_nodes = [tree.right.id]
    else:
        # adjust the tree
        tree.right = reorder_tree(tree.right, cluster_means)
        right_nodes = get_ordered_tree_nodes(tree.right)
        
    # calculate average cluster means for each set
    left_cluster_mean = np.sum(cluster_means[left_nodes,:], axis=0)
    right_cluster_mean = np.sum(cluster_means[right_nodes,:], axis=0)
    
    # extract the max
    left_max_idx = np.argmax(left_cluster_mean)
    right_max_idx = np.argmax(right_cluster_mean)

    # if max is at the edges, calculate slope to nearest 0
    flip = False
    if left_max_idx != right_max_idx:
        # good to go, carry on
        if left_max_idx > right_max_idx:
            flip = True
    else:
        # if left edge:
        if left_max_idx == 0:
            left_slope = left_cluster_mean[0] - left_cluster_mean[3]
            right_slope = right_cluster_mean[0] - right_cluster_mean[3]
            if left_slope < right_slope:
                flip = False
        # if right:
        elif left_max_idx == cluster_means.shape[1] - 1:
            left_slope = left_cluster_mean[-1] - left_cluster_mean[-3]
            right_slope = right_cluster_mean[-1] - right_cluster_mean[-3]
            if left_slope > right_slope:
                flip = True
        # if middle:
        else:
            left_side_max_idx = np.argmax(left_cluster_mean[[left_max_idx-1, left_max_idx+1]])
            right_side_max_idx = np.argmax(right_cluster_mean[[right_max_idx-1, right_max_idx+1]])
            if left_side_max_idx > right_side_max_idx:
                flip = True


    #import ipdb
    #ipdb.set_trace()
    
    # reorder accordingly
    if flip == True:
        right_tmp = tree.right
        left_tmp = tree.left
        tree.right = left_tmp
        tree.left = right_tmp
    
    return tree


def reorder_clusters(cluster_file, cluster_mat, out_cluster_file):
    """Sort clusters by hclust similarity (hierarchical clustering in sklearn)
     """
    from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, to_tree
    #from scipy.spatial.distance import pdist, squareform
    from scipy.stats import zscore
    
    # first extract cluster means
    cluster_means, cluster_covariances, cluster_sizes, cluster_names = get_cluster_sufficient_stats(
        cluster_file, cluster_mat)
    
    means_z = zscore(np.array(cluster_means), axis=0)
    
    #cluster_dist = pdist(np.array(cluster_means), "euclidean")
    hclust = linkage(means_z, method="ward")

    # this is all the reordering code below
    if False:
        # using leaves_list and fcluster, determine split and reverse the FIRST half
        top_cut = fcluster(hclust, 2, criterion="maxclust")
        ordered_leaves = leaves_list(hclust)
        for i in xrange(ordered_leaves.shape[0]):
            current_leaf = ordered_leaves[i]
            if top_cut[current_leaf] == 2:
                # found the leaf
                split_point = i
                break

        # take left side of dendrogram and reverse
        # a recursive reordering of leaves by weight?
        ordered_leaves[0:split_point] = np.flip(ordered_leaves[0:split_point], axis=0)
        #ordered_leaves[split_point:] = np.flip(ordered_leaves[split_point:], axis=0)
        print ordered_leaves + 1
    else:
        # try a recursive reordering
        hclust_tree = to_tree(hclust)
        old_ordered_leaves = leaves_list(hclust)
        
        reordered_tree = reorder_tree(hclust_tree, np.array(cluster_means))
        ordered_leaves = np.array(get_ordered_tree_nodes(reordered_tree))

        print old_ordered_leaves
        print ordered_leaves
        
    # build renumbering dict
    renumbering = dict(zip((ordered_leaves+1).tolist(), range(1, len(ordered_leaves)+1)))
    print renumbering
    
    # read in cluster file
    cluster_list = pd.read_table(cluster_file)

    # renumber and sort
    cluster_list["cluster"].replace(renumbering, inplace=True)
    cluster_list = cluster_list.sort_values("cluster")
    cluster_list.to_csv(out_cluster_file, sep="\t", index=False)
    
    return None


def split_clusters(cluster_file):
    """Split cluster file into files of ids per cluster
    """

    cluster_data = pd.read_table(cluster_file)
    cluster_data.columns = ["cluster", "id"]
    cluster_names = cluster_data["cluster"].unique().tolist()

    for cluster_name in cluster_names:
        # get that subset
        single_cluster = cluster_data.loc[cluster_data["cluster"] == cluster_name]
        
        # and write out
        out_file = "{}.cluster_{}.txt.gz".format(
            cluster_file.split(".clustering")[0],
            cluster_name)
        single_cluster.to_csv(
            out_file,
            columns=["id"],
            compression="gzip",
            sep='\t',
            header=False,
            index=False)
    
    return None
