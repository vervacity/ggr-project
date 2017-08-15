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


def merge_cluster_files(cluster_files, out_file):
    """Given a list of cluster files, merge
    """
    clusters = pd.read_table(cluster_files[0], sep='\t')
    col_name = os.path.basename(cluster_files[0]).split("_optimal")[0]
    clusters.columns = [col_name, "region"]
    clusters = clusters[["region", col_name]]
    
    for cluster_file in cluster_files[1:]:
        clusters_tmp = pd.read_table(cluster_file, sep='\t')
        col_name = os.path.basename(cluster_file).split("_optimal")[0]
        clusters_tmp.columns = [col_name, "region"]

        clusters = clusters.merge(clusters_tmp, on="region")

    clusters.to_csv(out_file, sep='\t', index=False)

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
    
    # get cluster means
    for cluster_idx in range(len(cluster_nums)):
        cluster_num = cluster_idx + 1

        cluster_data = pd.DataFrame(merged_data[merged_data["cluster"] == cluster_num])
        del cluster_data["cluster"]
        
        cluster_means.append(
            cluster_data.mean(axis=0).as_matrix())
        
        cluster_covariances.append(
            cluster_data.cov().as_matrix())
        
    return cluster_means, cluster_covariances


def get_consistent_soft_clusters(
        cluster_means,
        cluster_covariances,
        timepoint_vectors,
        pdf_cutoff,
        corr_cutoff=0.05):
    """Given a list of timepoint vectors, compare to all clusters
    then perform an intersect to get consistent clusters
    utilizing a multivariate normal distribution
    """

    # use a chi square to figure out sig area?
    
    soft_cluster_sets = []
    
    for timepoint_vector in timepoint_vectors:
        cluster_set = []
        
        for cluster_idx in range(len(cluster_means)):
            cluster_num = cluster_idx + 1

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


def get_consistent_dpgp_trajectories_v1(
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
            cluster_files[rep_idx],
            timeseries_zipped_files[rep_idx]))
    
    # compute correlations, rep1-pooled and rep2-pooled
    true_reps = ["rep1", "rep2"]
    cluster_corrs = []
    cluster_max_matches = []
    for rep_idx in range(len(true_reps)):
        cluster_corrs.append(get_corr_mat(
            cluster_means[rep_idx], cluster_means[2]))
        cluster_max_matches.append(np.argmax(cluster_corrs[rep_idx], axis=1) + 1)
        
    # TODO: plot similarities out?
    # save out cluster corr matrices

    # now read out to clustering file
    cluster_repr_mat_file = "{}/{}.clusterings.reproducible.txt".format(out_dir, prefix)
    with open(cluster_repr_mat_file, 'w') as out:
        with open(cluster_mat_file, 'r') as fp:
            for line in fp:

                if line.startswith("region"):
                    out.write(line)
                    continue
                
                fields = line.strip().split('\t')
                region = fields[0]
                rep1_cluster = int(fields[1]) - 1
                rep2_cluster = int(fields[2]) - 1
                pooled_cluster = int(fields[3])
                
                rep1_pooled_cluster = cluster_max_matches[0][rep1_cluster]
                rep2_pooled_cluster = cluster_max_matches[1][rep2_cluster]

                if rep1_pooled_cluster == rep2_pooled_cluster:
                    if pooled_cluster == rep1_pooled_cluster:
                        out.write(line)
                else:
                    # inconsistent
                    pass
                    #out.write("{}\t{}\t{}\t{}\n".format(region, rep1_cluster, rep2_cluster, "NA"))


    quit()
                    
    return cluster_repr_mat_file


def get_consistent_dpgp_trajectories_v2(
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
    
    # NEW PART (for v2): go through each region and do matching

    # first get cluster means for pooled DP-GP
    pooled_cluster_means = get_cluster_means(pooled_cluster_file, pooled_timeseries_file)


    # now load up rep1 rep2 timeseries files and out file
    cluster_repr_mat_file = "{}/{}.clusterings.reproducible.test.txt.gz".format(out_dir, prefix)
    #if not os.path.isfile(cluster_repr_mat_file):
    if True:
        with gzip.open(rep1_timeseries_file, 'r') as rep1:
            with gzip.open(rep2_timeseries_file, 'r') as rep2:
                with gzip.open(cluster_repr_mat_file, 'w') as out:

                    current_region_idx = 0
                    
                    while True:

                        rep1_line = rep1.readline()
                        rep2_line = rep2.readline()
                        
                        rep1_fields = rep1_line.strip().split("\t")
                        rep2_fields = rep2_line.strip().split("\t")
                        
                        # break if end of file
                        if rep1_line == "":
                            break

                        # ignore header line
                        if rep1_fields[0].startswith("d"):
                            continue

                        # split out region and timepoints
                        rep1_region, rep1_timepoints = rep1_fields[0], np.array(rep1_fields[1:]).reshape((1, len(rep1_fields[1:]))) # might want to change shape
                        rep2_region, rep2_timepoints = rep2_fields[0], np.array(rep2_fields[1:]).reshape((1, len(rep2_fields[1:])))
                        
                        assert rep1_region == rep2_region
                        
                        # compare to means
                        # TODO but give some softness to it - soft clustering
                        
                        rep1_corr = np.argmax(
                            get_corr_mat(rep1_timepoints, pooled_cluster_means),
                            axis=1)
                        rep2_corr = np.argmax(
                            get_corr_mat(rep2_timepoints, pooled_cluster_means),
                            axis=1)

                        if rep1_corr == rep2_corr:
                            out.write("{}\n".format(rep1_region))  # TODO add in pooled cluster

                        current_region_idx += 1

                        if current_region_idx % 1000 == 0:
                            print current_region_idx

                        #import ipdb
                        #ipdb.set_trace()

                        # arg max
                        
                        
                        # compare

                        # keep if consistent

    
    #import ipdb
    #ipdb.set_trace()

    

    

        
    quit()
    
    # compute correlations, rep1-pooled and rep2-pooled
    true_reps = ["rep1", "rep2"]
    cluster_corrs = []
    cluster_max_matches = []
    for rep_idx in range(len(true_reps)):
        cluster_corrs.append(get_corr_mat(
            cluster_means[rep_idx], cluster_means[2]))
        cluster_max_matches.append(np.argmax(cluster_corrs[rep_idx], axis=1) + 1)

                    
    return cluster_repr_mat_file



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
    
        # for each cluster (from POOLED clusters), get the mean and covariance matrix
        cluster_means, cluster_covariances = get_cluster_sufficient_stats(
            pooled_cluster_file, pooled_timeseries_file)

        # for rep1, rep2, pooled:
        # Use the PDF set at (0.0573) which is the 95% confidence interval to see if it falls in cluster
        region_idx = 0
        confidence_interval = 0.95
        corr_cutoff = 0.05
        pdf_cutoff = chi2.pdf(1-confidence_interval, cluster_means[0].shape[0]) # the pdf val of a multivariate normal is the chi2 value at alpha=pval
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
                       
                        # split out region and timepoints
                        #rep1_region, rep1_timepoints = rep1_fields[0], np.array(rep1_fields[1:]).reshape((1, len(rep1_fields[1:]))) # might want to change shape
                        #rep2_region, rep2_timepoints = rep2_fields[0], np.array(rep2_fields[1:]).reshape((1, len(rep2_fields[1:])))
                        #pooled_region, pooled_timepoints = pooled_fields[0], np.array(pooled_fields[1:]).reshape((1, len(rep2_fields[1:])))

                        rep1_region, rep1_timepoints = rep1_fields[0], zscore(np.array(map(float, rep1_fields[1:])))
                        rep2_region, rep2_timepoints = rep2_fields[0], zscore(np.array(map(float, rep2_fields[1:])))
                        pooled_region, pooled_timepoints = pooled_fields[0], zscore(np.array(map(float, pooled_fields[1:])))

                        assert rep1_region == rep2_region
                        assert rep1_region == pooled_region
                        assert rep2_region == pooled_region

                        # now for each cluster, check consistency
                        consistent_clusters = get_consistent_soft_clusters(
                            cluster_means,
                            cluster_covariances,
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
    
    # sanity check - replot the soft clusters
    plot_dir = "{}/soft/plot".format(out_dir)
    if not os.path.isdir(plot_dir):
        os.system("mkdir -p {}".format(plot_dir))
        
        # use R to plot 
        plot_soft_clusters = ("plot_soft_trajectories.R {0}/soft/plot {0}/soft/{1}").format(
            out_dir, "*.gz")
        print plot_soft_clusters
        os.system(plot_soft_clusters)
    

    return
