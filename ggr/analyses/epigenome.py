# analyses related to epigenomic datasets

import os
import gzip

import pandas as pd

from ggr.util.utils import run_shell_cmd


def get_nearest_mark(intersect_bed_file, out_file, histone_assignment):
    """Assumes that you did an intersect 
    where you have two regions on a line
    """
    current_region_info = [""] # order: name, dict of distance:cluster TODO(dk) differential:cluster
    with gzip.open(out_file, 'w') as out:
        with gzip.open(intersect_bed_file, 'r') as fp:
            for line in fp:
                
                fields = line.strip().split('\t')
                region = "{}:{}-{}".format(fields[0], fields[1], fields[2])
                cluster = fields[6]

                if current_region_info[0] == region:
                    # same region, save in the following info
                    distance = (int(fields[5]) - int(fields[4])) - (int(fields[2]) - int(fields[1]))
                    current_region_info[1][distance] = cluster
                    current_region_info[2].append(int(cluster))
                else:
                    # new region: save out old region, and replace with new
                    if current_region_info[0] != "":
                        if histone_assignment == "most_diff":
                            most_diff = np.abs(np.array(current_region_info[2]) - 4.)
                            most_diff[most_diff > 4] = 0 # ignore 9s
                            most_diff_idx = np.argmax(most_diff)
                            out.write("{}\t{}\n".format(
                                current_region_info[0],
                                current_region_info[2][most_diff_idx]))
                        elif histone_assignment == "nearest":
                            closest_dist = min(current_region_info[1].keys())
                            out.write("{}\t{}\n".format(
                                current_region_info[0],
                                current_region_info[1][closest_dist]))
                        else:
                            raise Exception("unrecognized method for histone assignment!")
                            
                    # replace with new
                    if fields[6] == ".":
                        current_region_info = [region, {0: "9"}, [9]]
                    else:
                        distance = (int(fields[5]) - int(fields[4])) - (int(fields[2]) - int(fields[1]))
                        current_region_info = [region, {distance: cluster}, [int(cluster)]]

            # don't forget to output the last region
            closest_dist = min(current_region_info[1].keys())
            out.write("{}\t{}\n".format(
                current_region_info[0],
                current_region_info[1][closest_dist]))

    return


def get_histone_overlaps(
        master_bed_file,
        overlap_params,
        out_dir,
        out_file,
        chromsizes_file,
        histone_assignment="nearest"):
    """Given a master bed file and a list of histone files,
    get the overlaps and merge back into a file with all clusters
    """
    prefix = os.path.basename(master_bed_file).split(".bed")[0]
    
    for i in range(len(overlap_params)):
        histone, overlap_bed_file, search_dist = overlap_params[i]
        out_tmp_file = "{}/{}.overlap.{}.tmp.bed.gz".format(
            out_dir, prefix, histone)
        intersect_bed = (
            "bedtools slop -i {0} -g {1} -b {2} | "
            "bedtools intersect -wb -loj -a {3} -b stdin | "
            "gzip -c > {4}").format(
                overlap_bed_file,
                chromsizes_file,
                search_dist,
                master_bed_file,
                out_tmp_file)
        run_shell_cmd(intersect_bed)

        # and then adjust to get nearest histone
        out_nearest_tmp = "{}/{}.overlap.{}.nearest.tmp.txt.gz".format(
            out_dir, prefix, histone)
        get_nearest_mark(
            out_tmp_file,
            out_nearest_tmp,
            histone_assignment)
    
        # and take these and overlap with each other to get out file
        mark_clusters = pd.read_table(out_nearest_tmp, header=None, sep='\t')
        mark_clusters.columns = ["id", histone]

        if i == 0:
            master = pd.DataFrame(mark_clusters)
        else:
            master = master.merge(mark_clusters, on="id")
    
    # sort and then mark by overall cluster
    master.index = master["id"]
    del master["id"]

    master["histone_cluster"] = [0 for i in range(master.shape[0])]
    num_cluster_cols = master.shape[1] - 1
    for i in range(num_cluster_cols):
        master["histone_cluster"] = master["histone_cluster"] + (10**i) * (master.iloc[:,num_cluster_cols-i-1])
    master["histone_cluster"] = master["histone_cluster"].astype(int)

    master_sorted = master.sort_values("histone_cluster")
    
    # save
    master_sorted.to_csv(out_file, sep='\t', compression="gzip")
            
    return


def cluster_by_chromatin_marks(
        cluster_mat_file,
        histone_names,
        state_name,
        mark_dir,
        state_dir,
        prefix,
        min_region_num=100,
        ignore_mark_levels=[9],
        ignore_state_levels=[99, 999]):
    """Take in a list of cluster files and sort by histone marks
    """
    assert os.path.isdir(mark_dir)
    assert os.path.isdir(state_dir)

    # read in mat file
    clusters = pd.read_table(cluster_mat_file)

    # start with high level clusters
    atac_clusters = list(set(clusters["cluster"].tolist()))

    for atac_cluster in atac_clusters:

        atac_cluster_prefix = "atac-cluster_{}".format(atac_cluster)
        atac_cluster = clusters[clusters["cluster"] == atac_cluster]

        # for marks:
        for histone in histone_names:
            
            # first get marks for that cluster
            mark_levels = list(set(atac_cluster[histone].tolist()))

            for mark_level in mark_levels:
                if mark_level in ignore_mark_levels:
                    continue

                # get the ids associated
                ids = atac_cluster[atac_cluster[histone] == mark_level]
                if ids.shape[0] > min_region_num:
                    # write out
                    cluster_file = "{}/{}.{}.{}-cluster_{}.txt.gz".format(
                        mark_dir, prefix, atac_cluster_prefix, histone, mark_level)
                    ids.to_csv(cluster_file, columns=["id"], sep="\t", header=False, index=False, compression="gzip")

        # then for states:
        state_levels = list(set(atac_cluster[state_name].tolist()))
        for state_level in state_levels:
            if state_level in ignore_state_levels:
                continue

            # get the ids associated
            ids = atac_cluster[atac_cluster[state_name] == state_level]
            if ids.shape[0] > min_region_num:
                # write out
                cluster_file = "{}/{}.{}.state-cluster_{}.txt.gz".format(
                    state_dir, prefix, atac_cluster_prefix, state_level)
                ids.to_csv(cluster_file, columns=["id"], sep="\t", header=False, index=False, compression="gzip")
    
    return None


def split_stable_atac_by_dynamic_marks(
        overlaps_file,
        dynamic_out_file,
        stable_out_file):
    """take the stable overlaps file and split out stable marks
    from dynamic marks (near the ATAC peak)
    """
    histone_overlaps = pd.read_table(overlaps_file)

    # first remove 999
    histone_overlaps_marked = histone_overlaps[
        histone_overlaps["histone_cluster"] != 999]

    # extract dynamic
    dynamic_histones = histone_overlaps_marked[~(
        ((histone_overlaps["H3K27ac"] == 9) | (histone_overlaps["H3K27ac"] == 4)) &
        ((histone_overlaps["H3K4me1"] == 9) | (histone_overlaps["H3K4me1"] == 4)) &
        ((histone_overlaps["H3K27me3"] == 9) | (histone_overlaps["H3K27me3"] == 4))
    )]
    dynamic_histones.to_csv(dynamic_out_file, sep='\t', header=True, index=False)
    print dynamic_histones.shape

    # and extract stable
    stable_histones = histone_overlaps_marked[(
        ((histone_overlaps["H3K27ac"] == 9) | (histone_overlaps["H3K27ac"] == 4)) &
        ((histone_overlaps["H3K4me1"] == 9) | (histone_overlaps["H3K4me1"] == 4)) &
        ((histone_overlaps["H3K27me3"] == 9) | (histone_overlaps["H3K27me3"] == 4))
    )]
    stable_histones.to_csv(stable_out_file, sep='\t', header=True, index=False)

    return


# this is old! throw out soon
def cluster_by_chromatin_state(
        soft_cluster_files,
        histone_overlap_mat_file,
        histone_names,
        out_dir,
        prefix,
        min_region_num=100):
    """Take in a list of cluster files and sort by histone marks
    """
    # set up directories for each type of split
    mark_dir = "{}/by_mark".format(out_dir)
    state_dir = "{}/by_state".format(out_dir)
    run_shell_cmd("mkdir -p {0} {0}/ids {0}/bed {1} {1}/ids {1}/bed".format(
        mark_dir, state_dir))

    for soft_cluster_file in soft_cluster_files:

        if "cluster" in soft_cluster_file:
            soft_prefix = "{}.{}".format(
                prefix,
                os.path.basename(soft_cluster_file).split('.')[2])
        else:
            soft_prefix = prefix
            
        # overlap files
        soft_cluster_data = pd.read_table(soft_cluster_file)
        soft_cluster_data = pd.DataFrame(soft_cluster_data.iloc[:,0])
        soft_cluster_data.columns = ["region"]
        histone_overlap_data = pd.read_table(histone_overlap_mat_file)

        cluster_w_histones = soft_cluster_data.merge(
            histone_overlap_data, how="left", on="region")
        
        # histone marks alone (dynamics)
        # TODO don't take 9s.
        for histone in histone_names:
            cluster_w_histones[histone] = cluster_w_histones[histone].astype(int)
            mark_states = list(set(cluster_w_histones[histone].tolist()))

            for mark_state in mark_states:
                if mark_state == 9:
                    continue
                
                subset = cluster_w_histones[cluster_w_histones[histone] == mark_state]
                if subset.shape[0] > min_region_num:
                    overlap_file = "{}/{}.{}-{}.txt.gz".format(
                        "{}/ids".format(mark_dir), soft_prefix, histone, mark_state)
                    pd.DataFrame(subset["region"]).to_csv(
                        overlap_file, sep='\t', header=False, index=False, compression="gzip")
            
        # chromatin state dynamics
        # don't take 99s
        dynamic_states = list(set(cluster_w_histones["histone_cluster"].tolist()))
        for dynamic_state in dynamic_states:
            if dynamic_state == 99:
                continue

            if dynamic_state == 999:
                continue
            
            subset = cluster_w_histones[cluster_w_histones["histone_cluster"] == dynamic_state]
            if subset.shape[0] > min_region_num:
                overlap_file = "{}/{}.state-{}.txt.gz".format(
                    "{}/ids".format(state_dir), soft_prefix, dynamic_state)
                pd.DataFrame(subset["region"]).to_csv(
                    overlap_file, sep='\t', header=False, index=False, compression="gzip")

    # and convert all to BED
    mark_files = glob.glob("{}/ids/*.txt.gz".format(mark_dir))
    for mark_file in mark_files:
        mark_bed_file = "{0}/bed/{1}.bed.gz".format(mark_dir, os.path.basename(mark_file).split('.txt')[0])
        id_to_bed(mark_file, mark_bed_file, sort=True)

    state_files = glob.glob("{}/ids/*.txt.gz".format(state_dir))
    for state_file in state_files:
        state_bed_file = "{0}/bed/{1}.bed.gz".format(state_dir, os.path.basename(state_file).split('.txt')[0])
        id_to_bed(state_file, state_bed_file, sort=True)
    
    return None
