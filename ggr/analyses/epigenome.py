# analyses related to epigenomic datasets

import os
import gzip

import numpy as np
import pandas as pd

from ggr.util.utils import run_shell_cmd

from ggr.analyses.bioinformatics import run_homer

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
        stable_out_file,
        min_group_size=300):
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

    # remove the small H3K27me3 groups
    dynamic_histones = dynamic_histones[~(
        ((dynamic_histones["H3K27me3"] != 9) & (dynamic_histones["H3K27ac"] != 9)) &
        ((dynamic_histones["H3K27me3"] != 9) & (dynamic_histones["H3K4me1"] != 9))
    )]

    # remove small groups
    dynamic_histones = dynamic_histones.groupby("histone_cluster").filter(lambda x: len(x) > min_group_size)
    
    dynamic_histones.to_csv(dynamic_out_file, sep='\t', header=True, index=False)
    
    # and extract stable
    stable_histones = histone_overlaps_marked[(
        ((histone_overlaps["H3K27ac"] == 9) | (histone_overlaps["H3K27ac"] == 4)) &
        ((histone_overlaps["H3K4me1"] == 9) | (histone_overlaps["H3K4me1"] == 4)) &
        ((histone_overlaps["H3K27me3"] == 9) | (histone_overlaps["H3K27me3"] == 4))
    )]
    stable_histones.to_csv(stable_out_file, sep='\t', header=True, index=False)

    return


def run_bioinformatics_on_bed(bed_dir, background_bed, regex, out_dir):
    """Given a directory of BED files, run HOMER/GREAT
    """
    # get the files
    bed_files = glob.glob("{}/*{}".format(bed_dir, regex))

    # run GREAT
    great_dir = "{}/great".format(out_dir)
    run_shell_cmd("mkdir -p {}".format(great_dir))
    for bed_file in bed_files:
        prefix = os.path.basename(bed_file).split(regex)[0]
        if not os.path.isdir("{}/{}".format(great_dir, prefix)):
            run_shell_cmd("mkdir -p {}/{}".format(great_dir, prefix))
            # TODO put in background set
            run_great = "run_rgreat.R {0} {1}/{2}".format(
                bed_file,
                great_dir,
                prefix)

    # run HOMER
    homer_dir = "{}/homer".format(out_dir)
    run_shell_cmd("mkdir -p {}".format(homer_dir))
    for bed_file in bed_files:
        prefix = os.path.basename(bed_file).split(regex)[0]
        if not os.path.isdir("{}/{}".format(homer_dir, prefix)):
            run_shell_cmd("mkdir -p {}/{}".format(homer_dir, prefix))
            run_homer(
                bed_file,
                background_bed,
                "{}/{}".format(homer_dir, prefix))

    return None


def choose_best_summit(intersect_output_file, out_file):
    """Given the results of an intersect (with -wao option)
    reduce to the best summit, and output a summit file
    """
    with open(intersect_output_file, "r") as fp:
        with open(out_file, "w") as out:

            current_region = ""
            current_qval = 0.0
            current_summit_pos = 0
            for line in fp:
                # compare to current region
                fields = line.strip().split("\t")
                region = "{}-{}-{}".format(fields[0], fields[1], fields[2])
                
                # if same region, check if pval stronger. if so, change summit
                if region == current_region:
                    qval = float(fields[12])
                    if qval > current_qval:
                        current_qval = qval
                        current_summit_pos = int(fields[13])
                else:
                    if current_region != "":
                        # write out current region
                        region_fields = current_region.split("-")
                        start = int(region_fields[1])
                        stop = int(region_fields[2])
                        out.write("{}\t{}\t{}\n".format(
                            region_fields[0],
                            start + current_summit_pos,
                            start + current_summit_pos + 1))

                    # save in new region
                    current_region = region
                    current_qval = float(fields[12])
                    current_summit_pos = int(fields[13])

            # finally write out the last bit
            region_fields = current_region.split("-")
            start = int(region_fields[1])
            stop = int(region_fields[2])
            out.write("{}\t{}\t{}\n".format(
                region_fields[0],
                start + current_summit_pos,
                start + current_summit_pos + 1))
                    
    return None


def get_best_summit(master_bed, peak_files, out_file):
    """From the peak files, choose the best summit per master bed file
    """
    print master_bed
    
    # intersect across all relevant peak files
    intersect = "bedtools intersect -wao -a {0} -b {1} > summit_overlap.tmp".format(
        master_bed, " ".join(peak_files))
    run_shell_cmd(intersect)

    # and then reduce to best summit
    choose_best_summit("summit_overlap.tmp", out_file)

    # rm tmp file
    run_shell_cmd("rm summit_overlap.tmp")
    
    return


def convert_overlaps_bed_to_id_mappings(overlaps_bed, extend_len=0):
    """convert overlaps bed file to id mappings
    """
    overlaps = pd.read_csv(overlaps_bed, sep="\t", header=None)
    overlaps = overlaps[overlaps.iloc[:,4] != -1]
    overlaps["atac_id"] = overlaps.iloc[:,0].astype(str) + ":" + overlaps.iloc[:,1].astype(str) + "-" + overlaps.iloc[:,2].astype(str)
    overlaps["histone_id"] = overlaps.iloc[:,3].astype(str) + ":" + (overlaps.iloc[:,4]+extend_len).astype(str) + "-" + (overlaps.iloc[:,5]-extend_len).astype(str)
    overlaps = overlaps[["atac_id", "histone_id"]]
    
    return overlaps


def get_aggregate_chromatin_state_summary(
        regions,
        trajectories,
        traj_idx,
        atac_mat,
        histones,
        index=None):
    """aggregate information
    """
    num_regions = float(regions.shape[0])

    # set up the trajectory info
    summary = pd.DataFrame(
        data=np.zeros((len(trajectories))),
        index=trajectories)
    summary.loc["TRAJ.{}".format(traj_idx),0] = 1
    summary = summary.transpose()
    summary.insert(0, "num_regions", num_regions) # TODO fix this for num bp covered!
    
    # first extract atac
    group_atac = atac_mat[atac_mat.index.isin(regions)]
    group_atac = group_atac.sum(axis=0) / num_regions # TODO some normalization first?
    group_atac = pd.DataFrame(group_atac).transpose()
    group_atac.columns = ["ATAC.{}".format(val) for val in group_atac.columns]
    summary = pd.concat([summary, group_atac], axis=1)
    
    # then for each histone, extract info
    for histone, mapping, histone_mat in histones:
        histone_regions = mapping[mapping["atac_id"].isin(regions)]["histone_id"].values
        group_histone = histone_mat[histone_mat.index.isin(histone_regions)]
        group_histone = group_histone.sum(axis=0) / num_regions # TODO some normalization first?
        group_histone = group_histone.fillna(0)
        group_histone = pd.DataFrame(group_histone).transpose()
        group_histone.columns = ["{}.{}".format(histone, val) for val in group_histone.columns]
        summary = pd.concat([summary, group_histone], axis=1)

    if index is not None:
        summary.index = [index]
        
    return summary
