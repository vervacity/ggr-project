"""Code to integrate atac w histone info
"""

import os
import glob
import gzip
import logging

import pandas as pd
import numpy as np

from ggr.util.bioinformatics import make_deeptools_heatmap
from ggr.util.bed_utils import id_to_bed


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
        print intersect_bed
        os.system(intersect_bed)

        # and then adjust to get nearest histone
        out_nearest_tmp = "{}/{}.overlap.{}.nearest.tmp.txt.gz".format(
            out_dir, prefix, histone)
        get_nearest_mark(
            out_tmp_file,
            out_nearest_tmp,
            histone_assignment)
    
        # and take these and overlap with each other to get out file
        mark_clusters = pd.read_table(out_nearest_tmp, header=None, sep='\t')
        mark_clusters.columns = ["region", histone]

        if i == 0:
            master = pd.DataFrame(mark_clusters)
        else:
            master = master.merge(mark_clusters, on="region")

    
    # sort and then mark by overall cluster
    master.index = master["region"]
    del master["region"]

    master["histone_cluster"] = [0 for i in range(master.shape[0])]
    num_cluster_cols = master.shape[1] - 1
    for i in range(num_cluster_cols):
        master["histone_cluster"] = master["histone_cluster"] + (10**i) * (master.iloc[:,num_cluster_cols-i-1])
    master["histone_cluster"] = master["histone_cluster"].astype(int)

    master_sorted = master.sort_values("histone_cluster")
    
    # save
    master_sorted.to_csv(out_file, sep='\t', compression="gzip")
            
    return


def make_deeptools_bed(bed_file, out_bed_file, cluster_col=3, midpoint=False):
    """Given a bed file with cluster column, put # in appropriate places
    and adjust to midpoint if necessary
    """
    current_cluster = ""
    with gzip.open(out_bed_file, 'w') as out:
        with gzip.open(bed_file, 'r') as fp:
            for line in fp:
                cluster = line.strip().split('\t')[cluster_col]

                if midpoint:
                    fields = line.strip().split('\t')
                    midpoint = (int(fields[1]) + int(fields[2])) / 2
                    fields[1], fields[2] = str(midpoint), str(midpoint + 1)
                    line = "{}\n".format("\t".join(fields))
                    
                if current_cluster == cluster:
                    out.write(line)
                else:
                    current_cluster = cluster
                    out.write("#\n")
                    out.write(line)
    
    return


def order_and_viz_dynamic_epigenome(
        args,
        out_dir,
        prefix,
        activating_histones,
        activating_histone_files,
        all_histones,
        histone_assignment="nearest"):
    """Order clusters by hclust and histones more simple
    """
    # get overlap with histones
    dynamic_histone_overlap_dir = "{}/histone_overlap".format(out_dir)
    os.system("mkdir -p {}".format(dynamic_histone_overlap_dir))

    args.integrative["atac_dynamic_w_histones_mat"] = (
        "{0}/{1}.overlaps.mat.txt.gz").format(
            dynamic_histone_overlap_dir,
            prefix)
    if not os.path.isfile(args.integrative["atac_dynamic_w_histones_mat"]):
        get_histone_overlaps(
            args.atac["counts_pooled_rlog_dynamic_bed"],
            activating_histone_files,
            dynamic_histone_overlap_dir,
            args.integrative["atac_dynamic_w_histones_mat"],
            args.annot["chromsizes"],
            histone_assignment=histone_assignment)
        
    # TODO: get atac clusters and sort by (hard) cluster
    # Now merge in histone clusters, and sort by ATAC + histone
    # this ordering is the FINAL ordering
    dynamic_plots_dir = "{}/plots".format(out_dir)
    os.system("mkdir -p {}".format(dynamic_plots_dir))
    dynamic_atac_plot_prefix = "{}/{}.atac.epigenome_ordered".format(
        dynamic_plots_dir, prefix)
    dynamic_atac_subsample_file = "{}/{}.atac.ordered.subsample.txt.gz".format(
        dynamic_plots_dir, prefix)

    if not os.path.isfile(dynamic_atac_subsample_file):
    
        order_and_plot_atac = (
            "plot_dynamic_atac.R {0} {1} {2} {3}").format(
                args.atac["final_hard_clusters"],
                args.integrative["atac_dynamic_w_histones_mat"],
                dynamic_atac_plot_prefix,
                dynamic_atac_subsample_file)
        print order_and_plot_atac
        os.system(order_and_plot_atac)
        
    # take subsample file and make BED
    dynamic_atac_subsample_bed = "{}.bed.gz".format(
        dynamic_atac_subsample_file.split(".txt")[0])
    if not os.path.isfile(dynamic_atac_subsample_bed):
        make_bed = (
            "zcat {0} | "
            "grep -v region | "
            "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
            "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
            "gzip -c > {1}").format(
                dynamic_atac_subsample_file,
                dynamic_atac_subsample_bed)
        print make_bed
        os.system(make_bed)

    # and then run through deeptools for histone marks
    histone_colors = args.chipseq["histones"]["ordered_deeptools_colors"]
    for histone_idx in range(len(all_histones)):
        histone = all_histones[histone_idx]
        histone_color = histone_colors[histone_idx]
        histone_bigwigs = sorted(
            glob.glob("{}/{}".format(
                args.chipseq["data_dir"],
                args.chipseq["histones"][histone]["pooled_bigwig_glob"])))

        out_prefix = "{}/{}.{}_overlap".format(
            dynamic_plots_dir,
            prefix,
            histone)
        out_file = "{}.heatmap.profile.png".format(out_prefix)

        if not os.path.isfile(out_file):
            make_deeptools_heatmap(
                dynamic_atac_subsample_bed,
                histone_bigwigs,
                out_prefix,
                sort=False,
                referencepoint="center",
                color=histone_color)
            
    return args


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
    os.system("mkdir -p {0} {0}/ids {0}/bed {1} {1}/ids {1}/bed".format(
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
    dynamic_histones.to_csv(dynamic_out_file, sep='\t', header=True, index=False, compression="gzip")
    

    # and extract stable
    stable_histones = histone_overlaps_marked[(
        ((histone_overlaps["H3K27ac"] == 9) | (histone_overlaps["H3K27ac"] == 4)) &
        ((histone_overlaps["H3K4me1"] == 9) | (histone_overlaps["H3K4me1"] == 4)) &
        ((histone_overlaps["H3K27me3"] == 9) | (histone_overlaps["H3K27me3"] == 4))
    )]
    stable_histones.to_csv(stable_out_file, sep='\t', header=True, index=False, compression="gzip")

    return


def order_and_viz_stable_epigenome(
        args,
        out_dir,
        prefix,
        histones,
        histone_files,
        histone_assignment="nearest"):
    """Analyze stable epigenome
    """
    # get histone overlap with stable ATAC
    stable_histone_overlap_dir = "{}/histone_overlap".format(out_dir)
    os.system("mkdir -p {}".format(stable_histone_overlap_dir))
    
    args.integrative["atac_stable_w_histones_mat"] = (
        "{0}/{1}.overlaps.mat.txt.gz").format(
            stable_histone_overlap_dir,
            prefix)

    if not os.path.isfile(args.integrative["atac_stable_w_histones_mat"]):
        get_histone_overlaps(
            args.atac["counts_pooled_rlog_stable_bed"],
            histone_files,
            stable_histone_overlap_dir,
            args.integrative["atac_stable_w_histones_mat"],
            args.annot["chromsizes"],
            histone_assignment=histone_assignment)

    # first separate out dynamic and non dynamic chromatin mark data
    args.integrative["atac_stable_w_histones_dynamic_mat"] = (
        "{0}/{1}.overlaps.histones_dynamic.mat.txt.gz").format(
            stable_histone_overlap_dir,
            prefix)
    args.integrative["atac_stable_w_histones_stable_mat"] = (
        "{0}/{1}.overlaps.histones_stable.mat.txt.gz").format(
            stable_histone_overlap_dir,
            prefix)
    if not os.path.isfile(args.integrative["atac_stable_w_histones_stable_mat"]):
        split_stable_atac_by_dynamic_marks(
            args.integrative["atac_stable_w_histones_mat"],
            args.integrative["atac_stable_w_histones_dynamic_mat"],
            args.integrative["atac_stable_w_histones_stable_mat"])

    # get the BED file for each of these subsets
    args.integrative["atac_stable_w_histones_dynamic_bed"] = "{}.bed.gz".format(
        args.integrative["atac_stable_w_histones_dynamic_mat"].split(".mat")[0])
    if not os.path.isfile(args.integrative["atac_stable_w_histones_dynamic_bed"]):
        id_to_bed(
            args.integrative["atac_stable_w_histones_dynamic_mat"],
            args.integrative["atac_stable_w_histones_dynamic_bed"])

    args.integrative["atac_stable_w_histones_stable_bed"] = "{}.bed.gz".format(
        args.integrative["atac_stable_w_histones_stable_mat"].split(".mat")[0])
    if not os.path.isfile(args.integrative["atac_stable_w_histones_stable_bed"]):
        id_to_bed(
            args.integrative["atac_stable_w_histones_stable_mat"],
            args.integrative["atac_stable_w_histones_stable_bed"])


    # plot
    stable_plots_dir = "{}/plots".format(out_dir)
    os.system("mkdir -p {}".format(stable_plots_dir))

    bed_to_plot = [
        args.integrative["atac_stable_w_histones_dynamic_bed"],
        args.integrative["atac_stable_w_histones_stable_bed"]
    ]
        
    for bed_file in bed_to_plot:

        prefix = os.path.basename(bed_file).split(".bed")[0]
        
        # TODO go into R, make fake heatmap to get color bars for histone marks
        histone_colors = ["Reds", "Blues", "Greens"]
        for histone_idx in range(len(histones)):
            
            histone = histones[histone_idx]
            histone_color = histone_colors[histone_idx]
            
            histone_bigwigs = sorted(
                glob.glob("{}/{}".format(
                    args.chipseq["data_dir"],
                    args.chipseq["histones"][histone]["pooled_bigwig_glob"])))

            out_prefix = "{}/{}.{}_overlap".format(
                stable_plots_dir, prefix, histone)
            out_file = "{}.heatmap.profile.png".format(out_prefix)

            if not os.path.isfile(out_file):
                make_deeptools_heatmap(
                    bed_file,
                    histone_bigwigs,
                    out_prefix,
                    kval=0,
                    referencepoint="center",
                    color=histone_color)

    return args


def run(args):
    """Pipeline for integrating epigenome datasets
    Assumes that ATAC analysis and ChIP analysis has
    already been run before this
    """
    args.integrative = {}

    # set up histones
    histones = args.chipseq["histones"]["ordered_names"]
    histone_files = [
        (histone,
         args.chipseq["histones"][histone]["clusters_bed"],
         args.params["histones"][histone]["overlap_extend_len"])
        for histone in histones]

    # dynamic ATAC
    logging.info("EPIGENOME: DYNAMIC: integrate dynamic ATAC w histones")
    dynamic_epigenome_prefix = "ggr.epigenome.dynamic"

    activating_marks = histones[0:2]
    activating_mark_files = histone_files[0:2]

    args = order_and_viz_dynamic_epigenome(
        args,
        args.folders["epigenome_dynamic_dir"],
        dynamic_epigenome_prefix,
        activating_marks,
        activating_mark_files,
        histones,
        histone_assignment=args.params["histone_overlap_assignment"])

    dynamic_state_beds = glob.glob("{}/*bed.gz".format(
        args.folders["epigenome_dynamic_state_bed_dir"]))
    if len(dynamic_state_beds) == 0:
        cluster_by_chromatin_state(
            args.atac["final_soft_clusters"],
            args.integrative["atac_dynamic_w_histones_mat"],
            activating_marks,
            args.folders["epigenome_dynamic_clusters_dir"],
            dynamic_epigenome_prefix,
            min_region_num=args.params["chrom_state_cluster_min_size"])
    
    # stable ATAC
    logging.info("EPIGENOME: STABLE: integrate dynamic ATAC w histones")
    stable_epigenome_prefix = "ggr.epigenome.stable"
    args = order_and_viz_stable_epigenome(
        args,
        args.folders["epigenome_stable_dir"],
        stable_epigenome_prefix,
        histones,
        histone_files,
        histone_assignment=args.params["histone_overlap_assignment"])
    
    # cluster_by_chromatin_state
    stable_state_beds = glob.glob("{}/*bed.gz".format(
        args.folders["epigenome_stable_state_bed_dir"]))
    if len(stable_state_beds) == 0:
        cluster_by_chromatin_state(
            [args.atac["counts_pooled_rlog_stable"]],
            args.integrative["atac_stable_w_histones_mat"],
            histones,
            args.folders["epigenome_stable_clusters_dir"],
            stable_epigenome_prefix,
            min_region_num=args.params["chrom_state_cluster_min_size"])

    # TODO now go through all folders and run GREAT and HOMER
    # key groups: ATAC final trajectorie (no histone info), dynamic atac (mark/state), stable atac (mark/state)
    # 5 folders

    label_dir_handles = [
        "atac_dp-gp_final_bed_dir",
        "epigenome_dynamic_mark_bed_dir",
        "epigenome_dynamic_state_bed_dir",
        "epigenome_stable_mark_bed_dir",
        "epigenome_stable_state_bed_dir"
    ]
    
    for label_dir_handle in label_dir_handles:
        label_bed_files = glob.glob("{}/*bed.gz".format(args.folders[label_dir_handle]))
        great_dir = "{}/great".format(args.folders[label_dir_handle].split("/bed")[0])
        os.system("mkdir -p {}".format(great_dir))
        
        for label_bed_file in label_bed_files:
            label_bed_prefix = os.path.basename(label_bed_file).split(".bed")[0]
            
            # run GREAT
            great_files = glob.glob("{0}/{1}*".format(great_dir, label_bed_prefix))
            run_great = "run_rgreat.R {0} {1}/{2}".format(
                label_bed_file,
                great_dir,
                label_bed_prefix)
            print run_great
            os.system(run_great)
        
        # run HOMER
        
        




    quit()
    
    # and generate group BED files
    stable_clusters = pd.read_table(
        args.integrative["atac_stable_histone_ordered_bed"], header=None)
    stable_clusters.columns = ["chrom", "start", "stop", "cluster"]
    cluster_ids = list(set(stable_clusters["cluster"].tolist()))

    for cluster_id in cluster_ids:
        cluster_subset = stable_clusters[stable_clusters["cluster"] == cluster_id]
        
        if cluster_subset.shape[0] > 1000: # manually chosen
            cluster_prefix = "ggr.atac_stable.cluster_{0}".format(cluster_id)
            
            cluster_file = "{0}/{1}.bed.gz".format(
                args.folders["epigenome_stable_dir"],
                cluster_prefix)

            if not os.path.isfile(cluster_file):
                cluster_subset.to_csv(
                    cluster_file,
                    sep='\t', header=False, index=False,
                    compression="gzip")
    
            # run GREAT
            great_files = glob.glob("{}/{}*".format(
                args.folders["epigenome_stable_great_dir"],
                cluster_prefix))
            if len(great_files) == 0:
                run_great = "run_rgreat.R {0} {1}/{2}".format(
                    cluster_file,
                    args.folders["epigenome_stable_great_dir"],
                    cluster_prefix)
                print run_great
                os.system(run_great)

            # run HOMER
            

    


            
    

    return args
