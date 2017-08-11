"""Code to integrate atac w histone info
"""

import os
import glob
import gzip

import pandas as pd

from ggr.util.bioinformatics import make_deeptools_heatmap


def get_nearest_mark(intersect_bed_file, out_file):
    """Assumes that you did an intersect 
    where you have two regions on a line
    """
    current_region_info = [""] # order: name, dict of distance:cluster
    with gzip.open(out_file, 'w') as out:
        with gzip.open(intersect_bed_file, 'r') as fp:
            for line in fp:
                
                fields = line.strip().split('\t')
                region = "{}:{}-{}".format(fields[0], fields[1], fields[2])

                if current_region_info[0] == region:
                    # same region, save in the following info
                    distance = (int(fields[5]) - int(fields[4])) - (int(fields[2]) - int(fields[1]))
                    current_region_info[1][distance] = fields[6]
                else:
                    # new region: save out old region, and replace with new
                    if current_region_info[0] != "":
                        closest_dist = min(current_region_info[1].keys())
                        out.write("{}\t{}\n".format(
                            current_region_info[0],
                            current_region_info[1][closest_dist]))
                        
                    # replace with new
                    if fields[6] == ".":
                        current_region_info = [region, {0: "9"}]
                    else:
                        distance = (int(fields[5]) - int(fields[4])) - (int(fields[2]) - int(fields[1]))
                        current_region_info = [region, {distance: fields[6]}]
    
    return


def get_histone_overlaps(
        master_bed_file,
        overlap_bed_files,
        out_dir,
        out_file,
        chromsizes_file):
    """Given a master bed file and a list of histone files,
    get the overlaps and merge back into a file with all clusters
    """
    prefix = os.path.basename(master_bed_file).split(".bed")[0]
    
    for i in range(len(overlap_bed_files)):
        histone, overlap_bed_file = overlap_bed_files[i]
        out_tmp_file = "{}/{}.overlap.{}.tmp.bed.gz".format(
            out_dir, prefix, histone)
        intersect_bed = (
            "bedtools slop -i {0} -g {1} -b {2} | "
            "bedtools intersect -wb -loj -a {3} -b stdin | "
            "gzip -c > {4}").format(
                overlap_bed_file,
                chromsizes_file,
                1000,
                master_bed_file,
                out_tmp_file)
        print intersect_bed
        os.system(intersect_bed)

        # and then adjust to get nearest histone
        out_nearest_tmp = "{}/{}.overlap.{}.nearest.tmp.txt.gz".format(
            out_dir, prefix, histone)
        get_nearest_mark(out_tmp_file, out_nearest_tmp)
    
        # and take these and overlap with each other to get out file
        mark_clusters = pd.read_table(out_nearest_tmp, header=None, sep='\t')
        mark_clusters.columns = ["regions", histone]

        if i == 0:
            master = pd.DataFrame(mark_clusters)
        else:
            master = master.merge(mark_clusters, on="regions")

    
    # sort and then mark by overall cluster
    master.index = master["regions"]
    del master["regions"]

    master["joint_cluster"] = [0 for i in range(master.shape[0])]
    num_cluster_cols = master.shape[1] - 1
    for i in range(num_cluster_cols):
        master["joint_cluster"] = master["joint_cluster"] + (10**i) * (master.iloc[:,num_cluster_cols-i-1])
    master["joint_cluster"] = master["joint_cluster"].astype(int)

    master_sorted = master.sort_values("joint_cluster")
    
    # save
    master_sorted.to_csv(out_file, sep='\t', compression="gzip")
            
    return


def make_deeptools_bed(bed_file, out_bed_file, cluster_col=3, midpoint=True):
    """Given a bed file with cluster column, put # in appropriate places
    and also adjust to point val
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




def run(args):
    """Pipeline for integrating epigenome datasets
    Assumes that ATAC analysis and ChIP analysis has
    already been run before this
    """
    args.integrative = {}

    histones = ["H3K27ac", "H3K4me1", "H3K27me3"] # keep ordered
    histone_files = [
        (histone, args.chipseq["histones"][histone]["clusters_bed"])
        for histone in histones]

    # get overlap with stable ATAC
    args.integrative["atac_stable_w_histones_mat"] = (
        "{}/ggr.atac_stable.histone_diff.overlaps.mat.txt.gz").format(
            args.folders["epigenome_dir"])

    if not os.path.isfile(args.integrative["atac_stable_w_histones_mat"]):
        get_histone_overlaps(
            args.atac["counts_pooled_rlog_stable_bed"],
            histone_files,
            args.folders["epigenome_dir"],
            args.integrative["atac_stable_w_histones_mat"],
            args.annot["chromsizes"])

    # make the BED file
    args.integrative["atac_stable_histone_ordered_bed"] = "{}.bed.gz".format(
        args.integrative["atac_stable_w_histones_mat"].split(".mat")[0])
    if not os.path.isfile(args.integrative["atac_stable_histone_ordered_bed"]):
        make_bed = (
            "zcat {0} | "
            "grep -v regions | "
            "awk -F '\t' '{{ print $1\"\t\"$5 }}' | "
            "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
            "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
            "gzip -c > {1}").format(
                args.integrative["atac_stable_w_histones_mat"],
                args.integrative["atac_stable_histone_ordered_bed"])
        print make_bed
        os.system(make_bed)

    # make the deeptools bed
    args.integrative["atac_stable_histone_deeptools_bed"] = "{}.deeptools.bed.gz".format(
        args.integrative["atac_stable_w_histones_mat"].split(".mat")[0])
    if not os.path.isfile(args.integrative["atac_stable_histone_deeptools_bed"]):
        make_deeptools_bed(
            args.integrative["atac_stable_histone_ordered_bed"],
            args.integrative["atac_stable_histone_deeptools_bed"],
            cluster_col=3)
        
    # And plot
    histone_colors = ["Reds", "Blues", "Greens"]
    for histone_idx in range(len(histones)):

        histone = histones[histone_idx]
        histone_color = histone_colors[histone_idx]
        
        histone_bigwigs = sorted(
            glob.glob("{}/{}".format(
                args.chipseq["data_dir"],
                args.chipseq["histones"][histone]["pooled_bigwig_glob"])))

        out_prefix = "{}/ggr.atac_stable.{}_overlap".format(
            args.folders["epigenome_dir"], histone)
        out_file = "{}.heatmap.profile.png".format(out_prefix)

        if not os.path.isfile(out_file):
            make_deeptools_heatmap(
                args.integrative["atac_stable_histone_deeptools_bed"],
                histone_bigwigs,
                out_prefix,
                sort=False,
                color=histone_color)

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
