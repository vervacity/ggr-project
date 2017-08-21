"""Code to integrate atac w histone info
"""

import os
import glob
import gzip
import logging

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

            # don't forget to output the last region
            closest_dist = min(current_region_info[1].keys())
            out.write("{}\t{}\n".format(
                current_region_info[0],
                current_region_info[1][closest_dist]))

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


def order_and_viz_dynamic_epigenome_v1(
        args,
        out_dir,
        prefix,
        histones,
        histone_files):
    """Order clusters by hclust and histones
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
            histone_files,
            dynamic_histone_overlap_dir,
            args.integrative["atac_dynamic_w_histones_mat"],
            args.annot["chromsizes"])
     
    # run hclust on enumerated soft clusters to get initial ATAC ordering
    dynamic_ordering_dir = "{}/ordering".format(out_dir)
    os.system("mkdir -p {}".format(dynamic_ordering_dir))
    
    args.integrative["atac_dynamic_init_ordered"] = (
        "{0}/{1}.ordered.mat.txt.gz").format(
            dynamic_ordering_dir,
            prefix)
    if not os.path.isfile(args.integrative["atac_dynamic_init_ordered"]):
        plot_ordered_atac = ("init_order_dynamic_atac.R {0} {1} {2} {3}").format(
            args.atac["counts_pooled_rlog_dynamic"],
            args.atac["soft_cluster_mat"],
            args.integrative["atac_dynamic_w_histones_mat"],
            "{0}/{1}".format(
                dynamic_ordering_dir,
                prefix))
        print plot_ordered_atac
        os.system(plot_ordered_atac)

    # make a bed file, keeping epigenome cluster
    args.integrative["atac_dynamic_init_ordered_bed"] = (
        "{0}/{1}.ordered.init.bed.gz").format(
            dynamic_ordering_dir,
            prefix)
    if not os.path.isfile(args.integrative["atac_dynamic_init_ordered_bed"]):
        make_bed = (
            "zcat {0} | "
            "grep -v d00 | "
            "awk -F '\t' '{{ print $1\"\t\"$14}}' | "
            "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
            "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
            "gzip -c > {1}").format(
                args.integrative["atac_dynamic_init_ordered"],
                args.integrative["atac_dynamic_init_ordered_bed"])
        print make_bed
        os.system(make_bed)
    
    # here, generate a deeptools bed file to run deeptools
    deeptools_dynamic_dir = "{}/deeptools".format(out_dir)
    os.system("mkdir -p {}".format(deeptools_dynamic_dir))
    
    args.integrative["atac_dynamic_init_ordered_deeptools_bed"] = (
        "{0}/{1}.deeptools.bed.gz").format(
            deeptools_dynamic_dir,
            os.path.basename(
                args.integrative["atac_dynamic_init_ordered_bed"]).split(".bed")[0])
    if not os.path.isfile(args.integrative["atac_dynamic_init_ordered_deeptools_bed"]):
        make_deeptools_bed(
            args.integrative["atac_dynamic_init_ordered_bed"],
            args.integrative["atac_dynamic_init_ordered_deeptools_bed"],
            cluster_col=3)
        
    # run deeptools with ALL histones first to get a final ordering
    all_histone_bigwigs = []
    for histone in histones:
        all_histone_bigwigs = all_histone_bigwigs + sorted(
            glob.glob("{}/{}".format(
                args.chipseq["data_dir"],
                args.chipseq["histones"][histone]["pooled_bigwig_glob"])))

    args.integrative["atac_dynamic_final_ordered_bed"] = (
        "{0}/{1}.ordered.final.bed.gz").format(
            out_dir,
            prefix)
    point_mat_file = "{0}/{1}.deeptools.point.mat.tmp.gz".format(
        deeptools_dynamic_dir,
        prefix)
    if not os.path.isfile(args.integrative["atac_dynamic_final_ordered_bed"]):
        deeptools_compute_matrix = (
            "computeMatrix reference-point "
            "--referencePoint center "
            "-b 1000 -a 1000 "
            "-R {0} "
            "-S {1} "
            "--sortRegions=descend "
            "--outFileSortedRegions {2} "
            "-o {3} ").format(
                args.integrative["atac_dynamic_init_ordered_deeptools_bed"],
                " ".join(all_histone_bigwigs),
                args.integrative["atac_dynamic_final_ordered_bed"].split(".gz")[0],
                point_mat_file)
        print deeptools_compute_matrix
        os.system(deeptools_compute_matrix)
        os.system("gzip {}".format(
            args.integrative["atac_dynamic_final_ordered_bed"].split(".gz")[0]))

    # NOW can go back to ATAC and plot out with histogram
    dynamic_plotting_dir = "{}/plots".format(
        out_dir)
    os.system("mkdir -p {}".format(dynamic_plotting_dir))
    dynamic_atac_plot = "{}/{}.ordered.atac_plot.png".format(
        dynamic_plotting_dir,
        prefix)
    dynamic_atac_subsample_file = "{}/{}.ordered.subsample.txt.gz".format(
        dynamic_plotting_dir,
        prefix)
    if not os.path.isfile(dynamic_atac_plot):
        plot_dynamic_atac = (
            "plot_dynamic_atac.R {0} {1} {2} {3}").format(
                args.atac["counts_pooled_rlog_dynamic"],
                args.integrative["atac_dynamic_final_ordered_bed"],
                dynamic_atac_plot,
                dynamic_atac_subsample_file)
        print plot_dynamic_atac
        os.system(plot_dynamic_atac)

    # take subsample file and make BED
    dynamic_atac_subsample_bed = "{}.bed.gz".format(
        dynamic_atac_subsample_file.split(".txt")[0])
    if not os.path.isfile(dynamic_atac_subsample_bed):
        make_bed = (
            "zcat {0} | "
            "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
            "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
            "gzip -c > {1}").format(
                dynamic_atac_subsample_file,
                dynamic_atac_subsample_bed)
        print make_bed
        os.system(make_bed)

    # and then run through deeptools for histone marks
    histone_colors = ["Reds", "Blues", "Greens"]
    for histone_idx in range(len(histones)):

        histone = histones[histone_idx]
        histone_color = histone_colors[histone_idx]
        
        histone_bigwigs = sorted(
            glob.glob("{}/{}".format(
                args.chipseq["data_dir"],
                args.chipseq["histones"][histone]["pooled_bigwig_glob"])))

        out_prefix = "{}/{}.{}_overlap".format(
            dynamic_plotting_dir,
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

    return


def order_and_viz_dynamic_epigenome_v2(
        args,
        out_dir,
        prefix,
        histones,
        histone_files):
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
            histone_files,
            dynamic_histone_overlap_dir,
            args.integrative["atac_dynamic_w_histones_mat"],
            args.annot["chromsizes"])
        
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
    histone_colors = ["Reds", "Blues", "Greens"]
    for histone_idx in range(len(histones)):

        histone = histones[histone_idx]
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
        min_region_fraction=0.2):
    """Take in a list of cluster files and sort by histone marks
    """
    # set up directories for each type of split
    mark_dir = "{}/by_mark".format(out_dir)
    state_dir = "{}/by_state".format(out_dir)
    dynamic_mark_dir = "{}/by_dynamic_mark".format(out_dir)
    dynamic_state_dir = "{}/by_dynamic_state".format(out_dir)
    os.system("mkdir -p {} {} {} {}".format(
        mark_dir, state_dir, dynamic_mark_dir, dynamic_state_dir))

    for soft_cluster_file in soft_cluster_files:

        soft_prefix = "{}.{}".format(
            prefix,
            os.path.basename(soft_cluster_file).split('.')[2])
        
        # overlap files
        soft_cluster_data = pd.read_table(soft_cluster_file, header=None)
        soft_cluster_data = pd.DataFrame(soft_cluster_data[0])
        soft_cluster_data.columns = ["region"]
        histone_overlap_data = pd.read_table(histone_overlap_mat_file)

        cluster_w_histones = soft_cluster_data.merge(
            histone_overlap_data, how="left", on="region")

        # set up a min region num to remove really small sets
        min_region_num = int(min_region_fraction * cluster_w_histones.shape[0])
        
        # perform several overlap types:
        # each histone mark alone (presence)
        # TODO - dont take absence, its the exact opposite (negative set)
        for histone in histone_names:
            subset = cluster_w_histones[cluster_w_histones[histone] != 9]
            if subset.shape[0] > min_region_num:
                overlap_file = "{}/{}.{}-present.txt.gz".format(
                    mark_dir, soft_prefix, histone)
                pd.DataFrame(subset["region"]).to_csv(
                    overlap_file, sep='\t', header=False, index=False, compression="gzip")

            #subset = cluster_w_histones[cluster_w_histones[histone] == 9]
            #if subset.shape[0] > min_region_num:
            #    overlap_file = "{}/{}.{}-absent.txt.gz".format(
            #        mark_dir, soft_prefix, histone)
            #    pd.DataFrame(subset["region"]).to_csv(
            #        overlap_file, sep='\t', header=False, index=False, compression="gzip")
            
        
        # histone marks in combination (just presence)
        histone_a = "H3K27ac"
        histone_b = "H3K4me1"

        # both present
        subset = cluster_w_histones[
            (cluster_w_histones[histone_a] != 9) &
            (cluster_w_histones[histone_b] != 9)]
        if subset.shape[0] > min_region_num:
            overlap_file = "{}/{}.{}-present_{}-present.txt.gz".format(
                state_dir, soft_prefix, histone_a, histone_b)
            pd.DataFrame(subset["region"]).to_csv(
                overlap_file, sep='\t', header=False, index=False, compression="gzip")

        # a present, b missing
        subset = cluster_w_histones[
            (cluster_w_histones[histone_a] != 9) &
            (cluster_w_histones[histone_b] == 9)]
        if subset.shape[0] > min_region_num:
            overlap_file = "{}/{}.{}-present_{}-absent.txt.gz".format(
                state_dir, soft_prefix, histone_a, histone_b)
            pd.DataFrame(subset["region"]).to_csv(
                overlap_file, sep='\t', header=False, index=False, compression="gzip")

        # a missing, b present
        subset = cluster_w_histones[
            (cluster_w_histones[histone_a] == 9) &
            (cluster_w_histones[histone_b] != 9)]
        if subset.shape[0] > min_region_num:
            overlap_file = "{}/{}.{}-absent_{}-present.txt.gz".format(
                state_dir, soft_prefix, histone_a, histone_b)
            pd.DataFrame(subset["region"]).to_csv(
                overlap_file, sep='\t', header=False, index=False, compression="gzip")

        # a missing, b missing
        # TODO don't take these, these are negatives (for function)
        #subset = cluster_w_histones[
        #    (cluster_w_histones[histone_a] == 9) &
        #    (cluster_w_histones[histone_b] == 9)]
        #if subset.shape[0] > min_region_num:
        #    overlap_file = "{}/{}.{}-absent_{}-absent.txt.gz".format(
        #        state_dir, soft_prefix, histone_a, histone_b)
        #    pd.DataFrame(subset["region"]).to_csv(
        #        overlap_file, sep='\t', header=False, index=False, compression="gzip")

        # histone marks alone (dynamics)
        # TODO don't take 9s.
        # NOTE: as currently defined, not useful
        for histone in histone_names:
            mark_states = list(set(cluster_w_histones[histone].tolist()))

            for mark_state in mark_states:
                if mark_state == 9:
                    continue
                
                subset = cluster_w_histones[cluster_w_histones[histone] == mark_state]
                #if subset.shape[0] > min_region_num:
                if subset.shape[0] > 100:
                    overlap_file = "{}/{}.{}-{}.txt.gz".format(
                        dynamic_mark_dir, soft_prefix, histone, mark_state)
                    pd.DataFrame(subset["region"]).to_csv(
                        overlap_file, sep='\t', header=False, index=False, compression="gzip")
            
        # chromatin state dynamics
        # don't take 99s
        # TODO: ignore 4s, and 9s? these are subsets of presence/absence?
        # NOTE: as currently defined, not useful
        dynamic_states = list(set(cluster_w_histones["histone_cluster"].tolist()))
        for dynamic_state in dynamic_states:
            if dynamic_state == 99:
                continue
            
            subset = cluster_w_histones[cluster_w_histones["histone_cluster"] == dynamic_state]
            #if subset.shape[0] > min_region_num:
            if subset.shape[0] > 100:
                overlap_file = "{}/{}.dynamic-state-{}.txt.gz".format(
                    dynamic_state_dir, soft_prefix, dynamic_state)
                pd.DataFrame(subset["region"]).to_csv(
                    overlap_file, sep='\t', header=False, index=False, compression="gzip")

        # TODO consider taking the dynamic marks (even if tiny sets) and super downweight their tasks?
        # this will cover for the issues above.
        # make sure that these groups are at least 100?
        
        
        
    return None


def order_and_viz_stable_epigenome(
        args,
        out_dir,
        prefix,
        histones,
        histone_files):
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
            args.annot["chromsizes"])

    # make the BED file
    args.integrative["atac_stable_histone_ordered_bed"] = "{}.bed.gz".format(
        args.integrative["atac_stable_w_histones_mat"].split(".mat")[0])
    if not os.path.isfile(args.integrative["atac_stable_histone_ordered_bed"]):
        make_bed = (
            "zcat {0} | "
            "grep -v region | "
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

    stable_plots_dir = "{}/plots".format(out_dir)
    os.system("mkdir -p {}".format(stable_plots_dir))
    
    # And plot
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
                args.integrative["atac_stable_histone_deeptools_bed"],
                histone_bigwigs,
                out_prefix,
                sort=False,
                kval=0,
                color=histone_color)

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

    # dynamic ATAC
    logging.info("EPIGENOME: DYNAMIC: integrate dynamic ATAC w histones")
    dynamic_epigenome_prefix = "ggr.epigenome.dynamic"

    activating_marks = histones[0:2]
    activating_mark_files = histone_files[0:2]

    args = order_and_viz_dynamic_epigenome_v2(
        args,
        args.folders["epigenome_dynamic_dir_v2"],
        dynamic_epigenome_prefix,
        activating_marks,
        activating_mark_files)


    print args.atac["final_soft_clusters"]
    print len(args.atac["final_soft_clusters"])
    cluster_by_chromatin_state(
        args.atac["final_soft_clusters"],
        args.integrative["atac_dynamic_w_histones_mat"],
        activating_marks,
        args.folders["epigenome_dynamic_clusters_dir"],
        dynamic_epigenome_prefix)
    
    quit()

    # TODO(dk) now for dynamic sets use soft clusters, split, and only keep "large enough" clusters
    # few levels: 1) just is_marked (5 possibilities for two histone marks)
    # then granular chromatin state. for granular level, remove things that are not at least some cutoff val (empirical)
    

    

    
    
    # stable ATAC
    logging.info("EPIGENOME: STABLE: integrate dynamic ATAC w histones")
    stable_epigenome_prefix = "ggr.epigenome.stable"
    args = order_and_viz_stable_epigenome(
        args,
        args.folders["epigenome_stable_dir"],
        stable_epigenome_prefix,
        histones,
        histone_files)

    quit()

    
    # get overlap with stable ATAC
    logging.info("EPIGENOME: STABLE: integrate dynamic ATAC w histones")
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
                sort=True,
                kval=0,
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
