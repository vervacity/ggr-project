"""Contains functions to run ATAC-specific analyses
"""

import os
import glob
import signal
import logging

import pandas as pd

from ggr.timeseries import get_consistent_dpgp_trajectories

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel

from ggr.util.bed_utils import merge_regions
from ggr.util.filtering import filter_for_ids
from ggr.util.counting import make_count_matrix
from ggr.util.counting import split_count_matrix_by_replicate

from ggr.util.bioinformatics import make_deeptools_heatmap


def run(args):
    """Run all ATAC analyses
    """
    atac_prefix = 'ggr.atac'
    
    # 1) generate a master regions file
    logging.info("ATAC: Generating master file...")
    args.atac['master_bed'] = '{0}/{1}.idr.master.bed.gz'.format(
        args.folders["data_dir"], atac_prefix)
    if not os.path.isfile(args.atac["master_bed"]):
        atac_peak_files = sorted(
            glob.glob('{0}/{1}'.format(
                args.atac['data_dir'],
                args.atac['idr_peak_glob'])))
        merge_regions(atac_peak_files, args.atac['master_bed'])
    
    # 2) get read counts in these regions (for DESeq2 analysis)
    logging.info("ATAC: Generating read counts per master region...")
    adjustment = "ends"
    args.atac['counts'] = '{0}/{1}.{2}.counts.mat.txt.gz'.format(
        args.folders['data_dir'], atac_prefix, adjustment)
    if not os.path.isfile(args.atac['counts']):
        # glob tagalign files
        atac_tagalign_files = sorted(
            glob.glob('{0}/{1}'.format(
                args.atac['data_dir'],
                args.atac['tagalign_glob'])))
        # remove files from media influenced timepoints
        for timepoint_string in args.misc["media_timepoints"]:
            atac_tagalign_files = [
                filename for filename in atac_tagalign_files
                if timepoint_string not in filename]
        # make count matrix
        make_count_matrix(
            args.atac['master_bed'],
            atac_tagalign_files,
            args.atac['counts'],
            "ATAC",
            adjustment="ends",
            tmp_dir=args.folders["atac_dir"])

    # 3) with read counts run DESeq2 as timeseries mode
    logging.info("ATAC: Running DESeq2 on all pairs of timepoints...")
    args.atac['dynamic_ids'] = '{0}/{1}.dynamic.ids.txt.gz'.format(
        args.folders['atac_timeseries_dir'], atac_prefix)
    if not os.path.isfile(args.atac['dynamic_ids']):
	run_timeseries_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3} full".format(
            args.atac['counts'],
            '{0}/{1}'.format(
                args.folders['atac_deseq2_dir'],
                atac_prefix),
            args.params['deseq2_fdr'],
            args.atac['dynamic_ids'])
        print run_timeseries_deseq2
        os.system(run_timeseries_deseq2)

    # 4) trajectories: make normalized rep1, rep2, pooled files (for DP_GP analysis)
    counts_prefix = args.atac["counts"].split(".mat")[0]
    args.atac["counts_rep1"] = "{}.rep1.mat.txt.gz".format(counts_prefix)
    args.atac["counts_rep2"] = "{}.rep2.mat.txt.gz".format(counts_prefix)
    args.atac["counts_pooled"] = "{}.pooled.mat.txt.gz".format(counts_prefix)
    
    logging.info("ATAC: Separating count mat into reps/pooled...")
    if not os.path.isfile(args.atac["counts_pooled"]):
        split_count_matrix_by_replicate(
            args.atac["counts"],
            args.atac["counts_rep1"],
            args.atac["counts_rep2"],
            args.atac["counts_pooled"])

    logging.info("ATAC: normalizing count files...")
    if not os.path.isfile("{}.rlog.mat.txt.gz".format(
            args.atac["counts_pooled"].split(".mat")[0])):
        run_rlogs = "normalize_count_mats.R {0} {1} {2} {3}".format(
            args.atac["counts"],
            args.atac["counts_rep1"],
            args.atac["counts_rep2"],
            args.atac["counts_pooled"])
        print run_rlogs
        os.system(run_rlogs)
    
    # for each of the counts files:
    logging.info("ATAC: filtering dynamic/stable...")
    counts_files_to_normalize = ["counts_rep1", "counts_rep2", "counts_pooled"]
    for counts_handle in counts_files_to_normalize:

        rlog_handle = "{}_rlog".format(counts_handle)
        args.atac[rlog_handle] = "{}.rlog.mat.txt.gz".format(
            args.atac[counts_handle].split('.mat')[0])

        # 6) filter appropriate count matrices for the dynamic and stable ones
        dynamic_handle = "{}_dynamic".format(rlog_handle)
        args.atac[dynamic_handle] = "{}.dynamic.mat.txt.gz".format(
            args.atac[rlog_handle].split(".mat")[0])
        if not os.path.isfile(args.atac[dynamic_handle]):
            filter_for_ids(args.atac[rlog_handle],
                           args.atac["dynamic_ids"],
                           args.atac[dynamic_handle])

        # only filter stable for the POOLED set
        if "pooled" in counts_handle:
            stable_handle = "{}_stable".format(rlog_handle)
            args.atac[stable_handle] = "{}.stable.mat.txt.gz".format(
                args.atac[rlog_handle].split('.mat')[0])
            if not os.path.isfile(args.atac[stable_handle]):
                filter_for_ids(args.atac[rlog_handle],
                               args.atac["dynamic_ids"],
                               args.atac[stable_handle],
                               opposite=True)
                
            # also make a BED file for downstream analyses
            stable_bed_handle = "{}_bed".format(stable_handle)
            args.atac[stable_bed_handle] = "{}.bed.gz".format(
                args.atac[stable_handle].split(".mat")[0])
            if not os.path.isfile(args.atac[stable_bed_handle]):
                make_bed = (
                    "zcat {} | "
                    "awk -F '\t' '{{ print $1 }}' | "
                    "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
                    "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
                    "grep -v d0 | "
                    "awk 'NF' | "
                    "gzip -c > {}").format(
                        args.atac[stable_handle],
                        args.atac[stable_bed_handle])
                print make_bed
                os.system(make_bed)
            

    # break for now
    if True:
        return args
                
    # 7) Run trajectories (on dynamic set) using DP_GP clustering
    args.atac["consistent_clusters"] = get_consistent_dpgp_trajectories(
        args.atac["counts_rep1_rlog_dynamic"],
        args.atac["counts_rep2_rlog_dynamic"],
        args.atac["counts_pooled_rlog_dynamic"],
        args.folders["atac_dp-gp_dir"],
        atac_prefix)

    quit()
    
    # 8) Take trajectories, run hclust to renumber clusters in order of hclust
    reorder_clusters = ("reorder_clusters_w_hclust.R {0} {1} {2}").format()

    

    # and make into a BED file with renumbered clusters

    
    # make BED (to overlap with histones)
    args.atac["clusters_bed_orig"] = "{}.bed.gz".format(
        args.atac["consistent_clusters"].split('.txt')[0])
    if not os.path.isfile(args.atac["clusters_bed_orig"]):
        make_bed = (
            "cat {0} | "
            "awk -F '\t' '{{ print $1\"\t\"$4 }}' | "
            "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
            "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
            "grep -v regions | "
            "sort -k1,1 -k2,2n | "
            "gzip -c > {1}").format(
                args.atac["consistent_clusters"],
                args.atac["clusters_bed"])
        print make_bed
        os.system(make_bed)



    quit()

    # MOVE THIS TO INTEGRATIVE
    # 8) split the clusters into bed files and run GREAT/HOMER
    # TODO convert this into a _from_dir method to make easy
    if False:
        data = pd.read_table(args.atac["dp-gp_clusters"])
        print data.columns
        for i in range(len(list(set(data["cluster"].tolist())))):
            cluster_idx = i + 1
            out_id_file = "{}/{}.dp-gp.cluster_{}.txt".format(
                args.folders["atac_dp-gp_dir"], atac_prefix, cluster_idx)
            print i
            data_subset = data[data["cluster"] == cluster_idx]
            data_subset["gene"].to_csv(out_id_file, sep='\t', header=False, index=False)

            # make bed
            out_bed_file = "{}/{}.dp-gp.cluster_{}.bed".format(
                args.folders["atac_dp-gp_dir"], atac_prefix, cluster_idx)
            make_bed = ("cat {} | "
                        "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
                        "awk -F '-' '{{ print $1\"\t\"$2 }}' > "
                        "{}").format(out_id_file, out_bed_file)
            print make_bed
            os.system(make_bed)

            # run great
            run_great = "run_rgreat.R {} {}/{}.dp-gp_cluster_{}".format(
                out_bed_file,
                args.folders["atac_dp-gp_dir"],
                atac_prefix,
                cluster_idx)
            print run_great
            os.system(run_great)
        
        
    # Make plotting functions to plot trajectories nicely
    # TODO generate heatmap as well as inset trajectory patterns
    if True:
        plot_trajectories = ("plot_trajectories.R {0} {1} {2}").format(
            args.atac["counts_rlog_nomedia_mat"],
            args.atac["dp-gp_clusters"],
            "{0}/{1}.dp_gp".format(args.folders["atac_dp-gp_dir"], atac_prefix))
        print plot_trajectories
        os.system(plot_trajectories)

    quit()

    # save key file
    args.atac["dp-gp_subsample"] = "{}/{}.dp_gp.subsampled.ordered.txt".format(
        args.folders["atac_dp-gp_dir"], atac_prefix)


    
    # Run homer on each group


    # as well as GREAT

    

    # overlap with histone information

    # TODO have a BED file for each histone mark with clusters


    # for each histone mark, do an intersect with the ATAC peaks.
    # to reconcile multiple hits for a peak, take nearest peak to midpoint of ATAC peak.


    
    # TODO - actually use peak calls to determine overlaps and marks
    args.atac["dp-gp_subsample_bed"] = "{0}/{1}.dp_gp.subsampled.ordered.bed".format(
        args.folders["atac_dp-gp_dir"], atac_prefix)
    get_atac_midpoints = (
        "cat {0} | "
        "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '\t' 'BEGIN {{ OFS=\"\t\" }} {{ $2=$2+int(($3-$2)/2); $3=$2+1; print }}' "
        "> {1}").format(args.atac["dp-gp_subsample"],
                        args.atac["dp-gp_subsample_bed"])

    print get_atac_midpoints
    os.system(get_atac_midpoints)

    # and plot with the midpoint file (testing HJ3K27ac)
    histones = ["H3K27ac", "H3K4me1", "H3K27me3"]

    for histone in histones:
    
        histone_bigwigs = sorted(
            glob.glob("{}/{}".format(
                args.chipseq["data_dir"],
                args.chipseq["histones"][histone]["pooled_bigwig_glob"])))
        make_deeptools_heatmap(
            args.atac["dp-gp_subsample_bed"],
            histone_bigwigs,
            "{0}/{1}.{2}_profile".format(args.folders["atac_dp-gp_dir"],
                                         atac_prefix,
                                         histone),
        sort=False)
    
    
    quit()
    

    


    return None
