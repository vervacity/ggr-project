"""Contains functions to run ATAC-specific analyses
"""

import os
import glob
import signal

import pandas as pd

from ggr.util.bed_utils import merge_regions
from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel
from ggr.util.filtering import remove_media_timepoints, filter_for_ids
from ggr.util.counting import make_count_matrix

from ggr.util.bioinformatics import make_deeptools_heatmap


def run(args):
    """Run everything
    """
    atac_prefix = 'ggr.atac'
    
    # generate a master regions file
    args.atac['master_bed'] = '{0}/{1}.idr.master.bed.gz'.format(
        args.folders["data_dir"], atac_prefix)
    atac_peak_files = glob.glob('{0}/{1}'.format(
        args.atac['data_dir'], args.atac['idr_peak_glob']))
    merge_regions(atac_peak_files, args.atac['master_bed'])

    # TODO any plots here?
    # chromatin state annotations? enh/prom/etc?
    # ChromHMM enrichments
    
    
    # get read counts in these regions
    adjustment = "ends"
    args.atac['counts_mat'] = '{0}/{1}.cuts.counts.mat.txt.gz'.format(
        args.folders['data_dir'], atac_prefix)
    atac_tagalign_files = glob.glob('{0}/{1}'.format(args.atac['data_dir'],
                                                     args.atac['tagalign_glob']))
    if not os.path.isfile(args.atac['counts_mat']):
        make_count_matrix(args.atac['master_bed'],
                               atac_tagalign_files,
                               args.atac['counts_mat'],
                               adjustment="ends",
                               tmp_dir=args.folders["atac_dir"])

    # use DESeq to run rlog normalization (for visualization and timeseries)
    args.atac['counts_rlog_mat'] = '{0}/{1}.cuts.rlog.mat.txt.gz'.format(args.folders["data_dir"], atac_prefix)
    if not os.path.isfile(args.atac['counts_rlog_mat']):
        run_rlog = 'normalize_count_mat.R {0} {1}'.format(
            args.atac['counts_mat'], 
            args.atac['counts_rlog_mat'])
        print run_rlog
        os.system(run_rlog)
        
    # Remove media influenced timepoints first (both counts and rlog norm matrices)
    args.atac['counts_nomedia_mat'] = '{}.nomedia.mat.txt.gz'.format(
        args.atac['counts_mat'].split('.mat')[0])
    if not os.path.isfile(args.atac['counts_nomedia_mat']):
	remove_media_timepoints(
            args.atac['counts_mat'],
            args,
            args.atac['counts_nomedia_mat'])
    args.atac['counts_rlog_nomedia_mat'] = '{}.nomedia.mat.txt.gz'.format(
        args.atac['counts_rlog_mat'].split('.mat')[0])
    if not os.path.isfile(args.atac['counts_rlog_nomedia_mat']):
	remove_media_timepoints(
            args.atac['counts_rlog_mat'],
            args,
            args.atac['counts_rlog_nomedia_mat'])
    
    # with read counts run DESeq2 as timeseries mode
    args.atac['dynamic_region_ids'] = '{0}/{1}.dynamic.ids.txt.gz'.format(
        args.folders['atac_timeseries_dir'],
        atac_prefix)
    if not os.path.isfile(args.atac['dynamic_region_ids']):
	run_timeseries_deseq2 = "run_deseq2_timediff.R {0} {1} {2} {3}".format(
            args.atac['counts_nomedia_mat'],
            '{0}/deseq2/{1}'.format(args.folders['atac_timeseries_dir'], atac_prefix),
            args.params['deseq2_fdr'],
            args.atac['dynamic_region_ids'])
        print run_timeseries_deseq2
        os.system(run_timeseries_deseq2)
        
    # filter matrix for the dynamic and stable ones
    args.atac['dynamic_mat'] = '{0}/{1}.dynamic.mat.txt.gz'.format(args.folders['atac_dynamic_dir'], atac_prefix)
    os.system('mkdir -p {}'.format(args.folders['atac_dynamic_dir']))
    if not os.path.isfile(args.atac['dynamic_mat']):
        filter_for_ids(args.atac['counts_rlog_nomedia_mat'],
                       args.atac['dynamic_region_ids'],
                       args.atac['dynamic_mat'])

    args.atac['stable_mat'] = '{0}/{1}.stable.mat.txt.gz'.format(args.folders['atac_stable_dir'], atac_prefix)
    os.system('mkdir -p {}'.format(args.folders['atac_stable_dir']))
    if not os.path.isfile(args.atac['stable_mat']):
        filter_for_ids(args.atac['counts_rlog_nomedia_mat'],
                       args.atac['dynamic_region_ids'],
                       args.atac['stable_mat'],
                       opposite=True)

    # quick plot of numbers (dynamic vs stable)?
    
    
    # Run trajectories (on dynamic set) using DP_GP clustering
    if not os.path.isdir(args.folders['atac_dp-gp_dir']):
        os.system('mkdir -p {}'.format(args.folders['atac_dp-gp_dir']))
        # unzip into tmp file
        tmp_unzipped_mat = '{0}/{1}.dp_gp.tmp.txt'.format(args.folders['atac_dp-gp_dir'], atac_prefix)
        unzip_mat = ("zcat {0} > {1}").format(args.atac['dynamic_mat'], tmp_unzipped_mat)
        print unzip_mat
        os.system(unzip_mat)
        # cluster
        cluster_dp_gp = (
            "DP_GP_cluster.py "
            "-i {0} "
            "-o {1}"
            "-p png --plot").format(
                tmp_unzipped_mat,
                '{0}/{1}'.format(args.folders['atac_dp-gp_dir'], atac_prefix))
        print cluster_dp_gp
        os.system(cluster_dp_gp)
        # delete tmp file
        os.system('rm {}'.format(tmp_unzipped_mat))

    # save key file to args - the clustering
    args.atac["dp-gp_clusters"] = "{0}/{1}_optimal_clustering.txt".format(
        args.folders["atac_dp-gp_dir"], atac_prefix)

    # split the clusters into bed files
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
