# workflows integrating epigenomic datasets

import os
import glob
import logging



def run_dynamic_epigenome_workflow(
        args,
        prefix,
        #out_dir,
        #prefix,
        #activating_histones,
        #activating_histone_files,
        #all_histones,
        histone_assignment="nearest"):
    """Order clusters by hclust and histones more simple
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get dynamic epigenome")

    # assertions
    
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "dynamic"
    results_dir = "{}/{}".format(
        args.outputs["results"]["epigenome"]["dir"],
        results_dirname)
    args.outputs["results"]["epigenome"][results_dirname] = {
        "dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"]["epigenome"][results_dirname]

    # -----------------------------------------
    # ANALYSIS 0 - get overlap with histones
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    histone_overlap_dir = "{}/overlap_histone".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(histone_overlap_dir))

    histone_overlap_mat_key = "atac.dynamic.overlap_histones.mat"
    out_results[histone_overlap_mat_key] = "{}/{}.overlaps.mat.txt.gz".format(
        histone_overlap_dir, prefix)

    
    if not os.path.isfile(out_results[histone_overlap_mat_key]):
        get_histone_overlaps(
            out_data["counts.pooled.rlog.dynamic.bed"],

            
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
    run_shell_cmd("mkdir -p {}".format(dynamic_plots_dir))
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
        run_shell_cmd(order_and_plot_atac)
        
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
        run_shell_cmd(make_bed)

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





def runall(args, prefix):
    """Integrate epigenomic datasets
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("MASTER_WORKFLOW: run epigenomic analyses")

    # set up data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "epigenome"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]

    # -----------------------------------------
    # ANALYSIS 0 - integrate histones into dynamic ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    logger.info("ANALYSIS: integrate dynamic ATAC w histones")
    
    

    


    # -----------------------------------------
    # ANALYSIS 1 - integrate histones into stable ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------




    

    return args
