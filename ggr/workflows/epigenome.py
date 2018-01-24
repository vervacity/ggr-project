# workflows integrating epigenomic datasets

import os
import gzip
import glob
import logging

from ggr.util.utils import run_shell_cmd

from ggr.analyses.bioinformatics import make_deeptools_heatmap

from ggr.util.bed_utils import id_to_bed

from ggr.analyses.epigenome import get_histone_overlaps
from ggr.analyses.epigenome import cluster_by_chromatin_marks
from ggr.analyses.epigenome import split_stable_atac_by_dynamic_marks

from ggr.analyses.filtering import sort_by_clusters
from ggr.analyses.filtering import get_ordered_subsample

from ggr.analyses.timeseries import plot_clusters


def make_fake_cluster_file(id_list_file, out_file):
    """Make a fake cluster file
    """
    with gzip.open(id_list_file, "r") as fp:
        with gzip.open(out_file, "w") as out:
            out.write("cluster\tid\n")
            for line in fp:
                if "d00" in line:
                    continue
                out.write("0\t{}\n".format(line.strip().split("\t")[0]))
                
    return None


def run_dynamic_epigenome_workflow(
        args,
        prefix,
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
    logger.info("ANALYSIS: overlap ATAC with histones")
    atac_dynamic_bed_key = "atac.dynamic.bed"
    histone_overlap_dir = "{}/overlap_histone".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(histone_overlap_dir))

    histone_overlap_mat_key = "{}.overlap_histones.mat".format(atac_dynamic_bed_key.split(".bed")[0])
    out_results[histone_overlap_mat_key] = "{}/{}.overlaps.mat.txt.gz".format(
        histone_overlap_dir, prefix)

    # set up activating histone file sets
    histones = args.inputs["chipseq"][args.cluster]["histones"]["ordered_names"]
    histone_files = [
        (histone,
         args.outputs["results"]["histones"][histone]["timeseries"]["ggr.histone.{}.enumerated.bed".format(histone)],
         args.inputs["params"]["histones"][histone]["overlap_extend_len"])
        for histone in histones]
    activating_marks = histones[0:2]
    activating_mark_files = histone_files[0:2]
    
    if not os.path.isfile(out_results[histone_overlap_mat_key]):
        # NOTE: this currently requires bedtools 2.23.0
        # fix this for upward compatibility
        get_histone_overlaps(
            out_data[atac_dynamic_bed_key],
            activating_mark_files,
            histone_overlap_dir,
            out_results[histone_overlap_mat_key],
            args.inputs["annot"][args.cluster]["chromsizes"],
            histone_assignment=histone_assignment)

    # TODO factor this out to subworkflow?
    # -----------------------------------------
    # ANALYSIS 1 - reorder ATAC with all clusters
    # inputs: dynamic ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: reorder ATAC regions using ATAC/histone clusters")
    atac_clusters_key = "clusters.reproducible.hard.reordered.list"
    atac_cluster_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][atac_clusters_key]
    atac_cluster_col = "cluster"

    histone_clusters_key = histone_overlap_mat_key
    histone_cluster_file = out_results[histone_clusters_key]
    histone_cluster_col = "histone_cluster"

    cluster_files = [
        (atac_cluster_file, [atac_cluster_col]),
        (histone_cluster_file, [histone_cluster_col])]

    # reorder full list of ids
    atac_ordered_key = "atac.epigenome_ordered.list"
    out_results[atac_ordered_key] = "{}/{}.ordered.txt.gz".format(
        results_dir, prefix)
    atac_ordered_clusters_key = "atac.epigenome_ordered.mat"
    out_results[atac_ordered_clusters_key] = "{}/{}.ordered.clusters.txt".format(
        results_dir, prefix)

    if not os.path.isfile(out_results[atac_ordered_key]):
        sort_by_clusters(
            cluster_files,
            out_results[atac_ordered_clusters_key],
            out_results[atac_ordered_key])

    # and generate a subsample for plotting too
    atac_ordered_subsample_key = "{}.subsample.list".format(
        atac_ordered_key.split(".list")[0])
    out_results[atac_ordered_subsample_key] = "{}.subsampled.txt".format(
        out_results[atac_ordered_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_key]):
        get_ordered_subsample(
            out_results[atac_ordered_clusters_key],
            out_results[atac_ordered_subsample_key])

    # produce a subsample BED too
    atac_ordered_subsample_bed_key = "{}.bed".format(
        atac_ordered_subsample_key.split(".list")[0])
    out_results[atac_ordered_subsample_bed_key] = "{}.bed".format(
        out_results[atac_ordered_subsample_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_bed_key]):
        id_to_bed(
            out_results[atac_ordered_subsample_key],
            "{}.gz".format(out_results[atac_ordered_subsample_bed_key]))

        # and unzip
        run_shell_cmd("gunzip {}.gz".format(out_results[atac_ordered_subsample_bed_key]))

    # -----------------------------------------
    # ANALYSIS 2 - separate out cluster files (ids and beds)
    # inputs: cluster file
    # outputs: separated id and bed files
    # -----------------------------------------
    logger.info("ANALYSIS: separate out cluster files")
    mark_dir = "{}/clusters/by_mark".format(results_dir)
    state_dir = "{}/clusters/by_state".format(results_dir)
    if not os.path.isdir(state_dir):
        run_shell_cmd("mkdir -p {0}/ids {0}/bed {1}/ids {1}/bed".format(
            mark_dir, state_dir))
        cluster_by_chromatin_marks(
            out_results[atac_ordered_clusters_key],
            ["H3K27ac", "H3K4me1"],
            "histone_cluster",
            "{}/ids".format(mark_dir),
            "{}/ids".format(state_dir),
            prefix,
            min_region_num=100)

        # and now convert files to beds and keep bed dir
        mark_files = glob.glob("{}/ids/*.txt.gz".format(mark_dir))
        for mark_file in mark_files:
            mark_bed_file = "{0}/bed/{1}.bed.gz".format(
                mark_dir, os.path.basename(mark_file).split('.txt')[0])
            id_to_bed(mark_file, mark_bed_file, sort=True)
        
        state_files = glob.glob("{}/ids/*.txt.gz".format(state_dir))
        for state_file in state_files:
            state_bed_file = "{0}/bed/{1}.bed.gz".format(
                state_dir, os.path.basename(state_file).split('.txt')[0])
            id_to_bed(state_file, state_bed_file, sort=True)

    args.outputs["results"]["label_dirs"].append(mark_dir)
    args.outputs["results"]["label_dirs"].append(state_dir)
        
    # -----------------------------------------
    # ANALYSIS 3 - plot this out
    # inputs: dynamic ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: plot out the heatmaps with the ordering in the subsample")
    mat_key = "atac.counts.pooled.rlog.dynamic.mat"
    plot_dir = "{}/plots".format(results_dir)
    
    # plot ATAC
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))
        plot_clusters(
            atac_cluster_file,
            out_results[atac_ordered_subsample_key],
            out_data[mat_key],
            plot_dir,
            prefix)
        
        # TODO need to plot color bars too?
        
        # plot the histone signal profiles with deeptools
        histone_colors = args.inputs["chipseq"][args.cluster]["histones"]["ordered_deeptools_colors"]
        for histone_idx in range(len(histones)):
            histone = histones[histone_idx]
            histone_color = histone_colors[histone_idx]
            histone_bigwigs = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_bigwig_glob"])))

            out_prefix = "{}/{}.{}_overlap".format(
                plot_dir,
                prefix,
                histone)
            out_file = "{}.heatmap.profile.pdf".format(out_prefix)
        
            if not os.path.isfile(out_file):
                make_deeptools_heatmap(
                    out_results[atac_ordered_subsample_bed_key],
                    histone_bigwigs,
                    out_prefix,
                    sort=False,
                    referencepoint="center",
                    color=histone_color)

    # don't forget to store outputs back
    args.outputs["data"] = out_data    
    args.outputs["results"]["epigenome"][results_dirname] = out_results
            
    return args


def run_stable_epigenome_workflow(
        args,
        prefix,
        histone_assignment="nearest"):
    """Same thing as above, except all ATAC is one big cluster (no splits)
    """
    # logging and folder set up
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: get stable epigenome")

    # assertions
    
    
    # setup data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "stable"
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
    logger.info("ANALYSIS: overlap ATAC with histones")
    atac_stable_bed_key = "atac.stable.bed"
    histone_overlap_dir = "{}/overlap_histone".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(histone_overlap_dir))

    histone_overlap_mat_key = "{}.overlap_histones.mat".format(atac_stable_bed_key.split(".bed")[0])
    out_results[histone_overlap_mat_key] = "{}/{}.overlaps.mat.txt.gz".format(
        histone_overlap_dir, prefix)

    # set up activating histone file sets
    histones = args.inputs["chipseq"][args.cluster]["histones"]["ordered_names"]
    histone_files = [
        (histone,
         args.outputs["results"]["histones"][histone]["timeseries"]["ggr.histone.{}.enumerated.bed".format(histone)],
         args.inputs["params"]["histones"][histone]["overlap_extend_len"])
        for histone in histones]
    activating_marks = histones[0:2]
    activating_mark_files = histone_files[0:2]
    
    if not os.path.isfile(out_results[histone_overlap_mat_key]):
        # NOTE: this currently requires bedtools 2.23.0
        # fix this for upward compatibility
        get_histone_overlaps(
            out_data[atac_stable_bed_key],
            histone_files,
            histone_overlap_dir,
            out_results[histone_overlap_mat_key],
            args.inputs["annot"][args.cluster]["chromsizes"],
            histone_assignment=histone_assignment)

    # -----------------------------------------
    # ANALYSIS 1 - reorder ATAC with all clusters
    # inputs: stable ATAC ids, and all clusters
    # outputs: ordered ATAC file and lists of chromatin states
    # -----------------------------------------
    logger.info("ANALYSIS: reorder ATAC regions using histone clusters")

    # make a fake ATAC clusters file
    atac_stable_mat_key = "atac.counts.pooled.rlog.stable.mat"
    atac_fake_clusters_key = "atac.stable.fake_cluster.list"
    out_results[atac_fake_clusters_key] = "{}/{}.atac.clusters.txt.gz".format(results_dir, prefix)
    if not os.path.isfile(out_results[atac_fake_clusters_key]):
        make_fake_cluster_file(
            out_data[atac_stable_mat_key],
            out_results[atac_fake_clusters_key])
    atac_cluster_cols = ["cluster"]
        
    histone_clusters_key = histone_overlap_mat_key
    histone_cluster_file = out_results[histone_clusters_key]
    histone_cluster_cols = ["histone_cluster"]

    cluster_files = [
        (out_results[atac_fake_clusters_key], atac_cluster_cols),
        (histone_cluster_file, histone_cluster_cols)]

    # reorder full list of ids
    atac_ordered_key = "atac.epigenome_ordered.list"
    out_results[atac_ordered_key] = "{}/{}.ordered.txt.gz".format(
        results_dir, prefix)
    atac_ordered_clusters_key = "atac.epigenome_ordered.mat"
    out_results[atac_ordered_clusters_key] = "{}/{}.ordered.clusters.txt".format(
        results_dir, prefix)

    if not os.path.isfile(out_results[atac_ordered_key]):
        sort_by_clusters(
            cluster_files,
            out_results[atac_ordered_clusters_key],
            out_results[atac_ordered_key])

    # and generate a subsample for plotting too
    atac_ordered_subsample_key = "{}.subsample.list".format(
        atac_ordered_key.split(".list")[0])
    out_results[atac_ordered_subsample_key] = "{}.subsampled.txt".format(
        out_results[atac_ordered_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_key]):
        get_ordered_subsample(
            out_results[atac_ordered_clusters_key],
            out_results[atac_ordered_subsample_key])

    # produce a BED too
    atac_ordered_subsample_bed_key = "{}.bed".format(
        atac_ordered_subsample_key.split(".list")[0])
    out_results[atac_ordered_subsample_bed_key] = "{}.bed".format(
        out_results[atac_ordered_subsample_key].split(".txt")[0])
    if not os.path.isfile(out_results[atac_ordered_subsample_bed_key]):
        id_to_bed(
            out_results[atac_ordered_subsample_key],
            "{}.gz".format(out_results[atac_ordered_subsample_bed_key]))

        # and unzip
        run_shell_cmd("gunzip {}.gz".format(out_results[atac_ordered_subsample_bed_key]))
        
    # -----------------------------------------
    # ANALYSIS 2 - separate out cluster files (ids and beds)
    # inputs: cluster file
    # outputs: separated id and bed files
    # -----------------------------------------
    logger.info("ANALYSIS: separate out cluster files")
    mark_dir = "{}/clusters/by_mark".format(results_dir)
    state_dir = "{}/clusters/by_state".format(results_dir)
    if not os.path.isdir(state_dir):
        run_shell_cmd("mkdir -p {0}/ids {0}/bed {1}/ids {1}/bed".format(
            mark_dir, state_dir))
        cluster_by_chromatin_marks(
            out_results[atac_ordered_clusters_key],
            ["H3K27ac", "H3K4me1", "H3K27me3"],
            "histone_cluster",
            "{}/ids".format(mark_dir),
            "{}/ids".format(state_dir),
            prefix,
            min_region_num=100)

        # and now convert files to beds
        mark_files = glob.glob("{}/ids/*.txt.gz".format(mark_dir))
        for mark_file in mark_files:
            mark_bed_file = "{0}/bed/{1}.bed.gz".format(
                mark_dir, os.path.basename(mark_file).split('.txt')[0])
            id_to_bed(mark_file, mark_bed_file, sort=True)
        
        state_files = glob.glob("{}/ids/*.txt.gz".format(state_dir))
        for state_file in state_files:
            state_bed_file = "{0}/bed/{1}.bed.gz".format(
                state_dir, os.path.basename(state_file).split('.txt')[0])
            id_to_bed(state_file, state_bed_file, sort=True)

    args.outputs["results"]["label_dirs"].append(mark_dir)
    args.outputs["results"]["label_dirs"].append(state_dir)

    # -----------------------------------------
    # ANALYSIS 3 - separate out dynamic and non-dynamic chromatin mark data,
    # mostly for plotting
    # inputs: clusters
    # outputs: separate id lists
    # -----------------------------------------
    cluster_file_key = atac_ordered_clusters_key
    stable_atac_dynamic_histones_key = "{}.dynamic_histones.clusters".format(
        cluster_file_key.split(".mat")[0])
    out_results[stable_atac_dynamic_histones_key] = "{}/{}.dynamic_histones.clusters.txt".format(
        results_dir, os.path.basename(out_results[cluster_file_key]).split(".clusters")[0])
    stable_atac_stable_histones_key = "{}.stable_histones.clusters".format(
        cluster_file_key.split(".mat")[0])
    out_results[stable_atac_stable_histones_key] = "{}/{}.stable_histones.clusters.txt".format(
        results_dir, os.path.basename(out_results[cluster_file_key]).split(".clusters")[0])
    if not os.path.isfile(out_results[stable_atac_stable_histones_key]):
        # TODO: check that it does not reorder
        split_stable_atac_by_dynamic_marks(
            out_results[cluster_file_key],
            out_results[stable_atac_dynamic_histones_key],
            out_results[stable_atac_stable_histones_key])
        
    # keep keys for future use
    stable_set_keys = [
        stable_atac_dynamic_histones_key,
        stable_atac_stable_histones_key]

    # get subsamples, and subsample bed files for each of these
    subsample_bed_keys = []
    for stable_set_key in stable_set_keys:
        # and generate a subsample for plotting too
        subsample_mat_key = "{}.subsample.mat".format(
            stable_set_key.split(".list")[0])
        out_results[subsample_mat_key] = "{}.subsampled.mat".format(
            out_results[stable_set_key].split(".txt")[0])
        if not os.path.isfile(out_results[subsample_mat_key]):
            get_ordered_subsample(
                out_results[stable_set_key],
                out_results[subsample_mat_key])

        # extract just the id
        subsample_key = "{}.subsample.list".format(
            stable_set_key.split(".list")[0])
        out_results[subsample_key] = "{}.subsampled.txt".format(
            out_results[stable_set_key].split(".txt")[0])
        if not os.path.isfile(out_results[subsample_key]):
            get_ids = "cat {} | awk '{{print $2}}' | grep -v id > {}".format(
                out_results[subsample_mat_key],
                out_results[subsample_key])
            run_shell_cmd(get_ids)
        
        # produce a BED too
        subsample_bed_key = "{}.bed".format(
            subsample_key.split(".list")[0])
        out_results[subsample_bed_key] = "{}.bed".format(
            out_results[subsample_key].split(".txt")[0])
        subsample_bed_keys.append(subsample_bed_key)
        if not os.path.isfile(out_results[subsample_bed_key]):
            id_to_bed(
                out_results[subsample_key],
                "{}.gz".format(out_results[subsample_bed_key]))

            # and unzip
            run_shell_cmd("gunzip {}.gz".format(out_results[subsample_bed_key]))
        
    # -----------------------------------------
    # ANALYSIS 4 - plot this out
    # inputs: stable groups
    # outputs: heatmaps with histone marks
    # -----------------------------------------
    logger.info("ANALYSIS: plot out the heatmaps with the ordering in the subsamples")
    plot_dir = "{}/plots".format(results_dir)
    
    # plot ATAC
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))

        for subsample_bed_key in subsample_bed_keys:
        
            # TODO need to plot color bars?
            
            # plot the histone signal profiles with deeptools
            histone_colors = args.inputs["chipseq"][args.cluster]["histones"]["ordered_deeptools_colors"]
            for histone_idx in range(len(histones)):
                histone = histones[histone_idx]
                histone_color = histone_colors[histone_idx]
                histone_bigwigs = sorted(
                    glob.glob("{}/{}".format(
                        args.inputs["chipseq"][args.cluster]["data_dir"],
                        args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_bigwig_glob"])))

                out_prefix = "{}/{}.{}_overlap".format(
                    plot_dir,
                    prefix,
                    histone)
                out_file = "{}.heatmap.profile.pdf".format(out_prefix)
                
                if not os.path.isfile(out_file):
                    make_deeptools_heatmap(
                        out_results[subsample_bed_key],
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
    args = run_dynamic_epigenome_workflow(
        args, "{}.dynamic".format(prefix))

    # -----------------------------------------
    # ANALYSIS 1 - integrate histones into stable ATAC
    # inputs: dynamic ATAC clusters, histone marks
    # outputs: chromatin states, plots
    # -----------------------------------------
    logger.info("ANALYSIS: integrate stable ATAC w histones")
    args = run_stable_epigenome_workflow(
        args, "{}.stable".format(prefix))
    
    
    # -----------------------------------------
    # ANALYSIS 2 - run bioinformatics tooks
    # inputs: BED files
    # outputs: TF enrichments (HOMER), GO enrichments (GREAT)
    # -----------------------------------------
    logger.info("ANALYSIS: running bioinformatics on chromatin mark/state files")
    
    


    

    return args
