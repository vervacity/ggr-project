"""Contains code to run histone analyses
"""

import os
import gzip
import glob
import logging
import signal

import pandas as pd

from ggr.analyses.counting import make_count_matrix
from ggr.analyses.counting import pull_single_replicate
from ggr.analyses.utils import build_id_matching_mat
from ggr.analyses.utils import plot_PCA

from ggr.workflows.timeseries import run_timeseries_enumeration_workflow

from ggr.util.bed_utils import merge_regions
from ggr.util.bed_utils import get_midpoint_and_extend
from ggr.util.diff import join_diff_region_lists_to_mat
from ggr.util.utils import run_shell_cmd



def runall(args, prefix):
    """workflow for histone analyses
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("MASTER_WORKFLOW: run histone analyses")

    # set up inputs
    inputs = args.inputs["chipseq"][args.cluster]
    
    # set up histones
    histones = inputs["histones"]["ordered_names"]

    # set up data and results
    data_dir = args.outputs["data"]["dir"]
    out_data = args.outputs["data"]

    results_dirname = "histones"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]
    
    for histone in histones:

        logger.info("HISTONE: Working on {}".format(histone))
        histone_prefix = "{}.{}".format(prefix, histone)

        # set up histone specific dir
        histone_dir = "{}/{}".format(results_dir, histone)
        args.outputs["results"][results_dirname][histone] = {"dir": histone_dir}
        run_shell_cmd("mkdir -p {}".format(histone_dir))
        out_histone_results = args.outputs["results"][results_dirname][histone]

        # -------------------------------------------
        # ANALYSIS 0 - set up master regions
        # input: peak files
        # output: master regions
        # -------------------------------------------
        logger.info("ANALYSIS:{}: generate master regions file".format(histone))
        master_regions_key = "{}.master.bed".format(histone)
        out_data[master_regions_key] = "{0}/{1}.overlap.master.bed.gz".format(
            data_dir, histone_prefix)
        histone_overlap_files = sorted(
            glob.glob("{0}/{1}".format(
                inputs["data_dir"],
                inputs["histones"][histone]["overlap_glob"])))
        if not os.path.isfile(out_data[master_regions_key]):
            # always generate a merged region set
            merge_regions(histone_overlap_files, out_data[master_regions_key])

        # build a master file based on ATAC midpoints
        atac_master_key = "atac.master.bed"
        extend_len = args.inputs["params"]["histones"][histone]["overlap_extend_len"]
        atac_midpoints_key = "atac.midpoints.slop_{}.bed".format(histone, extend_len)
        out_data[atac_midpoints_key] = "{}.slop_{}bp.bed.gz".format(
            out_data[atac_master_key].split(".bed")[0],
            extend_len)
        if not os.path.isfile(out_data[atac_midpoints_key]):
            get_midpoint_and_extend(
                out_data[atac_master_key],
                args.inputs["annot"][args.cluster]["chromsizes"],
                extend_len,
                out_data[atac_midpoints_key])

        # now intersect ATAC with the naive overlap files and only keep region if has an overlap
        histone_atac_regions_key = "{}.overlap_{}.bed".format(
            atac_midpoints_key.split(".bed")[0], histone)
        out_data[histone_atac_regions_key] = "{}.overlap_{}.bed.gz".format(
            out_data[atac_midpoints_key].split(".bed")[0], histone)
        if not os.path.isfile(out_data[histone_atac_regions_key]):
            keep_marked = (
                "bedtools intersect -u -a {0} -b {1} | "
                "gzip -c > {2}").format(
                    out_data[atac_midpoints_key],
                    out_data[master_regions_key],
                    out_data[histone_atac_regions_key])
            run_shell_cmd(keep_marked)

        # choose master regions
        master_regions_key = histone_atac_regions_key

        # -------------------------------------------
        # ANALYSIS 1 - get read counts per master region
        # input: peak files
        # output: master regions
        # -------------------------------------------
        logger.info("ANALYSIS:{}: get read counts per master region".format(histone))
        adjustment = "midpoints"
        counts_key = "{}.counts.mat".format(histone)
        out_data[counts_key] = "{0}/{1}.midpoints.counts.mat.txt.gz".format(
            data_dir, histone_prefix)
        if not os.path.isfile(out_data[counts_key]):
            bedpe_files = sorted(
                glob.glob("{}/{}".format(
                    inputs["data_dir"],
                    inputs["histones"][histone]["bedpe_glob"])))
            make_count_matrix(
                out_data[master_regions_key],
                bedpe_files,
                out_data[counts_key],
                histone,
                adjustment=adjustment,
                tmp_dir=histone_dir)

        args.outputs["data"] = out_data
        args.outputs["results"][results_dirname][histone] = out_histone_results
        
        # -------------------------------------------
        # ANALYSIS 3 - run DESeq on forward sequential pairs of timepoints
        # input: count matrix
        # output: clusters of regions (enumeration
        # -------------------------------------------
        logger.info("ANALYSIS:{}: Running timeseries enumeration workflow...".format(histone))
        args = run_timeseries_enumeration_workflow(
            args,
            histone_prefix,
            datatype_key="histones",
            subtype_key=histone,
            master_regions_key=master_regions_key,
            mat_key=counts_key)

        # -------------------------------------------
        # ANALYSIS - make PCA plots
        # input: count matrix
        # output: pca plot
        # -------------------------------------------
        logger.info("ANALYSIS:{}: PCA plots...".format(histone))

        # pull the replicates apart
        mat_key = "{}.counts.mat".format(histone)
        counts_prefix = out_data[mat_key].split(".mat")[0]
        rep_handles = []
        
        # rep1
        handle = "{}.rep1.mat".format(mat_key.split(".mat")[0])
        out_data[handle] = "{}.rep1.mat.txt.gz".format(counts_prefix)
        if not os.path.isfile(out_data[handle]):
            pull_single_replicate(out_data[mat_key], out_data[handle], rep="b1")
        rep_handles.append(handle)
            
        # rep2
        handle = "{}.rep2.mat".format(mat_key.split(".mat")[0])
        out_data[handle] = "{}.rep2.mat.txt.gz".format(counts_prefix)
        if not os.path.isfile(out_data[handle]):
            pull_single_replicate(out_data[mat_key], out_data[handle], rep="b2")
        rep_handles.append(handle)

        # normalize
        rep1_rlog_mat_file = "{}.rlog.mat.txt.gz".format(out_data[rep_handles[0]].split(".mat")[0])
        if not os.path.isfile(rep1_rlog_mat_file):
            run_rlogs = "normalize.rlog.master_and_transfer.R {0} {1} {2}".format(
                out_data[mat_key], *[out_data[handle] for handle in rep_handles])
            run_shell_cmd(run_rlogs)

        # get PCA
        plot_dir = "{}/plots".format(histone_dir)
        if not os.path.isdir(plot_dir):
            run_shell_cmd("mkdir -p {}".format(plot_dir))

            # pull the two reps
            histone_mat_files = sorted(glob.glob(
                "{}/ggr.histone.{}.*.rep*rlog.mat.txt.gz".format(
                    args.outputs["data"]["dir"], histone)))
            
            # and build filtered mat files based on final dynamic set
            final_dynamic_list = args.outputs["results"]["histones"][histone]["timeseries"][
                "ggr.histone.{}.dynamic.ids.list".format(histone)]
            
            filt_mat_files = []
            for histone_mat_file in histone_mat_files:
                filt_histone_mat_file = "{}/{}.filt.mat.txt.gz".format(
                    plot_dir,
                    os.path.basename(histone_mat_file).split(".mat")[0])
                build_id_matching_mat(
                    final_dynamic_list,
                    histone_mat_file,
                    filt_histone_mat_file,
                    keep_cols=[],
                    primary_id_col=None)
                filt_mat_files.append(filt_histone_mat_file)

            # plot PCA/correlation (bio reps separately and pooled)
            pca_file = "{}/{}.pca.pdf".format(
                plot_dir,
                os.path.basename(histone_mat_files[0]).split(".rep")[0])
            plot_PCA(filt_mat_files, pca_file)
        
    return args
