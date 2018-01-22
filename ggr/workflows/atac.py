"""ggr workflow for atac analyses
"""

import os
import glob
import signal
import logging

from ggr.util.utils import run_shell_cmd
from ggr.util.utils import parallel_copy

from ggr.util.bed_utils import merge_regions
from ggr.util.bed_utils import id_to_bed
from ggr.util.filtering import filter_for_ids

from ggr.analyses.counting import make_count_matrix

from ggr.workflows.timeseries import run_timeseries_workflow


def runall(args, prefix):
    """all workflows for atac-seq data
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: run atac analyses")

    # set up inputs
    inputs = args.inputs["atac"][args.cluster]
    
    # assertions
    assert inputs.get("data_dir") is not None
    assert inputs.get("idr_peak_glob") is not None

    # set up data
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]
    
    results_dirname = "atac"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]

    # -------------------------------------------
    # ANALYSIS 0 - download timepoint BED files
    # input: peak files
    # output: same peak files, new location
    # -------------------------------------------
    logger.info("ANALYSIS: copy idr peak files to new dir")
    atac_peak_files = sorted(
        glob.glob('{0}/{1}'.format(
            inputs['data_dir'],
            inputs['idr_peak_glob'])))
    timepoint_dir = "{}/peaks.timepoints".format(results_dir)
    out_results["timepoint_region_dir"] = timepoint_dir
    timepoints_files = glob.glob("{}/*.narrowPeak.gz".format(
        timepoint_dir))
    if len(timepoints_files) != len(atac_peak_files):
        parallel_copy(
            atac_peak_files,
            timepoint_dir)
    timepoints_files = glob.glob("{}/*.narrowPeak.gz".format(
        timepoint_dir))

    # -------------------------------------------
    # ANALYSIS 1 - generate master regions file
    # input: peak files
    # output: master peak file (BED)
    # -------------------------------------------
    logger.info("ANALYSIS: generate master regions file")
    master_regions_key = "atac.master.bed"
    out_data[master_regions_key] = '{0}/{1}.idr.master.bed.gz'.format(
        data_dir, prefix)
    if not os.path.isfile(out_data[master_regions_key]):
        merge_regions(timepoints_files, out_data[master_regions_key])

    # -------------------------------------------
    # ANALYSIS 2 - get read counts in these regions
    # input: master regions, read files (BED/tagAlign format)
    # output: matrix of read counts per region
    # -------------------------------------------
    logger.info("ANALYSIS: get read counts per region")
    adjustment = "ends"
    counts_key = "atac.counts.mat"
    out_data[counts_key] = '{0}/{1}.{2}.counts.mat.txt.gz'.format(
        data_dir, prefix, adjustment)
    if not os.path.isfile(out_data[counts_key]):
        # glob tagalign files
        atac_tagalign_files = sorted(
            glob.glob('{0}/{1}'.format(
                inputs["data_dir"],
                inputs["tagalign_glob"])))
        # remove files from media influenced timepoints
        for timepoint_string in args.inputs["params"]["media_timepoints"]:
            atac_tagalign_files = [
                filename for filename in atac_tagalign_files
                if timepoint_string not in filename]
        # make count matrix
        make_count_matrix(
            out_data[master_regions_key],
            atac_tagalign_files,
            out_data[counts_key],
            "ATAC",
            adjustment=adjustment,
            tmp_dir=results_dir)

    # -------------------------------------------
    # ANALYSIS 3 - run timeseries analysis on these regions
    # input: count matrix of regions
    # output: region trajectories
    # -------------------------------------------
    args = run_timeseries_workflow(
        args,
        prefix,
        datatype_key="atac",
        mat_key=counts_key)

    # extract dynamic BED and stable BED files
    # also cluster BED files

    # also set up stable mat for ATAC pooled

    # bioinformatics: GREAT, HOMER



    # now plot everything - pdfs, to put into illustrator

    

    return args
