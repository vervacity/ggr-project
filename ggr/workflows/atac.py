"""ggr workflow for atac analyses
"""

import os
import glob
import signal
import logging

import pandas as pd

from ggr.util.utils import run_shell_cmd
from ggr.util.utils import parallel_copy

from ggr.util.bed_utils import merge_regions
from ggr.util.bed_utils import id_to_bed

from ggr.analyses.atac import compare_to_scATAC
from ggr.analyses.atac import get_consensus_summits_file
from ggr.analyses.filtering import filter_for_ids
from ggr.analyses.counting import make_count_matrix
from ggr.analyses.counting import count_present_regions_per_sample
from ggr.analyses.utils import build_id_matching_mat
from ggr.analyses.utils import plot_PCA

from ggr.workflows.timeseries import run_timeseries_workflow

from ggr.analyses.bioinformatics import run_bioinformatics_on_bed
from ggr.analyses.bioinformatics import aggregate_homer_results_h5


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
    timepoints_files = sorted(
        glob.glob("{}/*.narrowPeak.gz".format(
            timepoint_dir)))
    if len(timepoints_files) != len(atac_peak_files):
        parallel_copy(
            atac_peak_files,
            timepoint_dir,
            num_threads=args.threads)
    timepoints_files = sorted(
        glob.glob("{}/*.narrowPeak.gz".format(
            timepoint_dir)))
    args.outputs["results"]["label_dirs"].append(timepoint_dir)

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

    # generate negatives (opposite)
    negatives_bed = "{}/{}.idr.negatives.bed.gz".format(
        data_dir, prefix)
    if not os.path.isfile(negatives_bed):
        os.system("cat {} | sort -k1,1 -k2,2n > hg19.chromsizes".format(
            args.inputs["annot"][args.cluster]["chromsizes"]))
        complement = (
            "bedtools complement -i {} -g {} | gzip -c > {}").format(
                out_data[master_regions_key],
                "hg19.chromsizes",
                negatives_bed)
        run_shell_cmd(complement)
        os.system("rm hg19.chromsizes")

    # figure out how many regions are active per timepoint
    plot_dir = "{}/plots".format(results_dir)
    static_region_summary_file = "{}/plots/{}.idr.timepoint_count_summary.static.txt".format(
        results_dir, prefix)
    if not os.path.isfile(static_region_summary_file):
        os.system("mkdir -p {}".format(plot_dir))
        count_present_regions_per_sample(
            out_data[master_regions_key],
            timepoints_files,
            static_region_summary_file,
            assay="ATAC")
    plot_file = "{}.pdf".format(static_region_summary_file.split(".txt")[0])
    if not os.path.isfile(plot_file):
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.counts.global.R {} {} static".format(
            static_region_summary_file, plot_file)
        print plot_cmd
        os.system(plot_cmd)

    # -------------------------------------------
    # ANALYSIS - make an ATAC summits file
    # input: master regions, bed files
    # output: consensus summits file
    # -------------------------------------------
    summits_key = "atac.master.summits.bed"
    out_data[summits_key] = '{0}/{1}.idr.summits.bed.gz'.format(
        data_dir, prefix)
    if not os.path.isfile(out_data[summits_key]):
        get_consensus_summits_file(
            out_data[master_regions_key],
            timepoints_files,
            out_data[summits_key])
        
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
            tmp_dir=results_dir,
            parallel=args.threads)

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

    # plot dynamic summary: filter for just d0
    differential_summary = "{}/timeseries/deseq2/differential_summary.txt.gz".format(results_dir)
    d0_baseline_file = "{}.d0_baseline_only.txt.gz".format(differential_summary.split(".txt")[0])
    if not os.path.isfile(d0_baseline_file):
        header_cmd = 'echo "timepoint\tup\tdown" | gzip -c > {}'.format(d0_baseline_file)
        os.system(header_cmd)
        filter_cmd = 'zcat {} | grep -e "^d00" | gzip -c >> {}'.format(
            differential_summary, d0_baseline_file)
        print filter_cmd
        os.system(filter_cmd)
    plot_file = "{}/plots/{}.idr.timepoint_count_summary.dynamic.pdf".format(
        results_dir, prefix)
    if not os.path.isfile(plot_file):
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.counts.global.R {} {} dynamic".format(
            d0_baseline_file, plot_file)
        print plot_cmd
        os.system(plot_cmd)

    # sequential
    differential_summary = "{}/timeseries/deseq2.sequential/differential_summary.txt.gz".format(results_dir)
    plot_file = "{}/plots/{}.idr.timepoint_count_summary.dynamic.sequential.pdf".format(
        results_dir, prefix)
    if not os.path.isfile(plot_file):
        plot_cmd = "Rscript ~/git/ggr-project/R/plot.counts.global.R {} {} dynamic".format(
            differential_summary, plot_file)
        print plot_cmd
        os.system(plot_cmd)
    
    # filter for reproducible dynamic IDs and save to data dir
    dynamic_mat_key = "atac.counts.pooled.rlog.dynamic.mat"
    dynamic_traj_mat_key = "{}.traj.mat".format(dynamic_mat_key.split(".mat")[0])
    reproducible_dynamic_file = "{}/{}.traj.mat.txt.gz".format(
        data_dir,
        os.path.basename(args.outputs["data"][dynamic_mat_key]).split(".mat")[0])
    args.outputs["data"][dynamic_traj_mat_key] = reproducible_dynamic_file
    if not os.path.isfile(args.outputs["data"][dynamic_traj_mat_key]):
        tmp_id_file = "ids.tmp.txt"
        pd.read_csv(
            args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
                "clusters.reproducible.hard.reordered.list"], sep="\t")["id"].to_csv(
                    tmp_id_file, index=False, sep="\t")
        filter_for_ids(
            args.outputs["data"][dynamic_mat_key],
            tmp_id_file,
            args.outputs["data"][dynamic_traj_mat_key])
        os.system("rm {}".format(tmp_id_file))
    
    # get PCA
    plot_dir = "{}/plots".format(
        out_results["timeseries"]["dir"])
    if not os.path.isdir(plot_dir):
        run_shell_cmd("mkdir -p {}".format(plot_dir))

        # pull the two reps
        atac_mat_files = sorted(glob.glob(
            "{}/ggr.atac.*.rep*rlog.dynamic.mat.txt.gz".format(
                args.outputs["data"]["dir"])))

        # and build filtered mat files based on final dynamic set
        final_dynamic_list = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"]
        filt_mat_files = []
        for atac_mat_file in atac_mat_files:
            filt_atac_mat_file = "{}/{}.filt.mat.txt.gz".format(
                plot_dir,
                os.path.basename(atac_mat_file).split(".mat")[0])
            build_id_matching_mat(
                final_dynamic_list,
                atac_mat_file,
                filt_atac_mat_file,
                keep_cols=[],
                primary_id_col="id")
            filt_mat_files.append(filt_atac_mat_file)

        # plot PCA/correlation (bio reps separately and pooled)
        pca_file = "{}/{}.pca.pdf".format(
            plot_dir,
            os.path.basename(atac_mat_files[0]).split(".rep")[0])
        plot_PCA(filt_mat_files, pca_file)
    
    # get BED files for each cluster and add to label_dirs
    cluster_dir = "{}/timeseries/dp_gp/reproducible/hard/reordered".format(results_dir)
    cluster_bed_dir = "{}/bed".format(cluster_dir)
    if not os.path.isdir(cluster_bed_dir):
        run_shell_cmd("mkdir -p {}".format(cluster_bed_dir))
        id_files = glob.glob("{}/*cluster_*txt.gz".format(cluster_dir))
        for id_file in id_files:
            bed_file = "{}/{}.bed.gz".format(
                cluster_bed_dir,
                os.path.basename(id_file).split(".txt")[0])
            id_to_bed(
                id_file,
                bed_file,
                sort=True)
    args.outputs["results"]["label_dirs"].append(cluster_bed_dir)
            
    # -------------------------------------------
    # ANALYSIS 4 - get stable and dynamic BED files
    # input: dynamic ids
    # output: dynamic and stable BED files
    # -------------------------------------------
    logger.info("ANALYSIS: split dynamic and stable BEDs")
    rlog_mat_key = "atac.counts.pooled.rlog.mat"
    
    dynamic_mat_key = "atac.counts.pooled.rlog.dynamic.mat"
    dynamic_bed_key = "atac.dynamic.bed"
    out_data[dynamic_bed_key] = "{}.bed.gz".format(
        out_data[dynamic_mat_key].split(".mat")[0])
    if not os.path.isfile(out_data[dynamic_bed_key]):
        id_to_bed(
            out_data[dynamic_mat_key],
            out_data[dynamic_bed_key],
            sort=True)
    
    stable_mat_key = "atac.counts.pooled.rlog.stable.mat"
    out_data[stable_mat_key] = "{}.stable.mat.txt.gz".format(
        out_data[dynamic_mat_key].split(".dynamic")[0])
    if not os.path.isfile(out_data[stable_mat_key]):
        filter_for_ids(
            out_data[rlog_mat_key],
            out_results["timeseries"]["dynamic_ids.list"],
            out_data[stable_mat_key],
            opposite=True)
    
    stable_bed_key = "atac.stable.bed"
    out_data[stable_bed_key] = "{}.bed.gz".format(
        out_data[stable_mat_key].split(".mat")[0])
    if not os.path.isfile(out_data[stable_bed_key]):
        id_to_bed(
            out_data[stable_mat_key],
            out_data[stable_bed_key],
            sort=True)

    # -------------------------------------------
    # ANALYSIS 5 - bioinformatics
    # input: bed dir
    # output: HOMER and GREAT results
    # -------------------------------------------
    logger.info("ANALYSIS: run HOMER/GREAT on clusters")
    bed_files = glob.glob("{}/*bed.gz".format(cluster_bed_dir))
    background_bed_file = out_data[master_regions_key]
    if not os.path.isdir("{}/homer_HOCOMOCO".format(cluster_bed_dir)):
        for bed_file in bed_files:
            run_bioinformatics_on_bed(
                bed_file,
                background_bed_file,
                cluster_bed_dir,
                mknown=args.outputs["annotations"]["pwms.renamed.nonredundant.homer"],
                mknown_name="HOCOMOCO")

    # make a pvals file of ATAC trajectory results that can be run through tronn script
    pval_file = "{}/homer_HOCOMOCO/pvals.h5".format(cluster_bed_dir)
    out_results["homer.sig_motifs.pvals"] = pval_file
    if not os.path.isfile(pval_file):
        aggregate_homer_results_h5(
            args.outputs["annotations"]["pwms.renamed.nonredundant"],
            "{}/homer_HOCOMOCO".format(cluster_bed_dir),
            pval_file)

    # compare results to scATAC
    compare_dir = "{}/sc_vs_bulk".format(results_dir)
    if not os.path.isdir(compare_dir):
        os.system("mkdir -p {}".format(compare_dir))
    scATAC_file = "/mnt/lab_data3/dskim89/ggr/data_external/GSE116248_scATAC/GSE116248_Peak_counts_Keratinocyte_WT.txt.gz"

    if False:
        bulk_read_files = sorted(glob.glob(
            "{}/{}".format(
                args.inputs["atac"][args.cluster]["data_dir"],
                args.inputs["atac"][args.cluster]["tagalign_glob"])))
        compare_to_scATAC(
            scATAC_file,
            bulk_read_files,
            compare_dir,
            threads=args.threads)
        
    return args
