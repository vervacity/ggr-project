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

from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import bed_to_gene_set_by_proximity
from ggr.analyses.linking import build_confusion_matrix


def runall(args, prefix):
    """all workflows for atac-seq data
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: run linking analyses")

    # set up data
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]
    
    results_dirname = "linking"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]

    # -------------------------------------------
    # ANALYSIS - overlap ATAC trajectories with RNA
    # input: atac trajectories, tss file
    # output: gene sets, enrichments, confusion matrix
    # -------------------------------------------
    logger.info("ANALYSIS: try naive proximal linking")

    # set up dir
    atac_linking_dir = "{}/atac".format(results_dir)
    run_shell_cmd("mkdir -p {}".format(atac_linking_dir))

    # use the TSS file that only has expressed genes
    tss_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"]["tss.w_clusters"]
    
    # generate the gene sets for each ATAC trajectory
    atac_traj_dir = "{}/reproducible/hard/reordered/bed".format(
        args.outputs["results"]["atac"]["timeseries"]["dp_gp"]["dir"])
    traj_bed_files = sorted(glob.glob("{}/*cluster*bed.gz".format(atac_traj_dir)))
    gene_set_files = []
    for traj_bed_file in traj_bed_files:
        gene_set_file = "{}/{}.linked_genes.txt.gz".format(
            atac_linking_dir,
            os.path.basename(traj_bed_file).split(".bed")[0])
        if not os.path.isfile(gene_set_file):
            bed_to_gene_set_by_proximity(
                traj_bed_file,
                tss_file,
                gene_set_file,
                k_nearest=2,
                max_dist=100000)
        gene_set_files.append(gene_set_file)
        
    # run enrichments (gprofiler)
    background_gene_set_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    atac_linked_genes_enrich_dir = "{}/enrichments".format(atac_linking_dir)
    if not os.path.isdir(atac_linked_genes_enrich_dir):
        run_shell_cmd("mkdir -p {}".format(atac_linked_genes_enrich_dir))
        for gene_set_file in gene_set_files:
            run_gprofiler(
                gene_set_file,
                background_gene_set_file,
                atac_linked_genes_enrich_dir)
    
    # collect into a confusion matrix
    cluster_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"]["clusters.reproducible.hard.reordered.list"]
    confusion_matrix_file = "{}/{}.confusion_mat.txt.gz".format(atac_linking_dir, prefix)
    #if not os.path.isfile(confusion_matrix_file):
    if True:
        build_confusion_matrix(
            traj_bed_files,
            gene_set_files,
            cluster_file,
            confusion_matrix_file)



    # and then rerun select enrichments?
    

    quit()
    
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
    if not os.path.isdir("{}/homer".format(cluster_bed_dir)):
        for bed_file in bed_files:
            run_bioinformatics_on_bed(
                bed_file,
                background_bed_file,
                cluster_bed_dir)

    return args
