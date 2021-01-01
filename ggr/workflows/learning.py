"""Wrapper for utilizing tronn on datasets
"""

import os
import glob
import logging

from ggr.analyses.tronn_scripts import tronn_intersect_with_expression_cmd
from ggr.util.utils import run_shell_cmd


def runall(args, prefix):
    """all workflows for deep learning (TRONN)
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: run atac analyses")

    # set up data and results dir
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]

    results_dirname = "learning"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]
    
    # -------------------------------------------
    # ANALYSIS 0 - preprocess data to TRONN data formats
    # input: peak files
    # output: hdf5 files
    # -------------------------------------------
    dataset_dir = "{}/datasets/{}".format(results_dir, prefix)
    if not os.path.isdir(dataset_dir):
        
        datasets = []
        label_strings = []

        # TODO need full paths
        
        # atac labels
        atac_labels = sorted(
            glob.glob("{}/{}".format(
                args.inputs["atac"][args.cluster]["data_dir"],
                args.inputs["atac"][args.cluster]["idr_peak_glob"])))
        atac_labels = [os.path.abspath(filename) for filename in atac_labels]
        print len(atac_labels)
        label_strings.append("{}_LABELS={}".format(
            "ATAC", ",".join(atac_labels)))
        
        # trajectory labels
        traj_labels = sorted(
            glob.glob("./results/atac/timeseries/dp_gp/reproducible/hard/reordered/bed/*bed.gz"))
        traj_labels = [os.path.abspath(filename) for filename in traj_labels]
        print len(traj_labels)
        label_strings.append("{}_LABELS={}".format(
            "TRAJ", ",".join(traj_labels)))

        # histone labels
        histones = args.inputs["chipseq"][args.cluster]["histones"]["ordered_names"]
        for histone in histones:
            histone_labels = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["histones"][histone]["overlap_glob"])))
            histone_labels = [os.path.abspath(filename) for filename in histone_labels]
            print len(histone_labels)
            label_strings.append("{}_LABELS={}::method=histone_overlap".format(
                histone, ",".join(histone_labels)))

        # TF labels
        tfs = args.inputs["chipseq"][args.cluster]["tfs"]["ordered_names"]
        for tf in tfs:
            tf_labels = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["tfs"][tf]["idr_peak_glob"])))
            tf_labels = [os.path.abspath(filename) for filename in tf_labels]
            print len(tf_labels)
            label_strings.append("{}_LABELS={}".format(
                tf, ",".join(tf_labels)))

        # add other tf labels (external to GGR)
        tfs = args.inputs["chipseq"][args.cluster]["tfs"]["other"]["ordered_names"]
        for tf in tfs:
            tf_labels = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["tfs"]["other"]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["tfs"]["other"][tf]["narrowPeak_glob"])))
            tf_labels = [os.path.abspath(filename) for filename in tf_labels]
            print len(tf_labels)
            label_strings.append("{}_LABELS={}".format(
                tf, ",".join(tf_labels)))
        
        # chrom labels
        chrom_labels = sorted(
            glob.glob("./results/epigenome/dynamic/clusters/by_mark/bed/*bed.gz"))
        chrom_labels = [os.path.abspath(filename) for filename in chrom_labels]
        print len(chrom_labels)
        label_strings.append("{}_LABELS={}".format(
            "DYNAMIC_MARK", ",".join(chrom_labels)))

        chrom_labels = sorted(
            glob.glob("./results/epigenome/dynamic/clusters/by_state/bed/*bed.gz"))
        chrom_labels = [os.path.abspath(filename) for filename in chrom_labels]
        print len(chrom_labels)
        label_strings.append("{}_LABELS={}".format(
            "DYNAMIC_STATE", ",".join(chrom_labels)))

        chrom_labels = sorted(
            glob.glob("./results/epigenome/stable/clusters/by_mark/bed/*bed.gz"))
        chrom_labels = [os.path.abspath(filename) for filename in chrom_labels]
        print len(chrom_labels)
        label_strings.append("{}_LABELS={}".format(
            "STABLE_MARK", ",".join(chrom_labels)))

        chrom_labels = sorted(
            glob.glob("./results/epigenome/stable/clusters/by_state/bed/*bed.gz"))
        chrom_labels = [os.path.abspath(filename) for filename in chrom_labels]
        print len(chrom_labels)
        label_strings.append("{}_LABELS={}".format(
            "STABLE_STATE", ",".join(chrom_labels)))

        
        signal_strings = []
        
        # atac signals
        atac_signals = sorted(
            glob.glob("{}/{}".format(
                args.inputs["atac"][args.cluster]["data_dir"],
                args.inputs["atac"][args.cluster]["bigwig_pooled_glob"])))
        atac_signals = [os.path.abspath(filename) for filename in atac_signals]
        print len(atac_signals)
        signal_strings.append("{}_SIGNALS={}".format(
            "ATAC", ",".join(atac_signals)))

        # histone signals
        histones = args.inputs["chipseq"][args.cluster]["histones"]["ordered_names"]
        for histone in histones:
            histone_signals = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["histones"][histone]["pooled_bigwig_glob"])))
            histone_signals = [os.path.abspath(filename) for filename in histone_signals]
            signal_strings.append("{}_SIGNALS={}::window=1000".format(
                histone, ",".join(histone_signals)))

        # tf signals
        tfs = args.inputs["chipseq"][args.cluster]["tfs"]["ordered_names"]
        for tf in tfs:
            tf_signals = sorted(
                glob.glob("{}/{}".format(
                    args.inputs["chipseq"][args.cluster]["data_dir"],
                    args.inputs["chipseq"][args.cluster]["tfs"][tf]["pooled_bigwig_glob"])))
            tf_signals = [os.path.abspath(filename) for filename in tf_signals]
            signal_strings.append("{}_SIGNALS={}".format(
                tf, ",".join(tf_signals)))

        # run tronn preprocess
        preprocess = (
            "tronn preprocess \\\n"
            "--annotations {0} \\\n"
            "--labels \\\n \t{1} \\\n"
            "--signals \\\n \t{2} \\\n"
            "-o ./datasets/{3} \\\n"
            "--prefix {3} \\\n"
            "--master_label_keys ATAC_LABELS \\\n"
            "--genomewide \\\n"
            "--parallel {4}").format(
                args.inputs["annot"][args.cluster]["tronn_preprocess_annot"],
                " \\\n\t".join(label_strings),
                " \\\n\t".join(signal_strings),
                prefix,
                args.threads)

        run_shell_cmd("echo '{}' > {}/preprocess.tmp".format(preprocess, results_dir))
        #run_shell_cmd("source preprocess.tmp")

    # NOTE: in a linear run, the tronn processes (train, eval, interpret, etc) would all
    # be here in the code
        
    # run the baseline version of getting motifs/TFs using just homer and ATAC trajectory info
    out_dir = "{}/dp_gp/reproducible/hard/reordered/bed/homer_HOCOMOCO".format(
        args.outputs["results"]["atac"]["timeseries"]["dir"])
    if not os.path.isdir("{}/motifs.rna_filt".format(out_dir)):
        tronn_intersect_with_expression_cmd(
            args.outputs["results"]["atac"]["homer.sig_motifs.pvals"],
            args.outputs["annotations"]["pwms.renamed.nonredundant"],
            args.outputs["annotations"]["pwms.metadata.nonredundant.expressed"],
            "{}/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz".format(
                args.outputs["results"]["rna"]["timeseries"]["dir"]),
            out_dir)
                    
    return args
