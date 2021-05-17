"""Wrapper for utilizing tronn on datasets
"""

import os
import h5py
import glob
import logging

import numpy as np
import pandas as pd

from ggr.analyses.baselines import compare_sig_motifs
from ggr.analyses.baselines import compare_tf_motif_corrs
from ggr.analyses.bioinformatics import run_gprofiler
from ggr.analyses.linking import regions_to_genes
from ggr.analyses.nn import get_motif_region_set
from ggr.analyses.nn import get_sig_motif_list
from ggr.analyses.nn import compare_nn_gene_set_to_diff_expr
from ggr.analyses.nn import evaluate_enhancer_prediction
from ggr.analyses.nn import run_low_affinity_workflow
from ggr.analyses.nn import summarize_functional_enrichments
from ggr.analyses.nn import get_summit_centric_motif_maps
from ggr.analyses.nn import compare_aligned_motif_orientations
from ggr.analyses.nn import filter_summit_centric_maps
from ggr.analyses.nn import compare_homotypic_region_sets_to_chipseq
from ggr.analyses.tronn_scripts import tronn_intersect_with_expression_cmd

from ggr.util.utils import run_shell_cmd


def _run_motif_functional_enrichment_workflow(args, motif_name, bed_file, out_dir):
    """motif to region set and functional enrichments
    """
    # linking file
    links_dir = "{}/linking/proximity".format(args.outputs["results"]["dir"])
    #links_dir = "{}/linking/hichip".format(args.outputs["results"]["dir"])
    interactions_file = "{}/ggr.linking.ALL.overlap.interactions.txt.gz".format(links_dir)
        
    # link to nearby genes
    # TODO how are we up/downweighting things?
    genes_file = "{}/{}.genes.all.txt.gz".format(out_dir, os.path.basename(bed_file).split(".bed")[0])
    if not os.path.isfile(genes_file):
        regions_to_genes(
            bed_file,
            interactions_file,
            #args.outputs["annotations"]["tss.pc.bed"],
            args.outputs["results"]["rna"]["expression_filtering"]["tss.pc.expressed"],
            #args.outputs["data"]["tss.dynamic"],
            genes_file,
            filter_by_score=0.01)
            #filter_by_score=0.1) # proximity corr thresh

    # just look at top 1000?
    top_genes_file = "{}.max_1k.txt.gz".format(genes_file.split(".txt")[0])
    if not os.path.isfile(top_genes_file):
        genes_data = pd.read_csv(genes_file, sep="\t")
        if genes_data.shape[0] > 1000:
            genes_data = genes_data.iloc[:1000]
        genes_data.to_csv(top_genes_file, sep="\t", compression="gzip", index=False, header=True)
    genes_file = top_genes_file
    
    # then run gene enrichment
    gprofiler_output = "{}.go_gprofiler.txt".format(genes_file.split(".txt")[0])
    if not os.path.isfile(gprofiler_output):
        run_gprofiler(
            genes_file,
            args.outputs["data"]["rna.counts.pc.expressed.mat"],
            out_dir, ordered=True, header=True)
    
    return


def _run_motifs_functional_enrichment_workflow(args, pvals_file, nn_files, out_dir):
    """get a mapping of motifs to regions and functions
    """
    # first get motifs list
    motif_names = get_sig_motif_list(pvals_file, nn_files[0])

    # then get results per motif
    for motif_name in motif_names:
        bed_file = "{}/{}.bed.gz".format(out_dir, motif_name)
        if not os.path.isfile(bed_file):
            get_motif_region_set(nn_files, motif_name, bed_file)
        _run_motif_functional_enrichment_workflow(args, motif_name, bed_file, out_dir)

    # summarize results
    summary_dir = "{}/summary".format(out_dir)
    if not os.path.isdir(summary_dir):
        os.system("mkdir -p {}".format(summary_dir))
        summarize_functional_enrichments(out_dir, summary_dir)
    
    return


def _run_atac_summit_centric_workflow(args, pvals_file, nn_files, work_dir):
    """look at motif positions around atac summit
    """
    # set up
    enrichments_dir = "{}/enrichments.hits.homotypic_filt.genes_all".format(work_dir)
    if not os.path.isdir(enrichments_dir):
        os.system("mkdir -p {}".format(enrichments_dir))
    matrices_dir = "{}/data".format(work_dir)
    if not os.path.isdir(matrices_dir):
        os.system("mkdir -p {}".format(matrices_dir))
        
    # first get motifs list
    motif_names = get_sig_motif_list(pvals_file, nn_files[0])
    
    # then get results per motif
    for motif_name in motif_names:

        # debug
        if False:
            run_motifs = ["GRHL", "TP53", "CEBPA", "CEBPD", "TFAP2A"]
            run_motifs = ["FOX", "CEBP", "TP53", "KLF", "TFAP2", "MAF", "OVOL"]
            run_motifs = ["KLF"]
            dont_run = True
            for motif_str in run_motifs:
                if motif_str in motif_name:
                    dont_run = False
            if dont_run:
                continue
        
        # run FWD motifs        
        fwd_regions_file = "{}/{}.FWD.nn_filt.bed.gz".format(matrices_dir, motif_name)
        if not os.path.isfile(fwd_regions_file):
            get_summit_centric_motif_maps(
                motif_name,
                nn_files,
                args.outputs["data"]["atac.master.bed"],
                args.outputs["data"]["atac.master.summits.bed"],
                matrices_dir)

        # plotting - use these for final figs
        aligned_file = "{}/{}.FWD.aligned.hits.mat.txt.gz".format(matrices_dir, motif_name)
        plot_prefix = "{}/{}.FWD.aligned.hits".format(work_dir, motif_name)
        if not os.path.isfile("{}.sum.pdf".format(plot_prefix)):
            r_cmd = "~/git/ggr-project/R/plot.summit_aligned_motif_maps.R {} {} {}".format(
                aligned_file, plot_prefix, motif_name)
            print r_cmd
            os.system(r_cmd)

        aligned_file = "{}/{}.FWD.aligned.nn.mat.txt.gz".format(matrices_dir, motif_name)
        plot_prefix = "{}/{}.FWD.aligned.nn".format(work_dir, motif_name)
        if not os.path.isfile("{}.sum.pdf".format(plot_prefix)):
            r_cmd = "~/git/ggr-project/R/plot.summit_aligned_motif_maps.R {} {} {}".format(
                aligned_file, plot_prefix, motif_name)
            print r_cmd
            os.system(r_cmd)

        data_file = "{}/{}.FWD.affinities.hits.mat.txt.gz".format(matrices_dir, motif_name)
        plot_file = "{}/{}.FWD.num_v_affinity.hits.pdf".format(work_dir, motif_name)
        if not os.path.isfile(plot_file):
            r_cmd = "~/git/ggr-project/R/plot.num_motifs_vs_affinity.R {} {} {}".format(
                data_file, plot_file, motif_name)
            print r_cmd
            os.system(r_cmd)
            
        # run REV motifs - do not use, no asymmetry
        rev_regions_file = "{}/{}.REV.summit_center_filt.bed.gz".format(matrices_dir, motif_name)
        #if not os.path.isfile(rev_regions_file):
        if False:
            get_summit_centric_motif_maps(
                motif_name,
                nn_files,
                args.outputs["data"]["atac.master.bed"],
                args.outputs["data"]["atac.master.summits.bed"],
                matrices_dir,
                reverse_pwm=True)
            # do not use: if there were positional dependence based on orientation, should have seen it
            # when just looking at FWD version
            fwd_aligned_file = "{}/{}.FWD.aligned.mat.txt.gz".format(work_dir, motif_name) # note this is NN filt
            rev_aligned_file = "{}/{}.REV.aligned.mat.txt.gz".format(work_dir, motif_name) # note this is NN filt
            compare_aligned_motif_orientations(fwd_aligned_file, rev_aligned_file, work_dir, motif_name)

        # filter motif hits based on params: distance from summit, average affinity
        aligned_file = "{}/{}.FWD.aligned.hits.mat.txt.gz".format(matrices_dir, motif_name)
        homotypic_bed_file = "{}/{}.FWD.aligned.hits.homotypic_filt.bed.gz".format(enrichments_dir, motif_name)
        #aligned_file = "{}/{}.FWD.aligned.nn.mat.txt.gz".format(matrices_dir, motif_name)
        #homotypic_bed_file = "{}/{}.FWD.aligned.nn.homotypic_filt.bed.gz".format(enrichments_dir, motif_name)
        if not os.path.isfile(homotypic_bed_file):
            filter_summit_centric_maps(aligned_file, homotypic_bed_file, max_dist=50)

        # and then filter again against the NN filt bed file to get
        # regions with motifs of thresholded distance/affinity with NN importance
        homotypic_nn_filt_bed_file = "{}.nn_filt.bed.gz".format(homotypic_bed_file.split(".bed")[0])
        if not os.path.isfile(homotypic_nn_filt_bed_file):
            bedtools_cmd = "bedtools intersect -u -a {} -b {} | gzip -c > {}".format(
                homotypic_bed_file, fwd_regions_file, homotypic_nn_filt_bed_file)
            print bedtools_cmd
            os.system(bedtools_cmd)
            
        # how to compare back to ChIP-seq data?
        # q: how well does filtering for homotypic clusters improve prediction of ChIP-seq peaks?
        # i guess compare using the data in the nn files? compare motif hits to nn filt to homotypic rules
        # or do they show greater binding (agg plot)
        if False:
            compare_homotypic_region_sets_to_chipseq(
                work_dir, motif_name, nn_files)

        bed_file = homotypic_bed_file
        # run functional enrichments
        _run_motif_functional_enrichment_workflow(args, motif_name, bed_file, enrichments_dir)

    # then summarize functional enrichment results
    summary_dir = "{}/summary.hits.homotypic_filt.genes_all".format(work_dir)
    if not os.path.isdir(summary_dir):
        os.system("mkdir -p {}".format(summary_dir))
        summarize_functional_enrichments(enrichments_dir, summary_dir)
        
    return


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

    # results dir
    NN_DIR = "/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze"

    # -------------------------------------------
    # ANALYSIS - motifs to TF occupancy - low affinity sites
    # -------------------------------------------
    affinity_dir = "{}/low_affinity".format(results_dir)
    if not os.path.isdir(affinity_dir):
        os.system("mkdir -p {}".format(affinity_dir))
    prefix = "nn.low_affinity_hits"

    scanmotifs_dir = "/mnt/lab_data3/dskim89/ggr/nn/2020-01-13/scanmotifs/motifs.background.lite"
    #scanmotifs_dir = "/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze/motifs.input_x_grad.lite"
    scanmotifs_file = "{}/ggr.scanmotifs.h5".format(scanmotifs_dir)

    tfs = [
        #("TP53", 6, "TP63_LABELS", 1),
        ("TP53", 0, "TP63_LABELS", 0),
        ("CTCF", 6, "CTCF_LABELS", 1),
        ("KLF12", 6, "KLF4_LABELS", 1),
        ("TFAP2A", 6, "ZNF750_LABELS", 0)
    ]
    for motif_name, timepoint_idx, chipseq_key, chipseq_idx in tfs:
        plot_file = "{}/{}.{}.chipseq_pr_curves.pdf".format(
            affinity_dir, prefix, motif_name)
        if not os.path.isfile(plot_file):
            run_low_affinity_workflow(
                scanmotifs_file,
                motif_name,
                timepoint_idx,
                chipseq_key,
                chipseq_idx,
                affinity_dir,
                "{}.{}".format(prefix, motif_name))
    
    # -------------------------------------------
    # ANALYSIS - baselines: compare results to non-NN framework
    # -------------------------------------------
    baselines_dir = "{}/baselines".format(results_dir)
    if not os.path.isdir(baselines_dir):
        os.system("mkdir -p {}".format(baselines_dir))
    prefix = "nn_vs_homerATAC"

    # NN results - pvals file
    nn_sig_motifs_dir = "{}/motifs.sig".format(NN_DIR)
    nn_pvals_file = "{}/motifs.adjust.diff/pvals.h5".format(
        nn_sig_motifs_dir)
    
    # for first venn diagram, compare full homer results (not expr filt) with NN results
    plot_file = "{}/{}.homer_unfilt.compare_sig.plot.pdf".format(baselines_dir, prefix)
    if not os.path.isfile(plot_file):
        compare_sig_motifs(
            nn_pvals_file,
            args.outputs["results"]["atac"]["homer.sig_motifs.pvals"],
            ["NN", "Homer+ATAC"],
            baselines_dir, "{}.homer_unfilt".format(prefix))

    # switch to expression filt pvals file
    nn_pvals_file = "{}/motifs.adjust.diff.rna_filt.dmim/pvals.rna_filt.corr_filt.h5".format(
        nn_sig_motifs_dir)
        
    # 1) for venn diagram overlap, use the intersect file that matches correlation
    out_dir = "{}/corr_75".format(baselines_dir)
    if not os.path.isdir(out_dir):
        os.system("mkdir -p {}".format(out_dir))
    homer_pvals_file = "{}/motifs.rna_filt/pvals.rna_filt.corr_filt.h5".format(
        out_dir)
    if not os.path.isfile(homer_pvals_file):
        tronn_intersect_with_expression_cmd(
            args.outputs["results"]["atac"]["homer.sig_motifs.pvals"],
            args.outputs["annotations"]["pwms.renamed.nonredundant"],
            args.outputs["annotations"]["pwms.metadata.nonredundant.expressed"],
            "{}/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz".format(
                args.outputs["results"]["rna"]["timeseries"]["dir"]),
            out_dir)

    # plot
    plot_file = "{}/{}.compare_sig.plot.pdf".format(baselines_dir, prefix)
    if not os.path.isfile(plot_file):
        compare_sig_motifs(
            nn_pvals_file,
            homer_pvals_file,
            ["NN", "Homer+ATAC"],
            baselines_dir, prefix)
        
    # 2) for the correlation plot, use the intersect file that shows as many correlation values for homer+ATAC
    # TODO also for nn
    out_dir = "{}/corr_00.nn".format(baselines_dir)
    if not os.path.isdir(out_dir):
        os.system("mkdir -p {}".format(out_dir))
    nn_pvals_file = "{}/motifs.rna_filt/pvals.rna_filt.corr_filt.h5".format(
        out_dir)
    if not os.path.isfile(nn_pvals_file):
        tronn_intersect_with_expression_cmd(
            "{}/motifs.adjust.diff/pvals.h5".format(nn_sig_motifs_dir),
            args.outputs["annotations"]["pwms.renamed.nonredundant"],
            args.outputs["annotations"]["pwms.metadata.nonredundant.expressed"],
            "{}/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz".format(
                args.outputs["results"]["rna"]["timeseries"]["dir"]),
            out_dir,
            min_cor=0)
    
    out_dir = "{}/corr_00".format(baselines_dir)
    if not os.path.isdir(out_dir):
        os.system("mkdir -p {}".format(out_dir))
    homer_pvals_file = "{}/motifs.rna_filt/pvals.rna_filt.corr_filt.h5".format(
        out_dir)
    if not os.path.isfile(homer_pvals_file):
        tronn_intersect_with_expression_cmd(
            args.outputs["results"]["atac"]["homer.sig_motifs.pvals"],
            args.outputs["annotations"]["pwms.renamed.nonredundant"],
            args.outputs["annotations"]["pwms.metadata.nonredundant.expressed"],
            "{}/matrices/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.mat.txt.gz".format(
                args.outputs["results"]["rna"]["timeseries"]["dir"]),
            out_dir,
            min_cor=0)

    plot_file = "{}/{}.compare_corr.plot.pdf".format(baselines_dir, prefix)
    if not os.path.isfile(plot_file):
        compare_tf_motif_corrs(
            nn_pvals_file,
            homer_pvals_file,
            ["NN", "Homer+ATAC"],
            baselines_dir, prefix)

    # -------------------------------------------
    # ANALYSIS - predicting enhancers as defined by chromatin state
    # -------------------------------------------
    enhancer_dir = "{}/enh_state".format(results_dir)
    if not os.path.isdir(enhancer_dir):
        os.system("mkdir -p {}".format(enhancer_dir))
    
    evals_dir = "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03"
    model_prefix = "ggr.basset.clf.pretrained.folds.testfold-"
    eval_files = sorted(glob.glob(
        "{}/{}*/ggr.eval.h5".format(evals_dir, model_prefix)))

    plot_prefix = "{}/enh_state".format(enhancer_dir)
    if not os.path.isfile("{}.pdf".format(plot_prefix)):
        evaluate_enhancer_prediction(eval_files, plot_prefix)

    # -------------------------------------------
    # ANALYSIS - homotypic grammars
    # -------------------------------------------
    work_dir = "{}/homotypic".format(results_dir)
    if not os.path.isdir(work_dir):
        os.system("mkdir -p {}".format(work_dir))

    # sig motifs
    nn_pvals_file = (
        "{}/motifs.sig/motifs.adjust.diff.rna_filt.dmim/pvals.rna_filt.corr_filt.h5").format(
            NN_DIR)

    # scanmotif files
    NN_DIR = "/mnt/lab_data3/dskim89/ggr/nn/2019-03-12.freeze"
    nn_motif_files = [
        "{}/motifs.input_x_grad.early/ggr.scanmotifs.h5".format(NN_DIR),
        "{}/motifs.input_x_grad.mid/ggr.scanmotifs.h5".format(NN_DIR),
        "{}/motifs.input_x_grad.late/ggr.scanmotifs.h5".format(NN_DIR)]

    if False: # debug
        NN_DIR = "/mnt/lab_data3/dskim89/ggr/nn/2020-01-13/scanmotifs"
        nn_motif_files = [
            "{}/motifs.dynamic.early/ggr.scanmotifs.h5".format(NN_DIR),
            "{}/motifs.dynamic.mid/ggr.scanmotifs.h5".format(NN_DIR),
            "{}/motifs.dynamic.late/ggr.scanmotifs.h5".format(NN_DIR)
        ]

    if False:
        # run a functional enrichment mapping just based on presence of NN-active motifs
        _run_motifs_functional_enrichment_workflow(args, nn_pvals_file, nn_motif_files, work_dir)
    
    # do ATAC summit-centric view on motifs
    # also functional enrichment mapping here too
    if False:
        work_dir = "{}/homotypic.atac_summits.FIGS".format(results_dir)
        if not os.path.isdir(work_dir):
            os.system("mkdir -p {}".format(work_dir))
        
    _run_atac_summit_centric_workflow(args, nn_pvals_file, nn_motif_files, work_dir)

    # -------------------------------------------
    # ANALYSIS - compare homotypic v heterotypic
    # -------------------------------------------
    # numbers comparison of homotypic vs heterotypic
    work_dir = "{}/homotypic_v_heterotypic".format(results_dir)
    if not os.path.isdir(work_dir):
        os.system("mkdir -p {}".format(work_dir))

    # linking file
    links_dir = "{}/linking/proximity".format(args.outputs["results"]["dir"])
    interactions_file = "{}/ggr.linking.ALL.overlap.interactions.txt.gz".format(links_dir)

    # heterotypic
    heterotypic_bed_dir = "/mnt/lab_data3/dskim89/ggr/review/pre-review"
    heterotypic_bed_all = "{}/heterotypic.all_regions.bed.gz".format(work_dir)
    if not os.path.isfile(heterotypic_bed_all):
        copy_cmd = "cp {}/rules_all_regions.bed.gz {}".format(heterotypic_bed_dir, heterotypic_bed_all)
        print copy_cmd
        os.system(copy_cmd)

    heterotypic_genes_file = "{}/{}.genes.all.txt.gz".format(
        work_dir, os.path.basename(heterotypic_bed_all).split(".bed")[0])
    if not os.path.isfile(heterotypic_genes_file):
        regions_to_genes(
            heterotypic_bed_all,
            interactions_file,
            #args.outputs["annotations"]["tss.pc.bed"],
            #args.outputs["results"]["rna"]["expression_filtering"]["tss.pc.expressed"],
            args.outputs["data"]["tss.dynamic"],
            genes_file,
            filter_by_score=0.5)
    
    homotypic_bed_dir = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/results/learning/homotypic/enrichments.hits.homotypic_filt.genes_all"
    homotypic_bed_all = "{}/homotypic.all_regions.bed.gz".format(work_dir)
    if not os.path.isfile(homotypic_bed_all):
        bedtools_cmd = "zcat {}/*.homotypic_filt.bed.gz | sort -k1,1 -k2,2n | bedtools merge -i stdin | gzip -c > {}".format(
            homotypic_bed_dir, homotypic_bed_all)
        print bedtools_cmd
        os.system(bedtools_cmd)

    homotypic_genes_file = "{}/{}.genes.all.txt.gz".format(
        work_dir, os.path.basename(homotypic_bed_all).split(".bed")[0])
    if not os.path.isfile(homotypic_genes_file):
        regions_to_genes(
            homotypic_bed_all,
            interactions_file,
            #args.outputs["annotations"]["tss.pc.bed"],
            #args.outputs["results"]["rna"]["expression_filtering"]["tss.pc.expressed"],
            args.outputs["data"]["tss.dynamic"],
            genes_file,
            filter_by_score=0.5)

    total_genes_file = "{}/genes_all.txt.gz".format(work_dir)
    if not os.path.isfile(total_genes_file):
        heterotypic = pd.read_csv(heterotypic_genes_file, sep="\t")
        homotypic = pd.read_csv(homotypic_genes_file, sep="\t")
        total_genes = heterotypic.merge(homotypic, how="outer", on="gene_id")
        total_genes.to_csv(total_genes_file, sep="\t", compression="gzip")
    
    quit()
    
    # -------------------------------------------
    # ANALYSIS - predicting TF KD/KO gene sets from learning results
    # -------------------------------------------
    tf_kd_dir = "{}/tfs".format(results_dir)
    print tf_kd_dir
    if not os.path.isdir(tf_kd_dir):
        os.system("mkdir -p {}".format(tf_kd_dir))

    # motif inputs
    motifs = [
        ("TP53", 6, "/mnt/lab_data3/dskim89/labmembers/margaret/TP63_GSE33495_dex_results.csv"),
        ("TP53", 6, "/mnt/lab_data3/dskim89/labmembers/margaret/TP63_GSE67382_dex_results.csv"),
        ("KLF12", 6, "/mnt/lab_data3/dskim89/labmembers/margaret/KLF4_GSE32685_dex_results.csv"),
        #("KLF12", 6, "/mnt/lab_data3/dskim89/labmembers/margaret/KLF4_GSE111786_dex_results.csv"),
        #("KLF12", 6, "/mnt/lab_data3/dskim89/labmembers/margaret/KLF4_GSE140992_dex_results.csv"),
        ("TFAP2A", 6, "/mnt/lab_data3/dskim89/labmembers/margaret/ZNF750_GSE32685_dex_results.csv"),
        ("MAFB", 9, "/mnt/lab_data3/dskim89/labmembers/margaret/MAF_GSE52651_dex_results.csv"),
        #("FOXA1", 0, "/mnt/lab_data3/dskim89/labmembers/margaret/FOXC1_D0_GSE76647_dex_results.csv"),
        ("FOXA1", 9, "/mnt/lab_data3/dskim89/labmembers/margaret/FOXC1_D5_GSE76647_dex_results.csv")
    ]

    if False:
        # ATAC traj files
        atac_clusters_file = args.outputs["results"]["atac"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"]
        atac_mat_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.traj.mat"]
        atac_mat = pd.read_csv(atac_mat_file, sep="\t", index_col=0)
        atac_mat = atac_mat.drop("d05", axis=1)
        
        # RNA traj files
        rna_clusters_file = args.outputs["results"]["rna"]["timeseries"]["dp_gp"][
            "clusters.reproducible.hard.reordered.list"]
        rna_mat_file = args.outputs["data"][
            "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
        rna_mat = pd.read_csv(rna_mat_file, sep="\t", index_col=0)
    
    # linking file
    links_dir = "{}/linking/proximity".format(args.outputs["results"]["dir"])
    #links_dir = "{}/linking/hichip".format(args.outputs["results"]["dir"])
    interactions_file = "{}/ggr.linking.ALL.overlap.interactions.txt.gz".format(links_dir)

    # background atac file
    atac_background_file = args.outputs["data"]["atac.master.bed"]
    
    # nn files
    #nn_motif_files = glob.glob(
    #    "{}/motifs.input_x_grad.lite/ggr.scanmotifs.h5".format(NN_DIR))
    nn_motif_files = [
        "/mnt/lab_data3/dskim89/ggr/nn/2020-01-13/scanmotifs/motifs.background.lite/ggr.scanmotifs.h5"
    ]
    
    #if True:
    #    nn_motif_files = glob.glob(
    #        "{}/motifs.input_x_grad.mid/ggr.scanmotifs.h5".format(NN_DIR))

    if False:
        # run for each motif-timepoint
        for motif_name, timepoint_idx, diff_expr_file in motifs:
            print ">>", motif_name

            # first run with motif hits

            # get bed file of regions
            bed_file = "{}/{}-{}.bed.gz".format(tf_kd_dir, motif_name, timepoint_idx)
            #if not os.path.isfile(bed_file):
            if True:
                get_motif_region_set(
                    nn_motif_files,
                    motif_name,
                    bed_file,
                    timepoint_idx=timepoint_idx)

            # link to nearby genes
            # TODO how are we up/downweighting things?
            genes_file = "{}.genes.txt.gz".format(bed_file.split(".bed")[0])
            #if not os.path.isfile(genes_file):
            if True:
                regions_to_genes(
                    bed_file,
                    interactions_file,
                    args.outputs["annotations"]["tss.pc.bed"],
                    genes_file,
                    filter_by_score=0.00)
                    #filter_by_score=0.5) # proximity corr thresh

            # compare to differential dataset
            if True:
                compare_nn_gene_set_to_diff_expr(
                    genes_file,
                    diff_expr_file,
                    convert_table=args.outputs["annotations"]["geneids.mappings.mat"])


            # then run with NN-active motifs
            bed_file = "{}/{}-{}.bed.gz".format(tf_kd_dir, motif_name, timepoint_idx)
            #if not os.path.isfile(bed_file):
            if True:
                get_motif_region_set(
                    nn_motif_files,
                    motif_name,
                    bed_file,
                    timepoint_idx=timepoint_idx,
                    motif_key="sequence.active.pwm-scores.thresh.sum")

            # link to nearby genes
            # TODO how are we up/downweighting things?
            genes_file = "{}.genes.txt.gz".format(bed_file.split(".bed")[0])
            #if not os.path.isfile(genes_file):
            if True:
                regions_to_genes(
                    bed_file,
                    interactions_file,
                    args.outputs["annotations"]["tss.pc.bed"],
                    genes_file,
                    filter_by_score=0.00)
                    #filter_by_score=0.5) # proximity corr thresh

            # compare to differential dataset
            if True:
                compare_nn_gene_set_to_diff_expr(
                    genes_file,
                    diff_expr_file,
                    convert_table=args.outputs["annotations"]["geneids.mappings.mat"])
        
        
                    
    return args
