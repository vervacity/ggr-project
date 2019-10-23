"""description: workflow for inference downstream
"""

import os
import re
import h5py
import glob
import logging

import numpy as np
import pandas as pd

from ggr.analyses.linking import regions_to_genes_w_correlation_filtering
from ggr.analyses.bioinformatics import run_gprofiler
from ggr.util.utils import run_shell_cmd
from ggr.util.go_utils import GoParams
from ggr.util.go_utils import is_enriched

from tronn.util.utils import DataKeys
from tronn.util.formats import array_to_bed
from tronn.interpretation.syntax import analyze_multiplicity


def run_multiplicity_workflow(args, prefix, sig_pwms_file):
    """for sig pwms, look at multiplicity
    """
    _MAX_COUNT = 6
    _MIN_HIT_REGIONS = 50
    OUT_DIR = args.outputs["results"]["inference"]["dir"]
    TMP_DIR = "{}/tmp".format(
        args.outputs["results"]["inference"]["dir"])
    os.system("mkdir -p {}".format(TMP_DIR))
            
    # set up motifs files being used
    motifs_files = {}
    motifs_files["early"] = "{}/{}/ggr.scanmotifs.h5".format(
        args.inputs["inference"][args.cluster]["scanmotifs_dir"],
        args.inputs["inference"][args.cluster]["scanmotifs_early_dir"])
    motifs_files["mid"] = "{}/{}/ggr.scanmotifs.h5".format(
        args.inputs["inference"][args.cluster]["scanmotifs_dir"],
        args.inputs["inference"][args.cluster]["scanmotifs_mid_dir"])
    motifs_files["late"] = "{}/{}/ggr.scanmotifs.h5".format(
        args.inputs["inference"][args.cluster]["scanmotifs_dir"],
        args.inputs["inference"][args.cluster]["scanmotifs_late_dir"])
    
    # read in the sig pwms
    sig_pwms = list(pd.read_csv(
        sig_pwms_file, sep="\t", header=0, index_col=0).index.values)
    logging.info("{} sig pwms found".format(len(sig_pwms)))
    sig_pwms_trajs = pd.read_csv(sig_pwms_file, sep="\t", header=0, index_col=0)
    
    # read in the list of all pwm names
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    with h5py.File(motifs_files["early"], "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    num_pwms = len(all_pwms)
    print num_pwms

    # traj to group
    traj_to_group = args.inputs["inference"]["traj_to_group"]

    # read in dynamic genes
    filter_genes = pd.read_table(
        args.outputs["data"][
            "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"],
        index_col=0).index.values.tolist()

    # read in rna signals
    rna_signal_file = args.outputs["data"]["rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
    rna_signal_mat = pd.read_csv(rna_signal_file, sep="\t", index_col=0)
    rna_signal_mat[:] = np.subtract(
        rna_signal_mat.values,
        np.expand_dims(rna_signal_mat.values[:,0], axis=-1))

    # read in region signals
    atac_signal_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.traj.mat"]
    region_signal_mat = pd.read_csv(atac_signal_file, sep="\t", index_col=0)
    region_signal_mat[:] = np.subtract(
        region_signal_mat.values,
        np.expand_dims(region_signal_mat.values[:,0], axis=-1))
    region_signal_mat = region_signal_mat.drop("d05", axis=1)

    # other files
    background_rna_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    links_file = args.outputs["results"]["linking"]["links.proximity"]
    tss_file = args.outputs["annotations"]["tss.pc.bed"]
    
    # results arrays
    all_results = {}
    signal_keys = args.inputs["inference"][args.cluster]["signal_keys"]
    for key in signal_keys:
        all_results[key] = []
    all_results["filt"] = []
        
    # analyze each pwm
    keep_pwm_names = [] # attach to results h5 for plotting
    keep_sig_pwm_names = [] # use for filtered sig pwms pattern file
    for pwm_idx in range(len(sig_pwms)):
        
        # name and global index
        pwm_name = sig_pwms[pwm_idx]
        pwm_name_clean = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name_clean = re.sub(".UNK.0.A", "", pwm_name_clean)
        pwm_global_idx = np.where(
            [1 if pwm_name in global_name else 0
             for global_name in all_pwms])[0][0]
        print pwm_name_clean, pwm_global_idx
        
        # reverse is in second half of indices
        rc_idx = pwm_global_idx + num_pwms
        pwm_indices = [pwm_global_idx, rc_idx]
        print pwm_indices

        # figure out which motif files to actually load
        keys = []
        trajs = sig_pwms_trajs.loc[pwm_name]
        trajs = trajs[trajs != 0].index.values.tolist()
        keys = list(set([traj_to_group[traj] for traj in trajs]))
        sig_motifs_files = [motifs_files[key] for key in keys]

        # TODO consider: loading traj mutatemotifs? would then need to pull the correct
        # index to find the right mutate results
        
        # get multiplicity
        # filtering: looking for
        results = analyze_multiplicity(
            sig_motifs_files, pwm_indices, max_count=_MAX_COUNT, solo_filter=True)

        filt = 0
        if results is not None:
            keep_pwm_names.append(pwm_name_clean)
        
            # save to summary array
            for key in signal_keys:
                all_results[key].append(results[key]["count"])
        
            # per level, get region set and do gene set enrichments
            hits_per_region = results["hits_per_region"]
            max_enriched_thresh = 0
            for count_thresh in range(1, _MAX_COUNT+1):
                # get region ids
                thresholded = hits_per_region[hits_per_region["hits"] >= count_thresh]
                thresholded_metadata = thresholded.index.values
                
                # check how many regions, do not continue if not enough regions
                if thresholded_metadata.shape[0] < _MIN_HIT_REGIONS:
                    continue

                # convert to BED
                tmp_bed_file = "{}/{}.regions.count_thresh-{}.bed.gz".format(
                    TMP_DIR, pwm_name_clean, count_thresh)
                if not os.path.isfile(tmp_bed_file):
                    array_to_bed(thresholded_metadata, tmp_bed_file, merge=False)

                # get linked genes
                tmp_genes_file = "{}/{}.linked_genes.count_thresh-{}.txt.gz".format(
                    TMP_DIR, pwm_name_clean, count_thresh)
                if not os.path.isfile(tmp_genes_file):
                    regions_to_genes_w_correlation_filtering(
                        tmp_bed_file,
                        links_file,
                        tss_file,
                        tmp_genes_file,
                        region_signal_mat,
                        rna_signal_mat,
                        filter_by_score=0.5,
                        filter_genes=filter_genes,
                        corr_thresh=0,
                        pval_thresh=1)

                # do not continue if no linked genes
                linked_genes = pd.read_csv(tmp_genes_file, sep="\t", header=0)
                if linked_genes.shape[0] == 0:
                    continue

                # run enrichment calculation
                enrichment_file = "{}/{}.linked_genes.count_thresh-{}.go_gprofiler.txt".format(
                    TMP_DIR, pwm_name_clean, count_thresh)
                if not os.path.isfile(enrichment_file):
                    run_gprofiler(
                        tmp_genes_file,
                        background_rna_file,
                        TMP_DIR,
                        ordered=True)

                # check if any enrichment, if not then continue
                if not is_enriched(enrichment_file):
                    continue

                # if passed all of this, move results to not tmp
                os.system("cp {} {}".format(
                    enrichment_file, OUT_DIR))

                # update important variables
                filt = 1
                max_enriched_thresh = count_thresh

            if max_enriched_thresh >= 2:
                keep_sig_pwm_names.append(pwm_name)
                
        # add in filt
        all_results["filt"].append(filt)

        # check in
        print all_results["filt"], np.sum(all_results["filt"])
        print keep_sig_pwm_names, len(keep_sig_pwm_names)
        
    # save out summary array
    h5_results_file = "{}/genome.homotypic.multiplicity_v_activity.h5".format(OUT_DIR)
    if not os.path.isfile(h5_results_file):
        filt = np.array(all_results["filt"])
        for key in signal_keys:
            summary_results = np.stack(all_results[key], axis=0)
            with h5py.File(h5_results_file, "a") as out:
                out_key = "{}/counts".format(key)
                out.create_dataset(out_key, data=summary_results)
                out[out_key].attrs["pwm_names"] = keep_pwm_names
                
                out_key = "{}/filt".format(key)
                out.create_dataset(out_key, data=filt)
                out[out_key].attrs["pwm_names"] = keep_pwm_names

    # save out new sig pwms file
    new_sig_pwms_file = "{}/{}.multiplicity_filt.txt.gz".format(
        OUT_DIR,
        os.path.basename(sig_pwms_file).split(".txt")[0])
    new_sig_pwms = sig_pwms_trajs.loc[keep_sig_pwm_names]
    print new_sig_pwms
    new_sig_pwms.to_csv(new_sig_pwms_file, sep="\t", compression="gzip")
    
    # and plot
    plot_cmd = "{}/plot.results.multiplicity.summary.R {} {}/genome.homotypic.multiplicity TRUE {}".format(
        "~/git/ggr-project/figs/fig_4.homotypic", h5_results_file,
        args.outputs["results"]["inference"]["dir"],
        " ".join(signal_keys))
    
    return new_sig_pwms_file


def run_syntax_workflow(args, prefix, sig_pwms_file):
    """for sig pwms, look at multiplicity
    """
    OUT_DIR = args.outputs["results"]["inference"]["dir"]
    TMP_DIR = "{}/tmp".format(
        args.outputs["results"]["inference"]["dir"])
    os.system("mkdir -p {}".format(TMP_DIR))
            
    # set up motifs files being used
    motifs_files = {}
    motifs_files["early"] = "{}/{}/ggr.scanmotifs.h5".format(
        args.inputs["inference"][args.cluster]["scanmotifs_dir"],
        args.inputs["inference"][args.cluster]["scanmotifs_early_dir"])
    motifs_files["mid"] = "{}/{}/ggr.scanmotifs.h5".format(
        args.inputs["inference"][args.cluster]["scanmotifs_dir"],
        args.inputs["inference"][args.cluster]["scanmotifs_mid_dir"])
    motifs_files["late"] = "{}/{}/ggr.scanmotifs.h5".format(
        args.inputs["inference"][args.cluster]["scanmotifs_dir"],
        args.inputs["inference"][args.cluster]["scanmotifs_late_dir"])
    
    # read in the sig pwms
    sig_pwms = list(pd.read_csv(
        sig_pwms_file, sep="\t", header=0, index_col=0).index.values)
    logging.info("{} sig pwms found".format(len(sig_pwms)))
    sig_pwms_trajs = pd.read_csv(sig_pwms_file, sep="\t", header=0, index_col=0)
    
    # read in the list of all pwm names
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    with h5py.File(motifs_files["early"], "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    num_pwms = len(all_pwms)
    print num_pwms

    # traj to group
    traj_to_group = args.inputs["inference"]["traj_to_group"]

    # read in dynamic genes
    filter_genes = pd.read_table(
        args.outputs["data"][
            "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"],
        index_col=0).index.values.tolist()

    # read in rna signals
    rna_signal_file = args.outputs["data"]["rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
    rna_signal_mat = pd.read_csv(rna_signal_file, sep="\t", index_col=0)
    rna_signal_mat[:] = np.subtract(
        rna_signal_mat.values,
        np.expand_dims(rna_signal_mat.values[:,0], axis=-1))

    # read in region signals
    atac_signal_file = args.outputs["data"]["atac.counts.pooled.rlog.dynamic.traj.mat"]
    region_signal_mat = pd.read_csv(atac_signal_file, sep="\t", index_col=0)
    region_signal_mat[:] = np.subtract(
        region_signal_mat.values,
        np.expand_dims(region_signal_mat.values[:,0], axis=-1))
    region_signal_mat = region_signal_mat.drop("d05", axis=1)

    # other files
    background_rna_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    links_file = args.outputs["results"]["linking"]["links.proximity"]
    tss_file = args.outputs["annotations"]["tss.pc.bed"]
    
    # results arrays
    all_results = {}
    signal_keys = args.inputs["inference"][args.cluster]["signal_keys"]
    for key in signal_keys:
        all_results[key] = []
    all_results["filt"] = []
        
    # analyze each pwm
    keep_pwm_names = [] # attach to results h5 for plotting
    keep_sig_pwm_names = [] # use for filtered sig pwms pattern file
    for pwm_idx in range(len(sig_pwms)):
        
        # name and global index
        pwm_name = sig_pwms[pwm_idx]
        pwm_name_clean = re.sub("HCLUST-\\d+_", "", pwm_name)
        pwm_name_clean = re.sub(".UNK.0.A", "", pwm_name_clean)
        pwm_global_idx = np.where(
            [1 if pwm_name in global_name else 0
             for global_name in all_pwms])[0][0]
        print pwm_name_clean, pwm_global_idx
        
        # reverse is in second half of indices
        rc_idx = pwm_global_idx + num_pwms
        pwm_indices = [pwm_global_idx, rc_idx]
        print pwm_indices

        # figure out which motif files to actually load
        keys = []
        trajs = sig_pwms_trajs.loc[pwm_name]
        trajs = trajs[trajs != 0].index.values.tolist()
        keys = list(set([traj_to_group[traj] for traj in trajs]))
        sig_motifs_files = [motifs_files[key] for key in keys]

        # TODO consider: loading traj mutatemotifs? would then need to pull the correct
        # index to find the right mutate results
        
        # go through orientations
        possible_orientations = ["FWD", "REV"]
        for orientation_a in possible_orientations:
            for orientation_b in possible_orientations:
                # set up new pwm indices
                pwm_indices = []
                
                # figure out whether fwd or rev
                if orientation_a == "FWD":
                    pwm_indices.append(pwm_global_idx)
                elif orientation_a == "REV":
                    pwm_indices.append(rc_idx)
                
                # figure out whether fwd or rev
                if orientation_b == "FWD":
                    pwm_indices.append(pwm_global_idx)
                elif orientation_b == "REV":
                    pwm_indices.append(rc_idx)
                    
                assert len(pwm_indices) == 2
                results_key = "{}_{}".format(orientation_a, orientation_b)
                print results_key, pwm_indices
                
                # build aligned array around first pwm
                aligned_results = analyze_syntax(
                    motifs_files, pwm_indices, solo_filter=True)
                
                if aligned_results is None:
                    continue

                if orientation_b == "FWD":
                    results[results_key] = aligned_results[pwm_global_idx]
                elif orientation_b == "REV":
                    results[results_key] = aligned_results[rc_idx]

        # 1) AGNOSTIC:  both(anchor)_x_both - merge all
        # 2) AGNOS_IN: both(anchor)_x_in - fwd_REV(+) + rev_REV(+) + fwd_FWD(-) + rev_FWD(-)
        # 3) AGNOS_OUT: both(anchor)_x_out - fwd_FWD(+) + rev_FWD(+) + fwd_REV(-) + rev_REV(-)
        # 4) SAME_DIR: same_x_same - fwd_FWD + rev_REV
        # 5) IN: in_x_in - fwd_REV(+) + rev_FWD(-)
        # 6) OUT: out_x_out - rev_FWD(+) + fwd_REV(-)
        # NOTE: ignoring rev_BOTH and fwd_BOTH - same as 2,3
        adjustments = {
            "BOTH_BOTH": (
                ["FWD_FWD", "FWD_REV", "REV_FWD", "REV_REV"],
                [".", ".", ".", "."]),
            "BOTH_IN": (
                ["FWD_FWD", "FWD_REV", "REV_FWD", "REV_REV"],
                ["-", "+", "-", "+"]),
            "BOTH_OUT": (
                ["FWD_FWD", "FWD_REV", "REV_FWD", "REV_REV"],
                ["+", "-", "+", "-"]),
            "SAME": (
                ["FWD_FWD", "REV_REV"],
                [".", "."]),
            "IN": (
                ["FWD_REV", "REV_FWD"],
                ["+", "-"]),
            "OUT": (
                ["FWD_REV", "REV_FWD"],
                ["-", "+"])}
        
        results_adjusted = {}
        enrichments_summary = None
        for adj_key in sorted(adjustments.keys()):
            print adj_key
            results_adjusted[adj_key] = recombine_syntax_results(
                results,
                adjustments[adj_key][0],
                adjustments[adj_key][1],
                signal_keys)

            if len(results_adjusted[adj_key].keys()) == 0:
                continue
            
            # debug
            print adj_key, results_adjusted[adj_key]["scores"].shape[0]
            
            # plot results
            plot_prefix = "{}/{}.{}".format(
                OUT_DIR, pwm_name_clean, adj_key)

            # plot pwm scores
            score_spacing_distr_file = "{}.genome.active_pwm_scores.avg.txt.gz".format(
                plot_prefix)
            aligned_scores = results_adjusted[adj_key]["scores"]
            num_examples = aligned_scores.shape[0]
            aligned_scores = aligned_scores.mean(axis=0).transpose()
            aligned_scores.loc[0] = 0
            aligned_scores = aligned_scores.divide(aligned_scores.sum()).reset_index()
            aligned_scores.columns = ["position", "active"]
            aligned_scores.to_csv(
                score_spacing_distr_file,
                sep="\t", header=True, index=False, compression="gzip")
            plot_cmd = "{}/plot.spacing.freq.indiv.R {} {} {}".format(
                SCRIPT_DIR, score_spacing_distr_file, num_examples, plot_prefix)
            print plot_cmd
            os.system(plot_cmd)

            # plot signals
            for signal_key in signal_keys:
                signal_prefix = "{}.{}".format(plot_prefix, signal_key)
                tasks = results_adjusted[adj_key][signal_key].keys()
                for task in tasks:
                    task_prefix = "{}.{}".format(signal_prefix, task)
                    # save out
                    out_file = "{}.genome.txt.gz".format(
                        task_prefix)
                    results_adjusted[adj_key][signal_key][task].to_csv(
                        out_file,
                        sep="\t", header=True, index=False, compression="gzip")
                    # plot
                    plot_cmd = "Rscript {}/plot.spacing.signals.indiv.R {} {}".format(
                        SCRIPT_DIR, out_file, task_prefix)
                    print plot_cmd
                    
                    os.system(plot_cmd)
                    
            # and take region sets and get gene set enrichments
            # TODO - adjust this, maybe even make wrapper to make consistent with multiplicity?
            syntax_enrichments = _get_gene_set_enrichment(
                results_adjusted[adj_key]["scores"].index.values,
                links_file,
                tss_file,
                background_gene_set_file,
                plot_prefix,
                OUT_DIR)
            syntax_enrichments["syntax"] = adj_key
            
            # add to summary
            if enrichments_summary is None:
                enrichments_summary = syntax_enrichments.copy()
            else:
                enrichments_summary = pd.concat([enrichments_summary, syntax_enrichments], axis=0)
                enrichments_summary = enrichments_summary.sort_values("term.name")
                
        # clean up summary
        enrichments_summary["log10pval"] = -np.log10(enrichments_summary["p.value"].values)
        enrichments_summary = enrichments_summary.drop("p.value", axis=1)
        summary_file = "{}/{}.summary.txt.gz".format(tmp_dir, pwm_name_clean)
        enrichments_summary.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")
                
        # remove substrings
        for substring in REMOVE_SUBSTRINGS:
            keep = [False if substring in term else True
                    for term in enrichments_summary["term.name"]]
            keep_indices = np.where(keep)[0]
            enrichments_summary = enrichments_summary.iloc[keep_indices]

        # remove exact strings
        keep = [False if term in REMOVE_EXACT_STRINGS else True
                for term in enrichments_summary["term.name"]]
        keep_indices = np.where(keep)[0]
        enrichments_summary = enrichments_summary.iloc[keep_indices]
        enrichments_summary = enrichments_summary.sort_values(["term.name", "log10pval"])
        summary_file = "{}/{}.summary.filt.txt.gz".format(OUT_DIR, pwm_name_clean)
        enrichments_summary.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")

        # TODO - here we are able to consider multiplicity too
        # NOTE though that here we don't capture solo with 1 motif regions
        # so then need to do this in another script i think

    return sig_pwms_file



def runall(args, prefix):
    """workflows for nn postprocessing
    """
    # set up logging, files, folders
    logger = logging.getLogger(__name__)
    logger.info("WORKFLOW: run inference on NN models")

    # set up data and results dir
    data_dir = args.outputs["data"]["dir"]
    run_shell_cmd("mkdir -p {}".format(data_dir))
    out_data = args.outputs["data"]

    results_dirname = "inference"
    results_dir = "{}/{}".format(args.outputs["results"]["dir"], results_dirname)
    args.outputs["results"][results_dirname] = {"dir": results_dir}
    run_shell_cmd("mkdir -p {}".format(results_dir))
    out_results = args.outputs["results"][results_dirname]
    
    # -------------------------------------------
    # ANALYSIS - homotypic - look at multiplicity
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------
    sig_pwms = "{}/{}".format(
        args.inputs["inference"][args.cluster]["sig_pwms_dir"],
        args.inputs["inference"][args.cluster]["sig_pwms.rna_filt.corr_filt"])

    # multiplicity
    sig_pwms = run_multiplicity_workflow(args, prefix, sig_pwms)

    # orientation, spacing
    sig_pwms = run_syntax_workflow(args, prefix, sig_pwms)
    
    # -------------------------------------------
    # ANALYSIS - homotypic - produce motif to gene set to gene function plots
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------


    

    

    return args


