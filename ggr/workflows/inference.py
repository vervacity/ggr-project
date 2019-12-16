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
from tronn.interpretation.syntax import analyze_syntax
from tronn.interpretation.syntax import recombine_syntax_results


def _setup_motifs_files(args):
    """convenience fn, make sure setup is same across 
    multiplicity/orientation/spacing workflows
    """
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

    return motifs_files


def _setup_signal_mats(args):
    """convenience fn, load signals
    """
    # read in rna signals
    rna_signal_file = args.outputs["data"][
        "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"]
    rna_signal_mat = pd.read_csv(rna_signal_file, sep="\t", index_col=0)
    rna_signal_mat[:] = np.subtract(
        rna_signal_mat.values,
        np.expand_dims(rna_signal_mat.values[:,0], axis=-1))

    # read in region signals
    atac_signal_file = args.outputs["data"][
        "atac.counts.pooled.rlog.dynamic.traj.mat"]
    region_signal_mat = pd.read_csv(atac_signal_file, sep="\t", index_col=0)
    region_signal_mat[:] = np.subtract(
        region_signal_mat.values,
        np.expand_dims(region_signal_mat.values[:,0], axis=-1))
    region_signal_mat = region_signal_mat.drop("d05", axis=1)

    return rna_signal_mat, region_signal_mat


def run_multiplicity_workflow(
        args, prefix, sig_pwms_file, out_dir,
        solo_filter=False, enrichments=True, plot=True):
    """for sig pwms, look at multiplicity
    """
    # params/dirs
    _MAX_COUNT = 5
    _MIN_HIT_REGIONS = 50
    OUT_DIR = "{}/{}".format(
        args.outputs["results"]["inference"]["dir"], out_dir)
    os.system("mkdir -p {}".format(OUT_DIR))
    TMP_DIR = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(TMP_DIR))
            
    # set up motifs files being used
    motifs_files = _setup_motifs_files(args)
    
    # read in the sig pwms
    sig_pwms_trajs = pd.read_csv(sig_pwms_file, sep="\t", header=0, index_col=0)
    sig_pwms = list(sig_pwms_trajs.index.values)
    logging.info("{} sig pwms found".format(len(sig_pwms)))
    
    # read in the list of all pwm names
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    with h5py.File(motifs_files["early"], "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    num_pwms = len(all_pwms)
    logging.info("{} total pwms".format(num_pwms))

    # read in dynamic genes
    filter_genes = pd.read_table(
        args.outputs["data"][
            "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"],
        index_col=0).index.values.tolist()

    # read in signals, other files/params
    rna_signal_mat, region_signal_mat = _setup_signal_mats(args)
    background_rna_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    links_file = args.outputs["results"]["linking"]["links.proximity"]
    tss_file = args.outputs["annotations"]["tss.pc.bed"]
    traj_to_group = args.inputs["inference"]["traj_to_group"]
    signal_keys = args.inputs["inference"][args.cluster]["signal_keys"]
    
    # results arrays
    all_results = {}
    for key in signal_keys:
        all_results[key] = []
    all_results["filt"] = []
    all_results["num_regions"] = []
        
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
            sig_motifs_files, pwm_indices, max_count=_MAX_COUNT, solo_filter=solo_filter)

        filt = 0
        if results is not None:
            keep_pwm_names.append(pwm_name_clean)
        
            # save to summary array
            for key in signal_keys:
                all_results[key].append(results[key]["count"])
            all_results["num_regions"] = results["num_regions_per_count"]
                
            # continue if not doing enrichments
            if not enrichments:
                filt = 1
                max_enriched_thresh = _MAX_COUNT
            else:
                
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
                    if not os.path.isfile(tmp_genes_file):
                        continue
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
        num_regions_per_count = np.stack(all_results["num_regions_per_count"], axis=0)
        with h5py.File(h5_results_file, "a") as out:
            out_key = "{}/num_regions".format(key)
            out.create_dataset(out_key, data=num_regions_per_count)
            out[out_key].attrs["pwm_names"] = keep_pwm_names

    # save out new sig pwms file
    new_sig_pwms_file = "{}/{}.multiplicity_filt.txt.gz".format(
        OUT_DIR,
        os.path.basename(sig_pwms_file).split(".txt")[0])
    new_sig_pwms = sig_pwms_trajs.loc[keep_sig_pwm_names]
    print new_sig_pwms
    new_sig_pwms.to_csv(new_sig_pwms_file, sep="\t", compression="gzip")
    
    # and plot
    # TODO move this plot outside to make easier to adjust?
    if plot:
        plot_cmd = "{}/plot.results.multiplicity.summary.R {} {}/genome.homotypic.multiplicity TRUE {}".format(
            "~/git/ggr-project/figs/fig_4.homotypic", h5_results_file,
            args.outputs["results"]["inference"]["dir"],
            " ".join(signal_keys))
        print plot_cmd
        os.system(plot_cmd)
        
    return new_sig_pwms_file


def run_syntax_workflow(args, prefix, sig_pwms_file, out_dir):
    """for sig pwms, look at multiplicity
    """
    # params
    _MIN_HIT_REGIONS = 50

    # dirs
    OUT_DIR = "{}/{}".format(args.outputs["results"]["inference"]["dir"], out_dir)
    os.system("mkdir -p {}".format(OUT_DIR))
    TMP_DIR = "{}/tmp".format(OUT_DIR)
    os.system("mkdir -p {}".format(TMP_DIR))
                        
    # set up motifs files being used
    motifs_files = _setup_motifs_files(args)
    
    # read in the sig pwms
    sig_pwms_trajs = pd.read_csv(sig_pwms_file, sep="\t", header=0, index_col=0)
    sig_pwms = list(sig_pwms_trajs.index.values)
    logging.info("{} sig pwms found".format(len(sig_pwms)))
    
    # read in the list of all pwm names
    max_val_key = DataKeys.WEIGHTED_PWM_SCORES_POSITION_MAX_VAL
    with h5py.File(motifs_files["early"], "r") as hf:
        all_pwms = hf[max_val_key].attrs["pwm_names"]
    num_pwms = len(all_pwms)
    logging.info("{} total pwms".format(num_pwms))
    
    # read in dynamic genes
    filter_genes = pd.read_table(
        args.outputs["data"][
            "rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat"],
        index_col=0).index.values.tolist()

    # read in signals, other files/params
    rna_signal_mat, region_signal_mat = _setup_signal_mats(args)
    background_rna_file = args.outputs["data"]["rna.counts.pc.expressed.mat"]
    links_file = args.outputs["results"]["linking"]["links.proximity"]
    tss_file = args.outputs["annotations"]["tss.pc.bed"]
    traj_to_group = args.inputs["inference"]["traj_to_group"]
    signal_keys = args.inputs["inference"][args.cluster]["signal_keys"]
    
    # results arrays
    all_results = {}
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

        results = {}
        
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
                    sig_motifs_files, pwm_indices, solo_filter=True)
                
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

            # adjust for correct patterns
            results_adjusted[adj_key] = recombine_syntax_results(
                results,
                adjustments[adj_key][0],
                adjustments[adj_key][1],
                signal_keys)

            # don't continue if does not exist or does not have enough regions
            if len(results_adjusted[adj_key].keys()) == 0:
                continue
            if results_adjusted[adj_key]["scores"].shape[0] < _MIN_HIT_REGIONS:
                continue

            # TODO consider spacing here somehow?
            
            # debug
            print adj_key, results_adjusted[adj_key]["scores"].shape[0]

            # convert to BED
            tmp_bed_file = "{}/{}.regions.oriented.{}.bed.gz".format(
                TMP_DIR, pwm_name_clean, adj_key)
            if not os.path.isfile(tmp_bed_file):
                array_to_bed(
                    results_adjusted[adj_key]["scores"].index.values,
                    tmp_bed_file, merge=False)

            # get linked genes
            tmp_genes_file = "{}/{}.linked_genes.oriented.{}.txt.gz".format(
                TMP_DIR, pwm_name_clean, adj_key)
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
            enrichment_file = "{}/{}.linked_genes.oriented.{}.go_gprofiler.txt".format(
                TMP_DIR, pwm_name_clean, adj_key)
            if not os.path.isfile(enrichment_file):
                run_gprofiler(
                    tmp_genes_file,
                    background_rna_file,
                    TMP_DIR,
                    ordered=True)

            # check if any enrichment, if not then continue
            if not is_enriched(enrichment_file):
                continue

            # read in file and clean
            syntax_enrichments = pd.read_csv(enrichment_file, sep="\t")
            syntax_enrichments = syntax_enrichments[syntax_enrichments["domain"] == "BP"]
            syntax_enrichments = syntax_enrichments[
                ["term.id", "p.value", "term.name"]]
            syntax_enrichments["syntax"] = adj_key
            print "term count:", syntax_enrichments.shape[0]
            
            # add to summary
            if enrichments_summary is None:
                enrichments_summary = syntax_enrichments.copy()
            else:
                enrichments_summary = pd.concat([enrichments_summary, syntax_enrichments], axis=0)
                enrichments_summary = enrichments_summary.sort_values("term.name")
            
            # if passed all of this, move results to not tmp
            os.system("cp {} {}".format(
                enrichment_file, OUT_DIR))
            
            # and now plot: plot prefix
            plot_prefix = "{}/{}.{}".format(
                TMP_DIR, pwm_name_clean, adj_key)

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
            plot_cmd = "plot.homotypic.spacing.freq.indiv.R {} {} {}".format(
                score_spacing_distr_file, num_examples, plot_prefix)
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
                    plot_cmd = "plot.homotypic.spacing.signals.indiv.R {} {}".format(
                        out_file, task_prefix)
                    print plot_cmd
                    
                    os.system(plot_cmd)
                    
        # continue if nothing enriched
        if enrichments_summary is None:
            continue
                
        # clean up summary
        enrichments_summary["log10pval"] = -np.log10(enrichments_summary["p.value"].values)
        enrichments_summary = enrichments_summary.drop("p.value", axis=1)
        summary_file = "{}/{}.summary.txt.gz".format(OUT_DIR, pwm_name_clean)
        enrichments_summary.to_csv(summary_file, sep="\t", index=False, header=True, compression="gzip")

        # TODO adjust cleanup
        if False:
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

    return None



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
    # NN ANALYSIS - scanmotifs, get differential
    # input: dynamic traj data
    # output: differential motifs, scanned examples
    # -------------------------------------------

    # TRONN: scanmotifs
    # TRONN: call_differential_motifs
    # TRONN: intersect_pwm_x_rna
    
    # sig pwms file
    sig_pwms = "{}/{}".format(
        args.inputs["inference"][args.cluster]["sig_pwms_dir"],
        args.inputs["inference"][args.cluster]["sig_pwms.rna_filt.corr_filt"])
    
    # -------------------------------------------
    # ANALYSIS - homotypic - look at multiplicity
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------

    # observing multiplicity (homotypic clusters, no higher syntax) in the genome

    # 1) across dynamic instances, do we see homotypic clusters of motifs?
    # do not use solo filter, just care about if we see them dispersed in dynamic
    # regulatory regions. NOTE that this will give you GO enrichments
    # ANSWER: yes.
    # output: heatmap of counts of each instance
    analysis_dir = "homotypic.dynamic.all.multiplicity.TMP"
    if not os.path.isdir(analysis_dir):
        # TODO focus on results file (to split out plotting)?
        run_multiplicity_workflow(
            args, prefix, sig_pwms, analysis_dir,
            solo_filter=False, enrichments=True, plot=False)

    quit()
        
    # 2) can they drive accessibility alone?
    # ANSWER yes, in some cases
    analysis_dir = "homotypic.dynamic.single_motif_attribution.multiplicity"
    #if not os.path.isdir(multiplicity_dir):
    #    run_multiplicity_workflow(args, prefix, sig_pwms, analysis_dir, solo_filter=True, plot=False, enrichments=False)

    # 3) confirm with simulations
    # NN ANALYSIS - run multiplicity simulations
    # determine thresholds for activation AND which ones activate accessibility alone    
    # NOTE that simulations provide the threshold (manually collect) and ability to drive accessibility alone
    
    # -------------------------------------------
    # ANALYSIS - homotypic - look at orientation/spacing
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------

    # observing orientation/spacing in the genome

    # 1) across dynamic instances, do we see orientation/spacing constraints (leading
    # to accessibility/modification differences) on motifs? do they drive accessibility alone?
    # use solo filter, this is the only way to remove effect of other TFs
    # use enrichment, which then tells us if this orientation/spacing matters
    # ANSWER: no to orientation, soft spacing constraint. yes they drive accessibility alone in some cases
    analysis_dir = "homotypic.dynamic.single_motif_attribution.orientation_spacing"
    if not os.path.isdir(orientation_spacing_dir):
        run_syntax_workflow(args, prefix, sig_pwms, analysis_dir, solo_filter=True)

    # 2) confirm with simulations
    # NN ANALYSIS - run orientation/spacing simulations
    # use to determin spacing constraints AND which ones activate accessibility alone

    # 3) rerun with all regions (ex GRHL)
    
    
    # TODO build a manual summary: tells an aggregator fn which ones to grab
    
    
        
    # -------------------------------------------
    # ANALYSIS - homotypic - produce motif to gene set to gene function plots
    # input: scanmotifs files
    # output: gene function chart
    # -------------------------------------------
        
    # first, determine if multiplicity/orientation/spacing are enriched in the genome for each motif
    # given above, we conclude that multiplicity is a thing, orientation not, spacing yes

    # also use sims to give us predicted threshs and which ones activate on their own
    # (hence can drive downstream function)
    
    # build a manual summary file that will grab the appropriate gene function chart

    
    # -------------------------------------------
    # ANALYSIS - heterotypic - buildgrammars
    # input: mutatemotifs outputs
    # output: grammars
    # -------------------------------------------

    # tronn: run mutatemotifs using sig pwms above
    # tronn: build grammars
    # build 3s, but only consider 2s for main figure?

    # still need a way to filter 3s to fit in supplements
    # maybe higher threshold for 3s?
    # or need a different test for 3s to get enrichment
    
    
    # filter grammars
    

    

    return args


