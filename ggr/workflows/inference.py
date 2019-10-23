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
from ggr.util.utils import run_shell_cmd
from ggr.util.go_utils import GoParams

from tronn.util.utils import DataKeys
from tronn.util.formats import array_to_bed
from tronn.interpretation.syntax import analyze_multiplicity


def run_multiplicity_workflow(args, prefix, sig_pwms_file):
    """for sig pwms, look at multiplicity
    """
    _MAX_COUNT = 6
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

    # results arrays
    all_results = {}
    signal_keys = args.inputs["inference"][args.cluster]["signal_keys"]
    for key in signal_keys:
        all_results[key] = []

    # TODO
    # read in region signals and rna signals
    # also read in dynamic genes
        
    # analyze each pwm
    keep_pwm_names = []
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

        if results is not None:
            keep_pwm_names.append(pwm_name_clean)
        
            # save to summary array
            for key in signal_keys:
                all_results[key].append(results[key]["count"])
        
            # per level, get region set and do gene set enrichments
            hits_per_region = results["hits_per_region"]
            for count_thresh in range(1, _MAX_COUNT+1):
                thresholded = hits_per_region[hits_per_region["hits"] >= count_thresh]
                thresholded_metadata = thresholded.index.values
                tmp_bed_file = "{}/{}.regions.count_thresh-{}.bed.gz".format(
                    TMP_DIR, pwm_name_clean, count_thresh)
                array_to_bed(thresholded_metadata, tmp_bed_file, merge=False)
                import ipdb
                ipdb.set_trace()
                
                regions_to_genes_w_correlation_filtering(
                    tmp_bed_file,
                    args.results["linking"]["links.proximity"],
                    args.annotations["tss.pc.bed"],
                    out_file,
                    region_signal,
                    rna_signal,
                    filter_by_score=0.5,
                    filter_genes=None,
                    corr_thresh=0,
                    pval_thresh=1) # TODO fix this

                # also run enrichment calculation
                
                
        
    # save out summary array
    h5_results_file = "{}/genome.homotypic.multiplicity_v_activity.h5".format(
        args.outputs["results"]["inference"]["dir"])
    if not os.path.isfile(h5_results_file):
        for key in signal_keys:
            summary_results = np.stack(all_results[key], axis=0)
            with h5py.File(h5_results_file, "a") as out:
                out_key = "{}/counts".format(key)
                out.create_dataset(out_key, data=summary_results)
                out[out_key].attrs["pwm_names"] = keep_pwm_names
    
    # and plot
    plot_cmd = "{}/plot.results.multiplicity.summary.R {} {}/genome.homotypic.multiplicity FALSE {}".format(
        "~/git/ggr-project/figs/fig_4.homotypic", h5_results_file,
        args.outputs["results"]["inference"]["dir"],
        " ".join(signal_keys))
    
    return






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
    # ANALYSIS - homotypic - look at multiplicity/spacing/orientation
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------
    sig_pwms = "{}/{}".format(
        args.inputs["inference"][args.cluster]["sig_pwms_dir"],
        args.inputs["inference"][args.cluster]["sig_pwms.rna_filt.corr_filt"])
    run_multiplicity_workflow(args, prefix, sig_pwms)

    # -------------------------------------------
    # ANALYSIS - homotypic - produce motif to gene set to gene function plots
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------


    

    

    return args


