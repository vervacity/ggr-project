#!/usr/bin/env python

import h5py

import networkx as nx
import numpy as np

from tronn.plot.visualization import plot_weights_group
from tronn.plot.visualization import scale_scores
from tronn.util.utils import DataKeys


def main():
    """help selecting a region
    """
    # check grammar file
    #grammar_file = "/srv/scratch/dskim89/ggr/inference.2019-03-12/dmim.shuffle/grammars.annotated.manual_filt.merged.final.viz/ggr.TRAJ_LABELS-2.grammar-60.annot-274.gml"
    #grammar_file = "/srv/scratch/dskim89/ggr/inference.2019-03-12/dmim.shuffle/grammars.annotated.manual_filt.merged.final.viz/ggr.TRAJ_LABELS-2.grammar-59.annot-272.gml"
    grammar_file = "/srv/scratch/dskim89/ggr/inference.2019-03-12/dmim.shuffle/grammars.annotated.manual_filt.merged.final.viz/ggr.TRAJ_LABELS-1.grammar-99.annot-167.gml"
    #grammar_file = "/srv/scratch/dskim89/ggr/inference.2019-03-12/dmim.shuffle/grammars.annotated.manual_filt.merged.final.viz/ggr.TRAJ_LABELS-2.grammar-99.annot-316.gml"

    
    # read in grammar file
    graph = nx.read_gml(grammar_file)
    examples = graph.graph["examples"].split(",")

    # filter for within locus
    #keep_examples = [example for example in examples if "chr1:15" in example]
    keep_examples = [example for example in examples if "chr12:5" in example]
    keep_examples = list(set(keep_examples))
    print keep_examples

    # get region id
    #region_id = "region=chr1:152538683-152539769;active=chr1:152539433-152539633;features=chr1:152539033-152540033" # traj 1
    region_id = "region=chr12:53092412-53093646;active=chr12:53092812-53093012;features=chr12:53092412-53093412" # traj 1
    
    # plot
    INFER_DIR = "/srv/scratch/dskim89/ggr/inference.2019-03-12/dmim.shuffle"
    TRAJ_DIR = "{}/TRAJ_LABELS-1".format(INFER_DIR)
    data_file = "{}/ggr.dmim.h5".format(TRAJ_DIR)

    LEFT_CLIP = 420
    RIGHT_CLIP = 580
    FINAL_LEFT_CLIP = 20
    FINAL_RIGHT_CLIP = 120

    # pull example idx based on region id
    with h5py.File(data_file, "r") as hf:
        metadata = hf[DataKeys.SEQ_METADATA][:,0]
    example_idx = np.where(np.isin(metadata, [region_id]))[0][0]
    print example_idx
    
    # pull data and clip to match
    with h5py.File(data_file, "r") as hf:
        data = hf[DataKeys.WEIGHTED_SEQ][example_idx][:,LEFT_CLIP:RIGHT_CLIP]
        sig = hf[DataKeys.WEIGHTED_SEQ_ACTIVE][example_idx]

    # normalize
    adjusted_scores = scale_scores(data, sig)

    # plot full active region
    out_file = "importances.pdf"
    plot_weights_group(adjusted_scores, out_file, sig_array=sig)

    # clip again and plot
    out_file = "importances.clipped.pdf"
    clipped_scores = adjusted_scores[:,FINAL_LEFT_CLIP:FINAL_RIGHT_CLIP]
    clipped_sig = sig[:,FINAL_LEFT_CLIP:FINAL_RIGHT_CLIP]
    plot_weights_group(clipped_scores, out_file, sig_array=clipped_sig)
    
    

    
    return

main()
