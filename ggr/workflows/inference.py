"""description: workflow for inference downstream
"""

import os
import glob
import logging

from ggr.util.utils import run_shell_cmd

from tronn.interpretation.syntax import analyze_multiplicity


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
    

    
    print args.inputs["inference"]
    

    # -------------------------------------------
    # ANALYSIS - homotypic - produce motif to gene set to gene function plots
    # input: scanmotifs files
    # output: syntax results
    # -------------------------------------------


    

    

    return args


