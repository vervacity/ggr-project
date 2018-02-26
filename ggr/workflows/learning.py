"""Wrapper for utilizing tronn on datasets
"""

import os
import glob
import logging

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
    
        label_files = []
        for label_dir in args.outputs["results"]["label_dirs"]:

            logger.info("label directory: {}".format(label_dir))
            label_dir_files = sorted(glob.glob("{}/*.bed.gz".format(label_dir)))
            label_dir_files += sorted(glob.glob("{}/*.narrowPeak.gz".format(label_dir)))
            label_dir_files = [filename for filename in label_dir_files if "background" not in filename]
            logger.info("total files: {}".format(len(label_dir_files)))
            label_files += label_dir_files

        label_files_string = " ".join(label_files)
        logger.info("total label sets: {}".format(len(label_files)))

        # run tronn preprocess
        preprocess = (
            "tronn preprocess "
            "--labels {0} "
            "--parallel 24 "
            "-o {1}/datasets/{2} --prefix {2}").format(
                label_files_string,
                results_dir,
                prefix)
        run_shell_cmd("echo {} > {}/preprocess.tmp".format(preprocess, results_dir))
        run_shell_cmd("source preprocess.tmp")

    # -------------------------------------------
    # ANALYSIS 1 - train tronn model
    # input: tronn datasets
    # output: model checkpoint
    # -------------------------------------------
    train_dir = "{}/models/{}".format(results_dir, prefix)
    if not os.path.isdir(train_dir):
        # train
        train = (
            "tronn train "
            "--data_dir {0}/h5 "
            "--cvfold {1} "
            "--model basset "
            "--transfer_model_checkpoint {2} "
            "-o {3}/{4} --prefix {4}").format(
                dataset_dir,
                0,
                "",
                train_dir,
                prefix)
                    
    return args
