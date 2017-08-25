"""Wrapper for utilizing tronn on datasets
"""

import glob


def run(args):
    """Run TRONN functions
    """

    # preprocess data
    label_files = []
    for label_dir_handle in args.params["label_dir_handles"]:

        label_dir_files = sorted(glob.glob("{}/*.bed.gz".format(args.folders[label_dir_handle])))
        label_files += label_dir_files

    print "Found {} label files".format(label_files)
        

    return args
