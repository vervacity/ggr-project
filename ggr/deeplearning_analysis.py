"""Wrapper for utilizing tronn on datasets
"""

import glob


def run(args):
    """Run TRONN functions
    """

    # preprocess data
    label_files = []
    label_files_string = ""
    for label_dir_handle in args.params["label_dir_handles"]:

        label_dir_files = sorted(glob.glob("{}/*.bed.gz".format(args.folders[label_dir_handle])))
        label_dir_files += sorted(glob.glob("{}/*.narrowPeak.gz".format(args.folders[label_dir_handle])))
        label_files += label_dir_files

        if "timepoints" in label_dir_handle:
            label_files_string += "{0}/*narrowPeak.gz ".format(args.folders[label_dir_handle])
        else:
            label_files_string += "{0}/*bed.gz ".format(args.folders[label_dir_handle])
                        
        
    print "Found {} label files".format(len(label_files))
    print "for tronn: {}".format(label_files_string)

    return args
