"""General code pieces to ease pipelining
"""

import os

def add_folder(args, handle, parent_folder, child_folder):
    """Makes folder_name inside parent_folder and stores in args
    """
    args.folders[handle] = "{}/{}".format(
        parent_folder, child_folder)
    os.system("mkdir -p {}".format(
        args.folders[handle]))
    
    return None
