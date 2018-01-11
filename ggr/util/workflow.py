# Description: simple workflow code

import os
import logging


class Analysis(object):

    def __init__(self, function, in_files, out_files):
        self.function = function
        self.in_files = in_files
        self.out_files = out_files

        
    def should_run(self):
        """This checks what conditions should be met for task to run
        """

        return

    
    def run(self):
        """Perform the check, and run if needed
        """
        

        return




    

def run_parallel_step(function, in_files, out_files, args):
    """Run the function in n parallel steps?
    """
    
    return None


def run_step(
        function,
        in_files,
        out_files,
        step_args=[],
        step_kwargs={},
        description=""):
    """Simple workflow control tool: checks mod dates on files 
    and whether they exist or not to decide whether to run code
    or not.
    """
    # check if all input files exist and get latest mod date
    # if not, should stop because an input file was missing
    last_input_file_mtime = -float("inf")
    for in_file in in_files:
        assert os.path.isfile(in_file)
        mtime = os.path.getmtime(in_file)
        if mtime > last_input_file_mtime:
            last_input_file_mtime = mtime

    # check if all output files exist and get latest mod date
    first_output_file_mtime = float("inf")
    for out_file in out_files:
        if os.path.isfile(out_file):
            mtime = os.path.getmtime(in_file)
            if mtime < first_input_file_mtime:
                first_output_file_mtime = mtime

    # check mtimes, if output file was mod before input file, then re-run
    if last_input_file_mtime > first_output_file_mtime:
        function(*step_args, **step_kwargs)

    return None
