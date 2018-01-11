"""General code pieces to ease pipelining
"""

import os
import subprocess
import logging

from ggr.util.parallelize import setup_multiprocessing_queue
from ggr.util.parallelize import run_in_parallel

def add_folder(args, handle, parent_folder, child_folder):
    """Makes folder_name inside parent_folder and stores in args
    """
    args.folders[handle] = "{}/{}".format(
        parent_folder, child_folder)
    os.system("mkdir -p {}".format(
        args.folders[handle]))
    
    return None


def run_shell_cmd(cmd): 
    """Set up shell command runs
    """
    logger = logging.getLogger(__name__)
    logger.debug(cmd)
    os.system(cmd)

    if False:

        try:
            p = subprocess.Popen(cmd, shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                preexec_fn=os.setsid)
            pid = p.pid
            pgid = os.getpgid(pid)
            ret = ''
            while True:
                line = p.stdout.readline()
                if line=='' and p.poll() is not None:
                    break
                # log.debug('PID={}: {}'.format(pid,line.strip('\n')))
                print('PID={}: {}'.format(pid,line.strip('\n')))
                ret += line
            p.communicate() # wait here
            if p.returncode > 0:
                raise subprocess.CalledProcessError(
                    p.returncode, cmd)
            return ret.strip('\n')
        except:
            # kill all child processes
            os.killpg(pgid, signal.SIGKILL)
            p.terminate()
            raise Exception('Unknown exception caught. PID={}'.format(pid))


def parallel_copy(in_files, out_dir, num_threads=24):
    """Parallel copy - use more threads
    """
    copy_queue = setup_multiprocessing_queue()

    for in_file in in_files:
        copy_cmd = "rsync -avz --progress {} {}/".format(
            in_file, out_dir.rstrip("/"))
        copy_queue.put([run_shell_cmd, [copy_cmd]])
    
    run_in_parallel(copy_queue, parallel=num_threads)
    
    return None
