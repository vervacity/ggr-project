"""General code pieces to ease pipelining
"""

import os
import subprocess
import logging


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
