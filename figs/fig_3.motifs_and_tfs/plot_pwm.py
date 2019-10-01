#!/usr/bin/env python

# adjust matplotlib params
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['axes.linewidth'] = 0.1
matplotlib.rcParams['xtick.labelsize'] = 4
matplotlib.rcParams['xtick.major.width'] = 0.1
matplotlib.rcParams['xtick.major.size'] = 1
matplotlib.rcParams['ytick.labelsize'] = 4
matplotlib.rcParams['ytick.major.width'] = 0.1
matplotlib.rcParams['ytick.major.size'] = 1


import sys

import matplotlib.pyplot as plt

from tronn.util.pwms import MotifSetManager # TODO move this?
from tronn.plot.visualization import plot_pwm

def main():
    """
    """
    pwm_file = sys.argv[1]
    pwm_name = sys.argv[2]
    
    # read in pwm file, find pwm of interest and plot out
    pwm_list = MotifSetManager.read_pwm_file(pwm_file)
    print_pwms = [pwm for pwm in pwm_list if pwm_name in pwm.name]

    # for each pwm, plot
    for pwm in print_pwms:
        plot_file = "{}.weights.pdf".format(pwm.name)
        plot_pwm(pwm.get_probs().transpose(), plot_file)

    
    return

main()
