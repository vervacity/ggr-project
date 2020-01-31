#!/usr/bin/env python

import sys

from tronn.util.pwms import MotifSetManager
from tronn.plot.visualization import plot_pwm

def main():
    """plot pwm weights
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
