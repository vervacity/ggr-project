#!/usr/bin/env python

import sys
import h5py

import numpy as np

from tronn.plot.visualization import plot_pwm
from tronn.util.pwms import PWM



def main():
    """plot pwm weights
    """

    # args
    tfmodisco_file = sys.argv[1]
    metacluster_name = sys.argv[2]
    pattern_name = sys.argv[3]
    plot_prefix = sys.argv[4]
    background_freq = 0.25 # keep it simple
    
    # read in tfmodisco file, plot both fwd and reverse
    # any trimming??
    with h5py.File(tfmodisco_file, "r") as hf:
        # get probs from file
        probs = hf["metacluster_idx_to_submetacluster_results"][
            metacluster_name]["seqlets_to_patterns_result"]["patterns"][pattern_name]["sequence"]["fwd"][:]
        probs[probs == 0] = 0.0001

        # convert to weights
        weights = []
        for idx in range(probs.shape[0]):
            weights.append(
                np.log2(np.array(probs[idx,:]) / background_freq).tolist()
            )
        weights = np.array(weights).transpose(1,0)

        # load to PWM class, trim
        pwm = PWM(weights)
        pwm = pwm.chomp(ic_thresh=0.4)

        # plot fwd
        plot_file = "{}.{}.{}.fwd.pdf".format(plot_prefix, metacluster_name, pattern_name)
        plot_pwm(pwm.get_probs().transpose(), plot_file)

        # plot rev
        plot_file = "{}.{}.{}.rev.pdf".format(plot_prefix, metacluster_name, pattern_name)
        plot_pwm(pwm.reverse_complement().get_probs().transpose(), plot_file)
    
    return

main()
