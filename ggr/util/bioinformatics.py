"""Contains wrappers to quickly work with bioinformatics tools
"""

import os

def run_homer(positives_bed, background_bed, out_dir, parallel=12):
    """Generic wrapper to run homer on a set of regions (vs background)
    """
    run_homer = ("findMotifsGenome.pl {0} hg19 {1}"
                 "-bg {2}"
                 "-p {3}"
                 "-nomotif").format(positives_bed,
                                    out_dir,
                                    background_bed,
                                    parallel)
    print run_homer
    os.system(run_homer)

    return None
