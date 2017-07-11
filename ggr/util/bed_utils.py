"""Contains useful functions for working with BED files
"""

import os

def merge_regions(inputs_bed, out_bed):
    """Given a list of bed files, make into master file
    """
    merge_all = ("zcat {0} | "
                 "sort -k1,1 -k2,2n | "
                 "bedtools merge -i stdin | "
                 "gzip -c "
                 "> {1}").format(' '.join(inputs_bed), out_bed)
    if not os.path.isfile(out_bed):
        os.system(merge_all)

    return None


