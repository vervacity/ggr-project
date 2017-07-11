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

def id_to_bed(id_file, bed_file):
    """Convert id of form chr:start-stop into BED
    """
    if os.path.splitext(idfile)[1] == '.gz':
        pipe_in = 'zcat'
    else:
        pipe_in = 'cat'
        
    convert = ("{0} {1} | ",
               "awk -F ':' '{{ print $1\"\t\"$2 }}' | ",
               "awk -F '-' '{{ print $1\"\t\"$2 }}' | ",
               "awk -F '(' '{{ print $1 }}' | ",
               "gzip -c > {1}").format(pipe_in, id_file, bed_file)
    print convert
    os.system(convert)

    return None
