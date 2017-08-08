"""Contains useful functions for working with BED files
"""

import os

def merge_regions(bed_files, out_bed):
    """Given a list of bed files, make into master file

    Args:
      bed_files: list of bed files
      out_bed: name of output merged bed file
    """
    merge_all = ("zcat {0} | "
                 "sort -k1,1 -k2,2n | "
                 "bedtools merge -i stdin | "
                 "gzip -c "
                 "> {1}").format(' '.join(bed_files), out_bed)
    print merge_all
    os.system(merge_all)

    return None


def id_to_bed(id_file, bed_file):
    """Convert id of form chr:start-stop into BED
    
    Args:
      id_file: file with IDs as first column, 
        can include other data
      bed_file: output BED file
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
