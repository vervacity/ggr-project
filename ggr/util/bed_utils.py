"""Contains useful functions for working with BED files
"""

import os


from ggr.util.utils import run_shell_cmd

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


def id_to_bed(id_file, bed_file, sort=False, remove_extra_columns=True):
    """Convert id of form chr:start-stop into BED
    
    Args:
      id_file: file with IDs as first column, 
        can include other data
      bed_file: output BED file
    """
    if os.path.splitext(id_file)[1] == '.gz':
        pipe_in = 'zcat'
    else:
        pipe_in = 'cat'

    if sort:
        sort_file = "sort -k1,1 -k2,2n | "
    else:
        sort_file = ""

    if remove_extra_columns:
        remove = "awk -F '\t' '{{ print $1 }}' | "
    else:
        remove = ""
        
    convert = (
        "{0} {1} | "
        "grep -v region | "
        "{2}"
        "awk -F ':' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '-' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '(' '{{ print $1 }}' | "
        "{3}"
        "gzip -c > {4}").format(
            pipe_in, id_file, remove, sort_file, bed_file)
    print convert
    os.system(convert)

    return None


def get_midpoint_and_extend(bed_file, chrom_sizes_file, extend_len, out_file):
    """Given a bed file and extend length, get midpoints
    of regions and extend outwards
    """
    slop_bed = (
        "zcat {0} | "
        "awk -F '\t' 'BEGIN{{OFS=\"\t\"}} "
        "{{ midpoint=$2+int(($3-$2)/2); "
        "$2=midpoint; $3=midpoint+1; print }}' | "
        "bedtools slop -i stdin -g {1} -b {2} | "
        "gzip -c > {3}").format(
            bed_file,
            chrom_sizes_file,
            extend_len,
            out_file)
    run_shell_cmd(slop_bed)

    return None
