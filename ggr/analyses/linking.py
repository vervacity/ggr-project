"""code for enh-gene linking
"""

import os

def bed_to_gene_set_by_proximity(
        bed_file,
        tss_file,
        k_nearest=3,
        max_dist=250000):
    """assuming proximal linking captures the majority of linking
    """
    # bedtools closest
    tmp_file = "tss.overlap.tmp.txt"
    closest = "bedtools closest -d -k {} -a {} -b {} > {}".format(
        k_nearest, bed_file, tss_file, tmp_file)
    os.system(closest)

    # load results and use distance cutoff
    data = pd.read_table(tmp_file, header=None)
    data = data[data[9] < max_dist]
    print data
    #data = data[[6,9]]
    #data.columns = ["gene_id", "distance"]
    #data = data.set_index("gene_id")

    # save out?

    # also run gprofiler?
    
    # cleanup
    os.system("rm {}".format(tmp_file))
    
    return data
