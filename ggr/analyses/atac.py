"""support analyses for accessibility
"""

import os
import re

import numpy as np
import pandas as pd

from ggr.analyses.counting import make_count_matrix


def compare_to_scATAC_OLD(atac_counts, sc_atac_counts, out_dir):
    """compare results
    """
    # read in data
    bulk_atac = pd.read_csv(atac_counts, sep="\t", index_col=0)
    sc_atac = pd.read_csv(sc_atac_counts, sep="\t", index_col=0)
    days = ["d0", "d3", "d6"]
    
    # summarize bulk atac
    bulk_atac = bulk_atac[["d00", "d30", "d50"]]
    bulk_atac.columns = days
    
    # summarize sc ATAC
    days_sc = ["D0", "D3", "D6"]
    for day in days_sc:
        headers = [col for col in sc_atac.columns if day in col]
        sc_atac[day] = sc_atac[headers].sum(axis=1)
    sc_atac = sc_atac[days_sc]
    sc_atac.columns = days
    
    # get a region id mapping
    mapping_file = "{}/sc_to_bulk_atac.mapping.txt.gz".format(out_dir)
    if not os.path.isfile(mapping_file):
        # sc atac bed
        sc_atac_bed = "{}/sc_atac.tmp.bed.gz".format(out_dir)
        sc_atac_regions = sc_atac.reset_index()["Peak_position"].str.split("_", n=3, expand=True)
        sc_atac_regions["id"] = sc_atac.index
        sc_atac_regions["score"] = 100
        sc_atac_regions["strand"] = "."
        sc_atac_regions.to_csv(sc_atac_bed, sep="\t", header=False, index=False, compression="gzip")

        # bulk atac bed
        bulk_atac_bed = "{}/bulk_atac.tmp.bed.gz".format(out_dir)
        bulk_atac_regions = bulk_atac.reset_index()["index"].str.split(":", n=2, expand=True)
        bulk_atac_regions[["start", "stop"]] = bulk_atac_regions[1].str.split("-", n=2, expand=True)
        bulk_atac_regions = bulk_atac_regions.drop(1, axis=1)
        bulk_atac_regions["id"] = bulk_atac.index
        bulk_atac_regions["score"] = 100
        bulk_atac_regions["strand"] = "."
        bulk_atac_regions.to_csv(bulk_atac_bed, sep="\t", header=False, index=False, compression="gzip")

        # intersect - keep all A features (scATAC), overlap is 0 if none in B (bulk)
        mapping_cmd = "bedtools intersect -wao -a {} -b {} | gzip -c > {}".format(
            sc_atac_bed, bulk_atac_bed, mapping_file)
        os.system(mapping_cmd)
        
    # read in and clean up mapping
    mapping = pd.read_csv(mapping_file, sep="\t", header=None)
    overlap_thresh = 0.33
    mapping["sc_fract"] = mapping[12] / (mapping[2] - mapping[1])
    mapping["bulk_fract"] = mapping[12] / (mapping[8] - mapping[7])
    mapping["mapped"] = (
        (mapping["sc_fract"] > overlap_thresh) & (mapping["bulk_fract"] > overlap_thresh)).astype(int)
    mapping = mapping[[3, 9, "mapped"]]
    mapping.columns = ["sc_atac_ids", "bulk_atac_ids", "mapped"]
    
    # add in count data, keeping all ids
    sc_atac_merge = sc_atac.copy()
    sc_atac_merge.columns = ["sc.{}".format(day) for day in days]
    data = mapping.merge(sc_atac_merge, how="right", left_on="sc_atac_ids", right_index=True)
    bulk_atac_merge = bulk_atac.copy()
    bulk_atac_merge.columns = ["bulk.{}".format(day) for day in days]
    data = data.merge(bulk_atac_merge, how="right", left_on="bulk_atac_ids", right_index=True)

    # at this point, save out and send to R for plotting
    out_file = "{}/sc_vs_bulk.data.mat.txt.gz".format(out_dir)
    data.to_csv(out_file, sep="\t", header=True, index=False, compression="gzip")

    # plotting
    plot_cmd = "~/git/ggr-project/R/plot.sc_vs_bulk.R {} {}".format(out_file, out_dir)
    print plot_cmd
    os.system(plot_cmd)

    # TODO - get overlap numbers and plot out? - should really probably use an rlog thresh
    sc_days = ["sc.{}".format(day) for day in days]
    bulk_days = ["bulk.{}".format(day) for day in days]
    
    for sc_day in sc_days:
        for bulk_day in bulk_days:
            sc_ids = set(data[data[sc_day] > 0].index.values)
            bulk_ids = set(data[data[bulk_day] > 100].index.values)
            print sc_day, bulk_day, len(sc_ids.intersection(bulk_ids))

    
    quit()

    
    return


def compare_to_scATAC(sc_atac_counts, bulk_read_files, out_dir, threads=12):
    """compare results
    """
    # read in data
    sc_atac = pd.read_csv(sc_atac_counts, sep="\t", index_col=0)
    days = ["d0", "d3", "d6"]

    # make a bed file from sc atac data
    sc_bed_file = "{}/sc.regions.bed.gz".format(out_dir)
    sc_regions = sc_atac.reset_index()["Peak_position"].str.split("_", n=3, expand=True)
    sc_regions.to_csv(
        sc_bed_file, columns=[0, 1, 2],
        compression="gzip", sep="\t", header=False, index=False)
    sc_bed_sorted_file = "{}.sorted.bed.gz".format(sc_bed_file.split(".bed")[0])
    sort_cmd = "zcat {} | sort -k1,1 -k2,2n | gzip -c > {}".format(
        sc_bed_file, sc_bed_sorted_file)
    os.system(sort_cmd)
    
    # make a counts matrix
    bulk_counts_mat_file = "{}/bulk.counts_in_sc_regions.mat.txt.gz".format(out_dir)
    if not os.path.isfile(bulk_counts_mat_file):
        make_count_matrix(
            sc_bed_sorted_file,
            bulk_read_files,
            bulk_counts_mat_file,
            "ATAC",
            adjustment="ends",
            tmp_dir=out_dir,
            parallel=threads)

    # merge it into the sc atac data
    sc_v_bulk_mat_file = "{}/sc_and_bulk.counts.mat.txt.gz".format(out_dir)
    if not os.path.isfile(sc_v_bulk_mat_file):
        bulk_atac = pd.read_csv(bulk_counts_mat_file, sep="\t")
        # fix bulk region ids to merge correctly with sc atac
        bulk_atac.index = bulk_atac.index.str.replace(":", "_")
        bulk_atac.index = bulk_atac.index.str.replace("-", "_")
        # rename sc headers
        sc_headers = sc_atac.columns
        sc_headers = [re.sub("^.+scATAC-D", "sc.d", header) for header in sc_headers]
        sc_headers = [re.sub(".st.rmdup.flt.bam", "", header) for header in sc_headers]
        sc_headers = [re.sub("-AR-", "_b", header) for header in sc_headers]
        sc_atac.columns = sc_headers
        # merge
        all_atac = sc_atac.merge(bulk_atac, left_index=True, right_index=True)
        all_atac.to_csv(
            sc_v_bulk_mat_file, sep="\t", index=True, header=True, compression="gzip")

    # run rlog norm
    if False:
        sc_v_bulk_rlog_mat_file = "{}.rlog.mat.txt.gz".format(
            sc_v_bulk_mat_file.split(".mat")[0])
        rlog_cmd = "normalize.rlog.R {} {}".format(
            sc_v_bulk_mat_file, sc_v_bulk_rlog_mat_file)
        print rlog_cmd
        os.system(rlog_cmd)
    
    # send to R to plot
    r_cmd = "~/git/ggr-project/R/plot.sc_vs_bulk.dim_reduction.R {} {}".format(
        sc_v_bulk_mat_file, "{}/sc_vs_bulk".format(out_dir))
    print r_cmd
    #os.system(r_cmd)
    
    quit()
    
    return
