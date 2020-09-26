

import os
import sys

import pandas as pd

import networkx as nx


def get_bed_from_nx_graph(
        graph,
        bed_file,
        interval_key="active",
        merge=True,
        return_key="region"):
    """get BED file from nx examples
    """
    if isinstance(graph.graph["examples"], basestring):
        graph.graph["examples"] = graph.graph["examples"].split(",")
    examples = list(graph.graph["examples"])

    return_intervals = []
    with open(bed_file, "w") as fp:
        for region_metadata in examples:
            interval_types = region_metadata.split(";")
            interval_types = dict([
                interval_type.split("=")[0:2]
                for interval_type in interval_types])
            interval_string = interval_types[interval_key]
            return_intervals.append(interval_types[return_key])
            
            chrom = interval_string.split(":")[0]
            start = interval_string.split(":")[1].split("-")[0]
            stop = interval_string.split("-")[1]
            fp.write("{}\t{}\t{}\n".format(chrom, start, stop))

    if merge:
        tmp_bed_file = "{}.tmp.bed".format(bed_file.split(".bed")[0])
        os.system("mv {} {}".format(bed_file, tmp_bed_file))
        os.system("cat {} | sort -k1,1 -k2,2n | bedtools merge -i stdin | gzip -c > {}".format(
            tmp_bed_file, bed_file))
        os.system("rm {}".format(tmp_bed_file))
    
    return return_intervals


def main():
    """pull the region sets as bed files to visualize results
    """
    # input files
    mpra_file = sys.argv[1]
    rule_dir = sys.argv[2]

    # for rule in rule dir, pull the BED file from it
    rules = pd.read_csv(mpra_file, sep="\t")
    rules = rules[~rules["interaction"].str.contains("FAILED")]
    for rule_idx in range(rules.shape[0]):
        rule_name = rules["grammar"].iloc[rule_idx]
        pwm1_name = rules["pwm1_clean"].iloc[rule_idx]
        pwm2_name = rules["pwm2_clean"].iloc[rule_idx]

        # get gml file
        gml_file = "{}/{}.gml".format(rule_dir, rule_name)
        graph = nx.read_gml(gml_file)

        # make bed file
        bed_file = "{}.{}_x_{}.bed.gz".format(rule_name, pwm1_name, pwm2_name)
        print bed_file
        get_bed_from_nx_graph(
            graph,
            bed_file,
            interval_key="active",
            merge=True,
            return_key="region")
        
    # set up json file with appropriate names
    bed_files = sorted(glob.glob("./*bed.gz"))
    
    
    
    return

main()
