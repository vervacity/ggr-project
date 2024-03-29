

import os
import re
import sys
import glob
import json

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
            fp.write("{0}\t{1}\t{2}\n".format(chrom, start, stop))

    if merge:
        # merge
        tmp_bed_file = "{}.tmp.bed".format(bed_file.split(".bed")[0])
        os.system("mv {} {}".format(bed_file, tmp_bed_file))
        os.system("cat {} | sort -k1,1 -k2,2n | bedtools merge -i stdin | bgzip > {}".format(
            tmp_bed_file, bed_file))
        os.system("rm {}".format(tmp_bed_file))

        # renumber
        tmp_bed_file = "{}.tmp.bed.gz".format(bed_file.split(".bed")[0])
        os.system("mv {} {}".format(bed_file, tmp_bed_file))
        #os.system("zcat {} | awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$3\"\t\"$1\":\"$2\"-\"$3\"\t\"NR\"\t+\" }}' | bgzip -c > {}".format(
        #    tmp_bed_file, bed_file))
        os.system("zcat {} | awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$3\"\ttest\"NR\"\t\"NR\"\t.\" }}' | bgzip > {}".format(
            tmp_bed_file, bed_file))
        os.system("rm {}".format(tmp_bed_file))

        # bgzip
        #new_tmp_bed_file = "{}.bed".format(tmp_bed_file.split(".tmp")[0])
        #os.system("mv {} {}".format(tmp_bed_file, new_tmp_bed_file))
        #os.system("bgzip {}".format(new_tmp_bed_file))
        #os.system("rm {}".format(tmp_bed_file))
    
    return return_intervals


def _make_json_bed_entry(bed_file, display_name, dir_url):
    """add a bed file to json
    """
    entry = {}
    entry["type"] = "bed"
    entry["url"] = "{}/{}".format(dir_url, bed_file)
    entry["name"] = display_name
    entry["mode"] = "full"
    
    return entry


def main():
    """pull the region sets as bed files to visualize results
    """
    # input files
    mpra_file = sys.argv[1]
    rule_dir = sys.argv[2]
    dir_url = "http://mitra.stanford.edu/kundaje/dskim89/ggr/paper"
    
    # for rule in rule dir, pull the BED file from it
    rules = pd.read_csv(mpra_file, sep="\t")
    rules = rules[~rules["interaction"].str.contains("FAILED")]
    json_entries = []
    for rule_idx in range(rules.shape[0]):
        rule_name = rules["grammar"].iloc[rule_idx]
        pwm1_name = rules["pwm1_clean"].iloc[rule_idx]
        pwm2_name = rules["pwm2_clean"].iloc[rule_idx]

        # get gml file
        gml_file = "{}/{}.gml".format(rule_dir, rule_name)
        graph = nx.read_gml(gml_file)

        # make bed file
        bed_file = "{}.{}_x_{}.bed.gz".format(rule_name, pwm1_name, pwm2_name)
        if not os.path.isfile(bed_file):
            print bed_file
            get_bed_from_nx_graph(
                graph,
                bed_file,
                interval_key="active",
                merge=True,
                return_key="region")
            os.system("chmod a+x {}".format(bed_file))
        display_name = bed_file.split(".")[-3]
        json_entry = _make_json_bed_entry(bed_file, display_name, dir_url)
        json_entries.append(json_entry)
            
        # also make a tbi file
        make_tbi = "tabix -p bed {}".format(bed_file)
        tbi_file = "{}.tbi".format(bed_file)
        if not os.path.isfile(tbi_file):
            os.system(make_tbi)
            
    # set up json file with appropriate names
    json_file = "combinatorial_rules.json"
    json_data = json.dumps(json_entries, indent=1)
    json_data = re.sub(r'"(.*?)"(?=:)', r'\1', json_data)
    with open(json_file, "w") as fp:
        fp.write(json_data)
        
    return

main()
