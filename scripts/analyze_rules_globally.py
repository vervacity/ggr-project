
import os
import sys

import pandas as pd
import networkx as nx

from ggr.analyses.linking import regions_to_genes

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



def main():
    """analyses across all rules and effects on landscape
    """
    # input files
    mpra_file = sys.argv[1]
    rule_dir = sys.argv[2]
    
    # use the validated rules only
    rules = pd.read_csv(mpra_file, sep="\t")
    rules = rules[~rules["interaction"].str.contains("FAILED")]

    # get bed files
    bed_files = []
    for rule_idx in range(rules.shape[0]):
        rule_name = rules["grammar"].iloc[rule_idx]
        pwm1_name = rules["pwm1_clean"].iloc[rule_idx]
        pwm2_name = rules["pwm2_clean"].iloc[rule_idx]
        
        # get gml file
        gml_file = "{}/{}.gml".format(rule_dir, rule_name)
        graph = nx.read_gml(gml_file)

        # get bed file
        bed_file = "{}.{}_x_{}.bed.gz".format(rule_name, pwm1_name, pwm2_name)
        if not os.path.isfile(bed_file):
            print bed_file
            get_bed_from_nx_graph(
                graph,
                bed_file,
                interval_key="region",
                merge=True,
                return_key="region")
        bed_files.append(bed_file)
        
    # summary info
    all_regions_file = "rules_all_regions.bed.gz"
    if not os.path.isfile(all_regions_file):
        merge_cmd = "zcat {} | sort -k1,1 -k2,2n | bedtools merge -i stdin | gzip -c > {}".format(
            " ".join(bed_files), all_regions_file)
        print merge_cmd
        os.system(merge_cmd)
        

    # how many genes covered
    links_file = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/results/linking/proximity/ggr.linking.ALL.overlap.interactions.txt.gz"
    tss_file = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/annotations/ggr.rna.tss.expressed.bed.gz"
    out_file = "./test.txt.gz"

    # make a tss file of just dynamic genes
    dynamic_genes_mat = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/data/ggr.rna.counts.pc.expressed.timeseries_adj.pooled.rlog.dynamic.traj.mat.txt.gz"
    dynamic_genes = pd.read_csv(dynamic_genes_mat, sep="\t", index_col=0)
    #print dynamic_genes
    #quit()

    
    regions_to_genes(
        all_regions_file,
        links_file,
        tss_file,
        out_file,
        filter_genes=dynamic_genes.index.values,
        filter_by_score=0.5)
    

    
    return

main()
