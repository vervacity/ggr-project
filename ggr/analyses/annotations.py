"""Contains methods for getting annotations
"""

import os
import gzip
import glob
import logging
import pandas as pd

from ggr.util.utils import run_shell_cmd


def get_proteincoding_gene_list_from_gtf(gtf_file, out_file):
    """Using a gtf file of protein coding genes, extract ids
    Note that it leaves the decimals on there

    Args:
      gtf_file: gtf annotation file with protein coding annotations
      out_file: output file of gene ids
    """
    assert gtf_file.endswith("gtf")
    assert out_file.endswith(".gz")
    extract_gene_ids = (
        "cat {0} | "
        "grep 'protein_coding' | "
        "awk -F '\"' '{{ print $2 }}' | "
        "awk -F '.' '{{ print $1 }}' | "
        "sort | "
        "uniq | "
        "gzip -c > {1}").format(
            gtf_file,
            out_file)
    run_shell_cmd(extract_gene_ids)
    
    return None


def get_proteincoding_tss_bed_from_gtf(gtf_file, out_file):
    """Make a TSS BED file from a gtf file

    Args:
      gtf_file: gtf annotation file
      out_file: BED file with TSS regions
    """
    assert gtf_file.endswith("gtf")
    assert out_file.endswith(".gz")
    make_tss_file = (
        "cat {0} | "
        "grep -P '\tgene\t' | "
        "grep 'protein_coding' | "
        "grep -v 'level 3' | "
        "awk -F '[\t|\"]' "
        "'{{ print $1\"\t\"$4\"\t\"$5\"\t\"$10\"\t0\t\"$7 }}' | "
        "awk -F '\t' 'BEGIN{{ OFS=\"\t\" }} "
        "{{ if ($6==\"+\") {{ $3=$2-1; $2=$2-2 }} "
        "else {{ $2=$3; $3=$3+1 }} print }}' | "
        
        "awk -F '[\t|\.]' 'BEGIN{{ OFS=\"\t\" }} "
        "{{ print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$6\"\t\"$7 }}' | "
        
        "sort -k1,1 -k2,2n | "
        "gzip -c "
        "> {1}").format(
            gtf_file,
            out_file)
    run_shell_cmd(make_tss_file)
    
    return None



def get_ensembl_to_geneid_mapping(ensembl_id_file, out_mapping_file):
    """uses biomart to produce a mapping table
    """
    get_mappings = ("annot.ensembl_to_mappings.R {} {}").format(
        ensembl_id_file, out_mapping_file)
    run_shell_cmd(get_mappings)
    
    return None


def hgnc_list_to_ensembl_ids(hgnc_file, mapping_file, out_file):
    """get ensembl IDs attached to HGNC ID list
    """
    print hgnc_file
    print mapping_file
    print out_file

    # read in files
    hgnc_ids = pd.read_csv(hgnc_file, sep="\t", header=None)
    hgnc_ids.columns = ["hgnc_symbol"]
    mappings = pd.read_csv(mapping_file, sep="\t")

    # overlap
    mapped_ids = hgnc_ids.merge(mappings, on="hgnc_symbol", how="inner")

    # keep only hgnc and ensemb
    mapped_ids = mapped_ids[["ensembl_gene_id", "hgnc_symbol"]]
    mapped_ids = mapped_ids.drop_duplicates()
    mapped_ids.to_csv(out_file, sep="\t", header=False, index=False, compression="gzip")
    
    return


def build_ordered_tss_file(ordered_gene_list, gtf_file, out_bed):
    """Using an ensembl list, build a TSS file
    """

    # make a tmp gene only file
    work_dir = os.path.dirname(ordered_gene_list)
    tmp_gene_file = '{0}/{1}.genes.bed.gz'.format(work_dir, os.path.basename(gtf_file).split('.gtf')[0])
    make_gene_file_from_gtf = ("cat {0} | "
        "grep -P \"\tgene\t\" | "
        "awk -F '\"' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '\t' '{{ print $1\"\t\"$4\"\t\"$5\"\t\"$10\"\t\"$6\"\t\"$7 }}' | "
        "gzip -c > {1}").format(gtf_file, tmp_gene_file)
    if not os.path.isfile(tmp_gene_file):
        print make_gene_file_from_gtf
        run_shell_cmd(make_gene_file_from_gtf)

    # read in ordered gene list and tmp list and sort
    ordered_genes = pd.read_table(ordered_gene_list)
    ordered_genes.columns = ['geneid']
    gtf_genes = pd.read_table(tmp_gene_file)
    gtf_genes.columns = ['chr', 'start', 'stop', 'geneid', 'score', 'strand']

    merged = ordered_genes.merge(gtf_genes, on='geneid')
    merged = merged[gtf_genes.columns]

    print merged.shape
    print merged['geneid'][0:10]
    merged.to_csv(out_bed, sep='\t', header=False, index=False, compression='gzip')

    return None

def build_grouped_tss_file(ordered_gene_list, gtf_file, out_bed):
    """Using an ensembl list, build a TSS file. Requires the second column of the 
    gene list to have clusters
    """

    # make a tmp gene only file
    work_dir = os.path.dirname(ordered_gene_list)
    tmp_gene_file = '{0}/{1}.genes.bed.gz'.format(work_dir, os.path.basename(gtf_file).split('.gtf')[0])
    make_gene_file_from_gtf = ("cat {0} | "
        "grep -P \"\tgene\t\" | "
        "awk -F '\"' '{{ print $1\"\t\"$2 }}' | "
        "awk -F '\t' '{{ print $1\"\t\"$4\"\t\"$5\"\t\"$10\"\t\"$6\"\t\"$7 }}' | "
        "gzip -c > {1}").format(gtf_file, tmp_gene_file)
    if not os.path.isfile(tmp_gene_file):
        print make_gene_file_from_gtf
        run_shell_cmd(make_gene_file_from_gtf)

    # read in ordered gene list and tmp list and sort
    ordered_genes = pd.read_table(ordered_gene_list)
    ordered_genes.columns = ['geneid', 'cluster']
    gtf_genes = pd.read_table(tmp_gene_file)
    gtf_genes.columns = ['chr', 'start', 'stop', 'geneid', 'score', 'strand']

    merged = ordered_genes.merge(gtf_genes, on='geneid')

    # for each cluster, save out to file and add a hash
    clusters = sorted(list(set(merged['cluster'])))
    for cluster_idx in range(len(clusters)):
        cluster = clusters[cluster_idx]
        print cluster

        cluster_bed = merged.loc[merged['cluster'] == cluster][gtf_genes.columns]
        cluster_bed.to_csv(out_bed, sep='\t', mode='a', header=False, index=False, compression='gzip')        

        if cluster_idx != len(clusters)-1:
            with gzip.open(out_bed, 'a') as fp:
                fp.write('#\n')

    return None


def run_homer_on_clusters(in_file, cluster_nums, prefix, background_bed):
    '''
    Homer wrapper
    '''
    
    # extract relevant clusters
    for cluster_num in cluster_nums:
        cluster_file = '{0}.cluster-{1}.bed.gz'.format(prefix, cluster_num)
        if not os.path.isfile(cluster_file):
            extract_cluster = ("cat {0} | "
                               "awk -F '\t' '{{ if ($13=={1}) {{ print $0 }} }}' | "
                               "awk -F '\t' '{{ print $1\"\t\"$2\"\t\"$3}}' |  "
                               "gzip -c > {2}").format(in_file, cluster_num, cluster_file)
            print extract_cluster
            os.system(extract_cluster)

    # Then merge these into one file of regions
    interesting_regions_bed = '{}.merged.bed'.format(prefix)
    if not os.path.isfile(interesting_regions_bed):
        cluster_files = glob.glob('{}.cluster*.bed.gz'.format(prefix))
        merge_bed = ("zcat {0} | "
                     "sort -k1,1 -k2,2n "
                     "> {1}").format(' '.join(cluster_files), interesting_regions_bed)
        print merge_bed
        os.system(merge_bed)
        
    # run homer
    out_dir = '{}.homer'.format(prefix)
    if not os.path.isdir(out_dir):
        run_homer = ("findMotifsGenome.pl {0} hg19 {1} "
                     "-bg {2} "
                     "-p {3} "
                     "-nomotif").format(interesting_regions_bed,
                                        out_dir,
                                        background_bed,
                                        20)
        print run_homer
        os.system(run_homer)

    return None


def build_tss_file(gene_list, gtf_file, out_bed):
    '''
    Using an ensembl list, build a TSS file
    '''

    # make a tmp gene only file
    work_dir = os.path.dirname(out_bed)
    tmp_gene_file = '{0}/{1}.genes.bed.gz'.format(work_dir, os.path.basename(gtf_file).split('.gtf')[0])
    make_gene_file_from_gtf = ("cat {0} | "
                               "grep -P \"\tgene\t\" | "
                               "awk -F '\"' '{{ print $1\"\t\"$2 }}' | "
                               "awk -F '\t' '{{ print $1\"\t\"$4\"\t\"$5\"\t\"$10\"\t\"$6\"\t\"$7 }}' | "
                               "awk -F '\t' -v OFS='\t' '{{ gsub(/\..*/, \"\", $4); print }}' | "
                               "gzip -c > {1}").format(gtf_file, tmp_gene_file)
    if not os.path.isfile(tmp_gene_file):
        print make_gene_file_from_gtf
        os.system(make_gene_file_from_gtf)

    # read in ordered gene list and tmp list and sort
    genes = pd.read_table(gene_list)
    genes.columns = ['geneid']
    gtf_genes = pd.read_table(tmp_gene_file)
    gtf_genes.columns = ['chr', 'start', 'stop', 'geneid', 'score', 'strand']

    merged = genes.merge(gtf_genes, on='geneid')
    merged = merged[gtf_genes.columns]

    print merged.shape
    print merged['geneid'][0:10]
    merged.to_csv(out_bed, sep='\t', header=False, index=False, compression='gzip')

    return None

def get_tss_atac_regions(gene_id_file, tss_file, atac_regions_file, out_dir):
    '''
    Given a list of gene ids, get the BED of TSS regions that belong to 
    gene ids. Then intersect with ATAC regions to get open promoters
    Returns: tss point file and tss+ATAC region file
    '''

    prefix = '{0}/{1}'.format(out_dir, os.path.basename(gene_id_file).split('.txt')[0])

    tss_point_tmp = '{}.tss.tmp.gz'.format(prefix)
    tss_point_bed = '{}.tss.bed.gz'.format(prefix)
    tss_atac_bed = '{}.tss.atac.bed.gz'.format(prefix)

    # get BED of TSS regions that belong to gene id
    gene_ids = []
    with open(gene_id_file, 'r') as fp:
        for line in fp:
            gene_ids.append(line.strip())

    gene_id_set = set(gene_ids)

    # get TSS file
    with gzip.open(tss_point_tmp, 'w') as out:
        with gzip.open(tss_file, 'r') as fp:
            for line in fp:
                fields = line.strip().split('\t')

                gene_id = fields[3].split('.')[0]
                fields[3] = gene_id
                
                if gene_id in gene_id_set:
                    out.write('{}\n'.format('\t'.join(fields)))

    sort_bed = ("zcat {0} | "
                "sort -k1,1 -k2,2n | "
                "gzip -c > {1}").format(tss_point_tmp, tss_point_bed)
    print sort_bed
    os.system(sort_bed)
    os.system('rm {}'.format(tss_point_tmp))

                    
    # get ATAC file
    intersect = ("bedtools intersect -u -a {0} -b {1} | "
                 "sort -k1,1 -k2,2n | "
                 "gzip -c > {2}").format(atac_regions_file, tss_point_bed, tss_atac_bed)
    print intersect
    os.system(intersect)
    
    return tss_point_bed, tss_atac_bed
