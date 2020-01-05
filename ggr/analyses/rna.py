# Description: RNA-specific analyses


import os
import glob
import gzip
import pandas as pd
import numpy as np

from ggr.analyses.filtering import filter_for_ids
from ggr.util.utils import run_shell_cmd


def filter_for_hgnc_ids(mat_file, keep_ids_file, mapping_file, out_file):
    """given an hgnc id list, using the mapping file to get ensembl ids and
    return the filtered matrix
    """
    # first, do the mapping
    hgnc_ids = pd.read_csv(keep_ids_file, header=None).iloc[:,0].values.tolist()
    mappings = pd.read_csv(mapping_file, sep="\t", header=0, index_col=False)
    matches = mappings[mappings["hgnc_symbol"].isin(hgnc_ids)]
    matches = matches.drop("entrezgene_id", axis=1)
    matches = matches.set_index("ensembl_gene_id")

    # then filter mat file
    data = pd.read_csv(mat_file, sep="\t", header=0, index_col=0)
    data_filt = matches.merge(data, how="inner", left_index=True, right_index=True).reset_index()

    # save out
    data_filt.to_csv(out_file, sep="\t", compression="gzip", header=True, index=False)
    
    return


def make_rsem_matrix(quant_files, out_file, colname="expected_count"):
    """Given a list of RSEM files, extract the desired column and make matrix
    column for counts: expected_counts
    column for tpm: TPM
    """
    signal_mat = pd.DataFrame()
    for quant_file in quant_files:
        fields = quant_file.split("/")[-2].split(".")
        prefix = '{0}_{1}'.format(fields[0].split('-')[1], fields[4])
        
        data = pd.read_table(quant_file, sep="\t")
        signal_mat[prefix] = data[colname]
        
    signal_mat["gene_id"], _ = data["gene_id"].str.split(".", 1).str
    signal_mat = signal_mat.set_index("gene_id")
    if colname == "expected_count":
        signal_mat = signal_mat.astype(int)
    signal_mat.to_csv(out_file, compression="gzip", sep='\t')
    
    return None


def threshold_empirically_off_genes(mat_file, out_file, thresh=1.0):
    """Given a matrix file, get a distribution and determine empirical threshold
    """
    data = pd.read_table(mat_file, index_col=0)
    data_thresh = data.loc[~(data < thresh).all(axis=1)]

    data_thresh.to_csv(out_file, compression="gzip", sep='\t')
    
    return None


def get_gene_sets_from_scRNA_dataset(
        sc_rna_file,
        out_file,
        pval_cutoff=0.05,
        logfc_cutoff=0.15): # 0.15
    """get results from Andrew's analysis file
    """
    diff_results = pd.read_csv(sc_rna_file, sep="\t")

    # get sig up and write out to file
    diff_sig_up = diff_results[
        (diff_results["p_val_adj"] < pval_cutoff) &
        (diff_results["avg_logFC"] > logfc_cutoff)
    ]

    with open(out_file, "w") as fp:
        fp.write("{}\t{}\n".format(
            "SINGLECELL_DIFFERENTIATED",
            "\t".join(list(diff_sig_up.index.values))))

    # get sig down and write out to file
    diff_sig_down = diff_results[
        (diff_results["p_val_adj"] < pval_cutoff) &
        (diff_results["avg_logFC"] < -logfc_cutoff)
    ]

    with open(out_file, "a") as fp:
        fp.write("{}\t{}\n".format(
            "SINGLECELL_PROGENITOR",
            "\t".join(list(diff_sig_down.index.values))))
        
    return


def run_gsea_on_series(
        deseq_files,
        genesets_file,
        out_file,
        id_conversion_file=None,
        gsea_jar_file="/users/dskim89/software/gsea/3.0/gsea-3.0.jar",
        tmp_dir="."):
    """gsea on sequence, aggregated to one output file
    """
    results_df = None
    for file_idx in range(len(deseq_files)):
        deseq_file = deseq_files[file_idx]
        deseq_file_prefix = os.path.basename(
            deseq_file).split("_over")[0].split("rna.")[1]
        
        # build rank file
        tmp_rank_file = "{}/{}.rnk".format(
            tmp_dir,
            os.path.basename(deseq_file).split(".txt")[0])
        build_gsea_rank_file(
            deseq_file, tmp_rank_file,
            id_conversion_file=id_conversion_file)

        # run GSEA
        deseq_out_dir = "{}/{}".format(
            tmp_dir,
            os.path.basename(deseq_file).split(".txt")[0])
        gsea_results_files = run_gsea_on_mat(
            tmp_rank_file, genesets_file, deseq_out_dir,
            gsea_jar_file=gsea_jar_file)

        # read in data and save to aggregate matrix
        # want to save p-vals and normalized enrichment scores
        for gsea_file_idx in range(len(gsea_results_files)):
            gsea_results_file = gsea_results_files[gsea_file_idx]
            gsea_results = pd.read_csv(gsea_results_file, sep="\t")
            gsea_results = gsea_results[["NAME", "NES", "FDR q-val"]]
            gsea_results["timepoints"] = deseq_file_prefix

            if (file_idx == 0) and (gsea_file_idx == 0):
                results_df = gsea_results.copy()
            else:
                results_df = pd.concat([results_df, gsea_results], axis=0)

    # save it all out
    results_df.to_csv(out_file, sep="\t", index=False, compression="gzip")

    # and plot
    dot_plot_file = "{}.dot_summary.pdf".format(out_file.split(".txt")[0])
    line_plot_file = "{}.line_summary.pdf".format(out_file.split(".txt")[0])
    plot_cmd = "plot.gsea_summary.R {} {} {}".format(
        out_file, dot_plot_file, line_plot_file)
    run_shell_cmd(plot_cmd)
    
    return None


def build_gsea_rank_file(deseq_results_file, out_file, id_conversion_file=None):
    """get a rank file from a DESeq2 results file
    """
    data = pd.read_csv(deseq_results_file, sep="\t")

    # use conversion if needed
    if id_conversion_file is not None:
        id_conversions = pd.read_csv(id_conversion_file, sep="\t")
        id_conversions = id_conversions[["ensembl_gene_id", "hgnc_symbol"]]
        id_conversions = id_conversions.set_index("ensembl_gene_id")
        id_conversions = id_conversions[~pd.isna(id_conversions["hgnc_symbol"])]
        data = data.merge(id_conversions, left_index=True, right_index=True)

    # extract just HGNC and logFC
    data = data[["hgnc_symbol", "log2FoldChange"]]
    data.to_csv(out_file, sep="\t", header=False, index=False)
    
    return None


def run_gsea_on_mat(
        rank_file, genesets_file, out_dir,
        min_geneset_size=15,
        max_geneset_size=1000,
        gsea_jar_file="/users/dskim89/software/gsea/3.0/gsea-3.0.jar"):
    """given a data mat, split to each column to run through prerank
    """
    #chip_file = "gseaftp.broadinstitute.org://pub/gsea/annotations/ENSEMBL_human_gene.chip"
    example_cmd = (
        "java -cp {} -Xmx1014m xtools.gsea.GseaPreranked "
        "-scoring_scheme classic "
        "-rnk {} "
        "-gmx {} "
        "-set_max {} "
        "-set_min {} "
        "-out {} "
        "-gui false").format(
            gsea_jar_file,
            rank_file,
            genesets_file,
            max_geneset_size,
            min_geneset_size,
            out_dir)
    if not os.path.isdir(out_dir):
        print example_cmd
        run_shell_cmd(example_cmd)

    # figure out which file has the results (neg and pos) and return
    gsea_results_files = sorted(glob.glob(
        "{}/*/gsea*xls".format(out_dir)))
    assert len(gsea_results_files) == 2
    
    return gsea_results_files


def filter_clusters_for_tfs(
        cluster_file,
        out_cluster_file,
        mapping_file,
        keep_list):
    """Filter cluster file for id strings in list
    """
    keep_list = [str(keep) for keep in keep_list]

    # read in mapping file
    mappings = {}
    with gzip.open(mapping_file, "r") as fp:
        for line in fp:
            fields = line.strip().split("\t")
            mappings[fields[0]] = fields[1]
    
    # open files and write out
    with open(cluster_file, "r") as fp:
        with open(out_cluster_file, "w") as out:

            for line in fp:
                fields = line.strip().split("\t")
                if "cluster" in fields[0]:
                    out.write(line)
                    continue

                current_id_string = fields[1]
                try:
                    mapped_symbol = mappings[current_id_string]
                except:
                    continue
                for keep_id in keep_list:
                    if keep_id in mapped_symbol:
                        out.write("{}\t{}\n".format(line.strip(), mapped_symbol))
                        
    return None


def filter_clusters_for_tss(
        cluster_file,
        out_prefix,
        tss_bed_file):
    """given cluster file, build TSS files for each cluster
    """
    # read in TSS file to dict
    geneid_to_tss = {}
    with gzip.open(tss_bed_file, "r") as fp:
        for line in fp:
            fields = line.strip().split()
            geneid_to_tss[fields[3]] = line
    
    # open files and write out
    with open(cluster_file, "r") as fp:
        for line in fp:
            if line.startswith("cluster"):
                continue
            
            fields = line.strip().split()
            cluster_file = "{}.cluster-{}.tmp.bed.gz".format(out_prefix, fields[0])
            try:
                with gzip.open(cluster_file, "a") as out:
                    out.write(geneid_to_tss[fields[1]])
            except:
                print "did not find {}".format(fields[1])

    # sort 
    cluster_tmp_files = glob.glob("{}*tmp.bed.gz".format(out_prefix))
    for cluster_tmp_file in cluster_tmp_files:
        cluster_file = "{}.bed.gz".format(cluster_tmp_file.split(".tmp")[0])
        sort_file = "zcat {} | sort -k1,1 -k2,2n | gzip -c > {}".format(
            cluster_tmp_file,
            cluster_file)
        run_shell_cmd(sort_file)

    # delete tmp files
    run_shell_cmd("rm {}*tmp*".format(out_prefix))
    
    return None


def get_neighborhood(
        bed_file,
        other_bed_file,
        out_bed_file,
        extend_len,
        chromsizes):
    """given a bed file and extend distance, call back regions
    from other bed file that are in that neighborhood
    """
    neighborhood = (
        "bedtools slop -i {0} -g {1} -b {2} | "
        "bedtools intersect -u -a {3} -b stdin | "
        "gzip -c > {4}").format(
            bed_file,
            chromsizes,
            extend_len,
            other_bed_file,
            out_bed_file)
    run_shell_cmd(neighborhood)

    return None


def get_nearest_regions(
        bed_file,
        other_bed_file,
        out_bed_file,
        opt_string=""):
    """use bedtools closest to get nearest region
    """
    closest = (
        "bedtools closest {0} -a {1} -b {2} | "
        "awk -F '\t' '{{ print $7\"\t\"$8\"\t\"$9\"\t\"$4\"\t\"$10\"\t\"$6 }}' | "
        "gzip -c > {3}").format(
            opt_string,
            bed_file,
            other_bed_file,
            out_bed_file)
    run_shell_cmd(closest)

    return None


def get_highly_expressed_genes(
        rna_mat_file,
        gene_list_file,
        out_file,
        percentile=90):
    """select genes that are most highly expressed and return gene list
    """
    # read in
    data = pd.read_table(rna_mat_file, index_col=0)
    gene_list = pd.read_table(gene_list_file, index_col=0)
    
    # filter
    data = data[data.index.isin(gene_list["id"])]

    # and then reduce
    data["max"] = data.max(axis=1)
    
    # and select top
    k = int(((100 - percentile) / 100.) * data.shape[0])
    data = data.nlargest(k, "max")

    # and save out gene list
    data.to_csv(out_file, columns=[], header=False, compression="gzip")
    
    return None


def convert_gene_list_to_tss(gene_list, tss_file, out_file):
    """given a gene list, get TSS
    """
    # read in TSS file to dict
    geneid_to_tss = {}
    with gzip.open(tss_file, "r") as fp:
        for line in fp:
            fields = line.strip().split()
            geneid_to_tss[fields[3]] = line
            
    # open files and write out
    tmp_file = "{}.tmp.bed.gz".format(out_file.split(".bed")[0])
    with gzip.open(gene_list, "r") as fp:
        for line in fp:
            fields = line.strip().split()
            try:
                with gzip.open(tmp_file, "a") as out:
                    out.write(geneid_to_tss[fields[0]])
            except:
                print "did not find {}".format(fields[0])
            
    # sort
    sort_file = "zcat {} | sort -k1,1 -k2,2n | gzip -c > {}".format(
        tmp_file, out_file)
    run_shell_cmd(sort_file)

    # delete tmp files
    run_shell_cmd("rm {}".format(tmp_file))

    return None


def add_clusters_to_tss_file(cluster_file, tss_file, out_file):
    """
    """
    # read in files
    tss_data = pd.read_csv(tss_file, sep="\t", header=None)
    cluster_data = pd.read_csv(cluster_file, sep="\t")

    # intersect
    tss_w_clusters = tss_data.merge(cluster_data, how="left", left_on=3, right_on="id")

    # clean up
    tss_w_clusters = tss_w_clusters.drop("id", axis=1)
    tss_w_clusters = tss_w_clusters.fillna(-1)
    tss_w_clusters[3] = "gene_id=" + tss_w_clusters[3].map(str) + ";traj=" + tss_w_clusters["cluster"].astype(int).map(str)
    tss_w_clusters = tss_w_clusters.drop("cluster", axis=1)

    # save out
    tss_w_clusters.to_csv(out_file, sep="\t", compression="gzip", header=False, index=False)
    
    return


def run_rdavid(gene_list, background_gene_list, out_dir):
    """Run DAVID using R
    """
    rdavid = ("bioinformatics.go.rdavid.R {} {} {}").format(
        gene_list,
        background_gene_list,
        out_dir)
    run_shell_cmd(rdavid)
    
    return None

