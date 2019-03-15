#!/user/bin/env python

import os
import glob
import gzip

import pandas as pd

from tronn.preprocess.variants import generate_new_fasta


def convert_gtex_to_vcf(
        variant_files,
        lookup_table,
        output_file,
        id_col=1,
        gene_col=2,
        filter_indels=True,
        tmp_dir="."):
    """convert GTEx files to VCF format
    1-index column choice
    """
    # set up filter cmd
    if filter_indels:
        filter_cmd = "| awk -F '_' '{{ if ( length($3) == 1 && length($4)==1 ) print $0 }}' "
    else:
        filter_cmd = ""
        
    # first read in variant files and keep as set
    tmp_all_variants_file = "{}/all_variants.tmp.txt".format(tmp_dir)
    if not os.path.isfile(tmp_all_variants_file):
        combine_variants = (
            "zcat {0} | "
            "awk -F '\t' '{{ print ${1}\"\t\"${2} }}' {3} | "
            "awk -F '.' '{{ print $1 }}' | "
            "sort | uniq "
            "> {4}").format(
                " ".join(variant_files),
                id_col,
                gene_col,
                filter_cmd,
                tmp_all_variants_file)
        print combine_variants
        os.system(combine_variants)

    # go through lookup table
    if not os.path.isfile(output_file):
        all_variants_df = pd.read_table(tmp_all_variants_file)
        variant_to_gene = dict(zip(
            all_variants_df.iloc[:,0],
            all_variants_df.iloc[:,1]))
        all_variants = set(all_variants_df.iloc[:,0])
        print len(all_variants)
        with open(output_file, "w") as out:
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            with open(lookup_table, "r") as fp:
                for line in fp:
                    fields = line.strip().split("\t")
                    variant_id = fields[2]
                    if variant_id in all_variants:
                        # write out
                        out_fields = [
                            "{}".format(fields[0]),
                            fields[1],
                            fields[6],
                            fields[3],
                            fields[4],
                            "1",
                            "1",
                            "gtex_id={};study=GTEx;study_type=eQTL;eGene={}".format(variant_id, variant_to_gene[variant_id])]
                        out.write("{}\n".format("\t".join(out_fields)))
                        
    return None


def convert_gwas_to_vcf(
        in_file,
        out_file,
        cleanup=True):
    """convert Yang's GWAS file to VCF
    """
    if not os.path.isfile(out_file):
        with open(out_file, "w") as out:
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

            # go through lines from input file
            with open(in_file, "r") as fp:
                for line in fp:
                    if line.startswith("Disease"):
                        fields = line.strip().split("\t")
                        continue

                    # read out line bits
                    fields = line.strip().split("\t")

                    # cleanup
                    if cleanup:
                        if len(fields[4]) > 1:
                            if "," in fields[4]:
                                fields[4] = fields[4][0]
                            else:
                                continue
                        if len(fields[5]) > 1:
                            if "," in fields[5]:
                                fields[5] = fields[5][0]
                            else:
                                continue
                    
                    # set up info field, keep lead snp disease and study id
                    metadata = "lead_snp={};study={};study_type=GWAS;reported_gene={};disease={}".format(
                        fields[17],
                        fields[19],
                        fields[24].replace(" ", ""),
                        fields[0].replace(" ", "_").replace("'", "").lower())

                    keep_fields = [
                        fields[2].replace("chr", ""),
                        fields[3],
                        fields[1],
                        fields[4],
                        fields[5],
                        "1",
                        "1",
                        metadata]
                    out.write("{}\n".format("\t".join(keep_fields)))

    return None


def merge_vcf_files(vcf_files, out_vcf_file):
    """merge vcf files, handling conflicts
    """
    # first cat together and sort?
    tmp_joint_file = "merge.tmp.vcf"
    if not os.path.isfile(tmp_joint_file):
        join_files = "cat {} | grep -v '^#' | sort -k1,1 -k2,2n | uniq > {}".format(
            " ".join(vcf_files), tmp_joint_file)
        print join_files
        os.system(join_files)
    
    # then go through and merge conflicts
    with open(out_vcf_file, "w") as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        num_overlaps = 0
        mismatch_alt = 0
        prev_id = ""
        prev_line = ""
        with open(tmp_joint_file, "r") as fp:
            for line in fp:
                fields = line.strip().split("\t")
                rs_id = fields[2]

                # fix for GTEX
                if rs_id == ".":
                    rs_id = dict([val.split("=") for val in fields[7].split(";")])["gtex_id"]
                
                if rs_id == prev_id:
                    # merge info
                    prev_fields = prev_line.strip().split("\t")
                    # check matches
                    if prev_fields[0] != fields[0]:
                        print "chromosome doesn't match!"
                    if prev_fields[1] != fields[1]:
                        print "pos doesn't match!"
                    if prev_fields[2] != fields[2]:
                        print "id doesn't match!"
                    if prev_fields[3] != fields[3]:
                        print "ref doesn't match!"
                    if prev_fields[4] != fields[4]:
                        print "alt doesn't match!"
                        mismatch_alt += 1
                    prev_fields[7] = "{}|{}".format(prev_fields[7], fields[7])
                    prev_line = "\t".join(prev_fields)
                    num_overlaps += 1
                else:
                    # write out
                    if len(prev_line) != 0:
                        out.write(prev_line)
                    # and replace
                    prev_id = rs_id
                    prev_line = line
                        

        # finally write out last line
        out.write(prev_line)

    print num_overlaps
    print mismatch_alt

    
    return None


def main():
    """take variants and generate useful files for downstream analysis
    """
    GWAS_DIR = "../GWAS"
    GTEX_DIR = "../GTEx"
    GGR_DIR = "/mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/data"
    
    # files
    # NOTE: the egenes file contains only the strongest pval hit AND also contains insignificant hits
    gtex_lookup_file = "{}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt".format(GTEX_DIR)
    gtex_skin_files = glob.glob("{}/GTEx_Analysis_v7_eQTL/Skin*signif_variant*".format(GTEX_DIR))
    gwas_expanded_file = "{}/causal-SNP-keratinocytes-all.txt".format(GWAS_DIR)
    gwas_vcf = "gwas.skin.vcf"
    atac_bed_file = "{}/ggr.atac.idr.master.bed.gz".format(GGR_DIR)
    orig_fasta = "/mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa"
    
    # convert GTex to VCF format
    gtex_vcf = "gtex.skin.vcf"
    convert_gtex_to_vcf(
        gtex_skin_files,
        gtex_lookup_file,
        gtex_vcf)
    
    # GWAS convert to VCF format
    convert_gwas_to_vcf(
        gwas_expanded_file,
        gwas_vcf)

    # merge GWAS file with GTEx file? and then handle conflicts appropriately
    # NOTE: too many conflicts, don't do this
    if False:
        merged_vcf = "gtex-gwas.skin.vcf"
        merge_vcf_files(
            [gtex_vcf, gwas_vcf],
            merged_vcf)

    # for each file, generate the ref/alt fastas and filter with ATAC
    vcf_files = [gtex_vcf, gwas_vcf]
    for vcf_file in vcf_files:
        # use seqtk to generate ref/alt fastas
        ref_fasta = "./{}.ref.fa".format(vcf_file.split(".vcf")[0])
        if not os.path.isfile(ref_fasta):
            generate_new_fasta(vcf_file, orig_fasta, ref_fasta, ref=True)
        alt_fasta = "./{}.alt.fa".format(vcf_file.split(".vcf")[0])
        if not os.path.isfile(alt_fasta):
            generate_new_fasta(vcf_file, orig_fasta, alt_fasta, ref=False)

        # and convert it into a BED file
        bed_file = "{}.bed.gz".format(vcf_file.split(".vcf")[0])
        if not os.path.isfile(bed_file):
            make_bed = (
                "cat {} | "
                "awk -F '\t' 'NR>1{{print \"chr\"$1\"\t\"$2\"\t\"$2+1\"\t\"$3\"\t1000\t.\"}}' | "
                "gzip -c > {}").format(vcf_file, bed_file)
            print make_bed
            os.system(make_bed)

        # filter bed file through ATAC
        bed_filt_file = "{}.atac_filt.bed.gz".format(bed_file.split(".bed")[0])
        if not os.path.isfile(bed_filt_file):
            filter_bed = (
                "bedtools intersect -wa -a {} -b {} | "
                "sort -k1,1 -k2,2n | "
                "gzip -c > {}").format(
                    bed_file, atac_bed_file, bed_filt_file)
            print filter_bed
            os.system(filter_bed)
    
    return


main()
