#!/usr/bin/env Rscript

# generate mapping table for ensembl id list

library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
mapping_file <- args[2]

# read table
ensembl_ids <- read.table(gzfile(in_file), header=TRUE, sep=",")

# get the annotations
grch37 <- useMart(
    biomart="ENSEMBL_MART_ENSEMBL",
    host="grch37.ensembl.org",
    path="/biomart/martservice",
    dataset="hsapiens_gene_ensembl")

# If looking for specific names, can use attributes
# grch37_attr <- listAttributes(grch37)
# grep('.*ntrez.*', grch37_attr$name)

conversion_table <- getBM(
    attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"),
    filters='ensembl_gene_id',
    values=ensembl_ids,
    mart=grch37)

write.table(
    conversion_table,
    file=gzfile(mapping_file),
    col.names=TRUE,
    row.names=FALSE,
    quote=FALSE,
    sep='\t')
