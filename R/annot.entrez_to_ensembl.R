#!/usr/bin/env Rscript

# convert entrez to ensembl
# specifically for fantom5 tfs

library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
id_file <- args[2]

# read table
fantom5_tfs <- read.table(in_file, header=TRUE, sep=",")

# get the annotations
grch37 <- useMart(
    biomart="ENSEMBL_MART_ENSEMBL",
    host="grch37.ensembl.org",
    path="/biomart/martservice",
    dataset="hsapiens_gene_ensembl")

# If looking for specific names, can use attributes
#grch37_attr <- listAttributes(grch37)
#print(grep('.*ntrez.*', grch37_attr$name))
#print(grch37_attr[c(55,56,57),])

# make conversion table
prefix <- basename(unlist(strsplit(in_file, "[.]"))[1])
out_dir <- dirname(id_file)
out_prefix <- paste(out_dir, prefix, sep="/")
conversion_file <- paste(out_prefix, "conversion_ids.txt.gz", sep=".")

conversion_table <- getBM(
    attributes=c('entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id'),
    filters='entrezgene_id',
    values=fantom5_tfs$EntrezGene,
    mart=grch37)
write.table(
    conversion_table,
    file=gzfile(conversion_file),
    col.names=TRUE,
    row.names=FALSE,
    quote=FALSE,
    sep='\t')

# write out simplified table
write.table(
    conversion_table$ensembl_gene_id,
    file=gzfile(id_file),
    col.names=FALSE,
    row.names=FALSE,
    quote=FALSE)



