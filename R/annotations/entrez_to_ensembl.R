#!/usr/bin/env Rscript

# convert entrez to ensembl

library(biomaRt)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
conversion_file <- args[2]
id_file <- args[3]

# read table
fantom5_tfs <- read.table(in_file, header=TRUE, sep=",")

# get the annotations
grch37 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# If looking for specific names, can use attributes
# grch37_attr <- listAttributes(grch37)
# grep('.*ntrez.*', grch37_attr$name)

# make conversion table
conversion_table <- getBM(attributes=c('entrezgene', 'hgnc_symbol', 'ensembl_gene_id'), filters='entrezgene', values=fantom5_tfs$EntrezGene, mart=grch37)
write.table(conversion_table, file=gzfile(conversion_file), col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

write.table(conversion_table$ensembl_gene_id, file=gzfile(id_file), col.names=FALSE, row.names=FALSE, quote=FALSE)



