#!/usr/bin/env Rscript

# quick script to take a data table and 
# convert gene names to HGNC

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
conversion_file <- args[2]
out_file <- args[3]

data <- read.table(gzfile(in_file), sep='\t', header=TRUE)
colnames(data)[1] <- 'ensembl_gene_id'

conversion_table <- read.table(gzfile(conversion_file), sep='\t', header=TRUE)

annotated_data <- merge(data, conversion_table, by='ensembl_gene_id', all.x=TRUE, all.y=FALSE, sort=FALSE)
write.table(annotated_data, gzfile(out_file), sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
