#!/usr/bin/env Rscript

# run gProfiler
# requires gene set and background gene set

library(gProfileR)

# folders files etc
args <- commandArgs(trailingOnly=TRUE)
gene_list_file <- args[1]
background_list_file <- args[2]
out_dir <- args[3]
header <- as.numeric(args[4]) == 1
ordered <- as.numeric(args[5]) == 1

# read in gene list and background gene list
if (header) {
    gene_list <- read.table(
        gzfile(gene_list_file), header=TRUE, stringsAsFactors=FALSE)[,1]
} else {
    gene_list <- read.table(
        gzfile(gene_list_file), header=FALSE, stringsAsFactors=FALSE)[,1]
}
background_list <- read.table(
    gzfile(background_list_file), header=TRUE, stringsAsFactors=FALSE)[,1]

prefix <- sub('\\.txt.gz$', "", basename(gene_list_file))
padj_cutoff <- 0.1

# run gProfileR
results <- gprofiler(
    gene_list,
    ordered_query=ordered,
    organism="hsapiens",
    #max_p_value=padj_cutoff,
    #correction_method="fdr",
    custom_bg=background_list)

# clean up
results <- results[results$domain != "tf",]
results <- results[results$domain != "mir",]

# save out
out_file <- paste(out_dir, "/", prefix, ".go_gprofiler.txt", sep="")
write.table(results, out_file, quote=FALSE, sep="\t")
