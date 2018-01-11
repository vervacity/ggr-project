#!/usr/bin/env Rscript

library('DESeq2')


# quick script using DESeq2 size estimation and rlog transform

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
out_file <- args[2]

data <- read.table(
    gzfile(data_file),
    header=TRUE,
    row.names=1,
    sep='\t')

# Set up conditions table
conditions <- sub("_b.", "", colnames(data))
col_data <- data.frame(condition=conditions)
rownames(col_data) <- colnames(data)

# Set up DESeq object
if (length(unique(conditions)) != 1) {
    dds <- DESeqDataSetFromMatrix(countData=data,
                                  colData=col_data,
                                  design = ~ condition)
} else {
    dds <- DESeqDataSetFromMatrix(countData=data,
                                  colData=col_data,
                                  design = ~ 1)
}

# Estimate size factors and dispersion
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Perform rlog norm
data_rlog <- rlog(dds, blind=FALSE)

# Save out new data
write.table(assay(data_rlog), file=gzfile(out_file), quote=FALSE, sep='\t')
