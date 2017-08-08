#!/usr/bin/env Rscript

library('DESeq2')

# Use DESeq2 size estimation and rlog transform on 
# count matrices. Requires a master count matrix
# on which to build the normalization models and then
# the model is applied to following matrices.

args <- commandArgs(trailingOnly=TRUE)
master_count_file <- args[1]
other_count_files <- args[2:length(args)]

master_counts <- read.table(
    gzfile(master_count_file),
    header=TRUE,
    row.names=1,
    sep='\t')

# Set up conditions table
conditions <- sub("_b.", "", colnames(master_counts))
col_data <- data.frame(condition=conditions)
rownames(col_data) <- colnames(master_counts)

# Set up DESeq object
dds <- DESeqDataSetFromMatrix(countData=master_counts,
                              colData=col_data,
                              design = ~ condition)

# Estimate size factors and dispersion
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# Perform rlog norm
data_rlog <- rlog(dds, blind=FALSE)

# transfer rlog norm to other count matrices
for (count_file_idx in 1:length(other_count_files)) {

    count_file <- other_count_files[count_file_idx]
    prefix <- unlist(strsplit(count_file, ".mat", fixed=TRUE))
    out_file <- paste(prefix, ".rlog.mat.gz", sep="")

    counts <- read.table(
        gzfile(count_file),
        header=TRUE,
        row.names=1,
        sep='\t')
    
    # set up conditions table
    conditions <- sub("_b.", "", colnames(counts))
    col_data <- data.frame(condition=conditions)
    rownames(col_data) <- colnames(counts)

    # set up DESeq object
    dds_out <- DESeqDataSetFromMatrix(
        countData=counts,
        colData=col_data,
        design = ~ condition)
    
    # transfer rlog norm
    mcols(dds_out)$dispFit <- mcols(dds)$dispFit
    betaPriorVar <- attr(data_rlog, "betaPriorVar")
    intercept <- mcols(data_rlog)$rlogIntercept
    data_out_rlog <- rlog(
        dds_out,
        blind=FALSE,
        intercept=intercept,
        betaPriorVar=betaPriorVar)

    # save out norm data
    write.table(assay(data_out_rlog), file=gzfile(out_file), quote=FALSE, sep='\t')

}
