#!/usr/bin/env Rscript

# Description: calculates all pairwise comparisons
# in a dataset and returns the union of them all

# script to run DESeq2 on sequential timepoints
# importantly, it only runs the noise model between
# two timepoints, since there may be significant
# cell type changes across time (ie the noise model
# does not hold between two very different cell types)

# assumes header is of d00_b1 type and there are 
# bio reps in order

library(DESeq2)
library(stringr)

print('running DESeq2...')

# folders files etc
args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
prefix <- args[2]
fdr_cutoff <- as.numeric(args[3])
out_file <- args[4]
is_sequential_only <- args[5]


mkdir <- paste("mkdir -p",
               dirname(prefix),
               sep=" ")
system(mkdir)

# make count matrix
count_data <- as.matrix(read.table(
    gzfile(in_file),
    sep='\t',
    header=TRUE,
    row.names=1))

# keep track of up/down numbers for each comparison
differential_summary <- data.frame()

# run pairwise comparisons for all
for (idx1 in seq(1, ncol(count_data)-1, 2)) {
    for (idx2 in seq(1, ncol(count_data)-1, 2)) {

        # do not evaluate if same indices
        if (idx1 == idx2) { next }

        # if only running sequential pairwise,
        # only run if idx2 == idx1 + 2
        # ie, if they are next to each other
        if (is_sequential_only == "sequential") {
            if (idx2 != idx1 + 2) { next }
        }
        
        print(paste(idx1, idx2))
        
        # make smaller count matrix
        desired_columns <- c( idx1, idx1+1, idx2, idx2+1  )
        print(desired_columns)
        pairwise_count_data <- count_data[, desired_columns]

        # make the condition table
        conditions <- sub("_b.", "", colnames(pairwise_count_data))
        unique_conditions <- unique(conditions)
        t_baseline <- unique_conditions[1]
        t_compare <- unique_conditions[2]
        col_data <- data.frame(condition=conditions)
        rownames(col_data) <- colnames(pairwise_count_data)

        # do not evaluate if files exist
        sigresultsall_file <- paste(
            prefix, '.', t_compare, '_over_', t_baseline, '_sigResultsAll.txt', sep='')
        if (file.exists(sigresultsall_file)) { next }

        # make DESeq dataset
        dds <- DESeqDataSetFromMatrix(countData=pairwise_count_data,
                                  colData=col_data,
                                  design = ~ condition)

        # build DESeq model
        dds <- DESeq(dds)

        # get results and save out to tables
        res <- results(dds, contrast=c('condition', t_compare, t_baseline), alpha=fdr_cutoff)

        # remove NA and filter for cutoff
        res_noNA <- res[!is.na(res$padj),]
        res_filt <- res_noNA[res_noNA$padj < fdr_cutoff,]

        # separate up and down sets
        res_filt_up <- res_filt[res_filt$log2FoldChange > 0,]
        res_filt_down <- res_filt[res_filt$log2FoldChange < 0,]

        # write everything out
        write.table(
            res_noNA,
            file=gzfile(
                paste(prefix, '.', t_compare, '_over_', t_baseline, '_resultsAll.txt.gz', sep='')),
            quote=FALSE, sep='\t')

        write.table(
            res_filt,
            file=gzfile(
                paste(prefix, '.', t_compare, '_over_', t_baseline, '_sigResultsAll.txt.gz', sep='')),
            quote=FALSE, sep='\t')
        
        up_names <- str_split_fixed(rownames(res_filt_up), '[.]', 2)[,1] # for ENSG genes
        write.table(
            up_names,
            file=gzfile(
                paste(prefix, '.', t_compare, '_over_', t_baseline, '_sigResultsUp.txt.gz', sep='')),
            quote=FALSE, row.names=FALSE, col.names=FALSE,
            sep='\t')
        
        down_names <- str_split_fixed(rownames(res_filt_down), '[.]', 2)[,1] # for ENSG genes
        write.table(
            down_names,
            file=gzfile(
                paste(prefix, '.', t_compare, '_over_', t_baseline, '_sigResultsDown.txt.gz', sep='')),
            quote=FALSE, row.names=FALSE, col.names=FALSE,
            sep='\t')

        # save to differential summary
        compare_summary <- data.frame(
            compare=paste(t_baseline, "_to_", t_compare, sep=""),
            up=length(up_names),
            down=length(down_names))
        differential_summary <- rbind(differential_summary, compare_summary)
        print(differential_summary)
        
    }
}

# save out the differential summary
differential_summary_file <- paste(
    dirname(prefix), "/differential_summary.txt.gz", sep="")
write.table(
    differential_summary,
    file=gzfile(differential_summary_file),
    quote=FALSE, row.names=FALSE, col.names=TRUE,
    sep="\t")

# now merge sigResults to get list of full regions that are significant
merge_sigresults <- paste("zcat ",
                          dirname(prefix),
                          "/*sigResultsAll.txt.gz | ",
                          "awk -F '\t' '{ print $1 }' | ",
                          "grep -v baseMean | ",
                          "sort | ",
                          "uniq | ",
                          "gzip -c > ",
                          out_file,
                          sep="")
system(merge_sigresults)

