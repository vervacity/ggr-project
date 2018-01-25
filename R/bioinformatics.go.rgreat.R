#!/usr/bin/env Rscript

# Description: wrapper for rGREAT

library(rGREAT)
library(qvalue)

# args
args <- commandArgs(trailingOnly=TRUE)
bed_file <- args[1]
background_bed_file <- args[2]
prefix <- args[3]

# setup bed file
bed <- read.table(gzfile(bed_file), header=FALSE, row.names=NULL)
colnames(bed) <- c('chr', 'start', 'stop')

# setup background file
background_bed <- read.table(gzfile(background_bed_file), header=FALSE, row.names=NULL)
colnames(background_bed) <- c("chr", "start", "stop")

# create GREAT job
job <- submitGreatJob(bed, species='hg19', version="3.0.0")
print(job)

# get enrichments
tb <- getEnrichmentTables(job) # this is where you would get other tables if you wanted them
print(names(tb))

# threshold and save out
for (name in names(tb)) {
    table_file_name <- paste(prefix,
                             gsub(" ", "_", name),
                             "txt",
                             sep=".")
    print(table_file_name)
    
    enrichment_table <- tb[[name]]
    #q_obj <- qvalue(p=enrichment_table$Hyper_Raw_PValue)
    q_obj <- qvalue(p=enrichment_table$Binom_Raw_PValue)
    enrichment_table$hyper_q_vals <- q_obj$qvalues
    enrichment_table_thresh <- enrichment_table[(enrichment_table$hyper_q_vals < 0.05) &
                                                   (enrichment_table$Binom_Fold_Enrichment > 2),]
    enrichment_table_thresh_sorted <- enrichment_table_thresh[with(enrichment_table_thresh, order(hyper_q_vals)),]
    write.table(enrichment_table_thresh_sorted, table_file_name, quote=FALSE, sep='\t')
}
