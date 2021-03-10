#!/usr/bin/env Rscript

library(ggplot2)

# description: plot distances

args <- commandArgs(trailingOnly=TRUE)
plot_prefix <- args[1]
distance_files <- args[2:length(args)]

# read in all data
for (distance_i in 1:length(distance_files)) {

    # read in and simplify
    distance_file <- distance_files[distance_i]
    assay <- strsplit(distance_file, ".", fixed=TRUE)[[1]]
    assay <- rev(assay)[3]
    file_data <- read.table(distance_file, sep="\t", header=FALSE)
    file_data <- data.frame(distance=file_data$V10)
    file_data$assay <- assay

    print(assay)
    print(dim(file_data[file_data$distance > 100000,])[1])
    print(dim(file_data[file_data$distance <= 100000,])[1])
    print(dim(file_data[file_data$distance <= 100000,])[1] / dim(file_data)[1])
    
    
    # merge
    if (distance_i == 1) {
        data <- file_data
    } else {
        data <- rbind(data, file_data)
    }

}

# plot out
ggplot(data, aes(x=distance, after_stat(density), colour=assay)) +
    #geom_density() +
    geom_freqpoly(binwidth=50) +
    theme_bw() +
    scale_x_continuous(limits=c(0,5000)) +
    scale_y_continuous(limits=c(0,0.002))        

plot_file <- paste(plot_prefix, ".pdf", sep="")
ggsave(plot_file)
