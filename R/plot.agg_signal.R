#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
mat_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(mat_file, sep="\t", header=TRUE)

# plot
ggplot(data, aes(x=position, y=mean, group=timepoint, color=timepoint)) +
    geom_line()
ggsave(plot_file)
