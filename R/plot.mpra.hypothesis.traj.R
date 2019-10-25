#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]


# load data
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)
data_melt <- melt(data, id.vars="example_id")
data_melt$day <- gsub("_.+", "", data_melt$variable)
data_melt$rep <- gsub(".+_", "", data_melt$variable)

# plot

ggplot(data_melt, aes(x=day, y=value, colour=rep)) +
    geom_line(alpha=0.7, size=0.115, colour="gray", aes(group=example_id)) +
    geom_boxplot(position=position_dodge()) +
    geom_point(position=position_jitterdodge(jitter.width=0)) +
    theme_bw()        

ggsave(plot_file)





