#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]

# load data
prefix <- strsplit(basename(data_file), ".txt.gz")[[1]]
data <- read.table(gzfile(data_file), sep="\t", header=TRUE, row.names=1)
data_melt <- melt(data, id.vars="counts")

# get means
data_means <- aggregate(
    data,
    by=list(data$counts),
    mean)
data_means$Group.1 <- NULL
data_means_melt <- melt(data_means, id.vars="counts")

# plot
# DON'T finetune here - just want to see data, and THEN choose good vignette!
plot_file <- paste(prefix, ".pdf", sep="")
ggplot(data_melt, aes(x=jitter(counts), y=value, colour=variable)) +
    #geom_jitter() +
    #geom_smooth()
    geom_line(data=data_means_melt, aes(x=counts, y=value, colour=variable))
ggsave(plot_file, useDingbats=FALSE)



