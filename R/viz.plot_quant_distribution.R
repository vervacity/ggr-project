#!/usr/bin/env Rscript

library(reshape)
library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
threshold <- as.numeric(args[2])
out_file <- args[3]

# melt and use ggplot
data <- read.table(gzfile(in_file), header=TRUE, row.names=1)

data_nonzero <- data[rowSums(data)>0, ]
print(dim(data_nonzero))

data_melted <- melt(data_nonzero)
ggplot(data_melted, aes(x=value, colour=variable)) + geom_density() + geom_vline(xintercept=threshold)
ggsave(out_file)
