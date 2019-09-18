#!/bin/usr/env Rscript

library(ggplot2)


args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
out_file <- args[2]

# read data
data <- read.table(gzfile(in_file), sep="\t", header=TRUE)
print(dim(data))

# plot
ggplot(data, aes(x=rep1, y=rep2)) + geom_point()
ggsave(out_file)

