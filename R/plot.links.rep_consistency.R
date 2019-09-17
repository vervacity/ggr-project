#!/bin/usr/env Rscript

library(ggplot2)


args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
out_file <- args[2]

# read data
data <- read.table(gzfile(in_file), sep="\t", header=TRUE)
print(head(data))

# filter rep2
#data <- data[data$r2 > 10,]


#thresh <- 1
#data <- data[(data$r1 > thresh | data$r2 > thresh),]


# normalize and filter?
#data <- sweep(data, 2, colSums(data), "/") * 1e4 # 1e4

#data[data < 0] <- 0
#thresh <- 0
#data <- data[(data$r1 > thresh & data$r2 > thresh),]
print(dim(data))

# log normal scale
#data <- log10(data)

# plot
ggplot(data, aes(x=rep1, y=rep2)) + geom_point()
ggsave(out_file)

