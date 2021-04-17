#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(reshape2)
library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(data_file, sep="\t", header=TRUE)

# filler data to make it a donut
filler_data <- data.frame(fill=0, count=1, name="center")
data <- rbind(data, filler_data)

# make NA to not show no data
data[data == 0] <- NA

# levels for ordering rings
ordered_levels <- c(
    "center", # nothing
    "atac_dynamic",
    #"rna_present",
    "H3K27ac_present",
    "H3K4me1_present",
    "H3K27me3_present",
    "rna_present")
data$name <- factor(data$name, levels=ordered_levels)

# plot
ggplot(data, aes(x=name, y=count, fill=fill)) +
    geom_bar(stat="identity") +
    coord_polar(theta="y") +
    theme_bw() +
    scale_fill_continuous(na.value=NA)
    #scale_fill_discrete(na.value=NA)
        
ggsave(plot_file)





