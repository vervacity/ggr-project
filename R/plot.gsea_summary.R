#!/usr/bin/env Rscript

# description: plot GSEA summary
library(ggplot2)
library(viridis)

# args
args <- commandArgs(trailingOnly=TRUE)
gsea_summary_file <- args[1]
plot_file <- args[2]

# read in data
gsea_summary <- read.table(gsea_summary_file, sep="\t", header=TRUE)

# build sig (for size)
gsea_summary$sig <- gsea_summary$FDR.q.val < 0.25
gsea_summary$sig[gsea_summary$sig] <- 2
gsea_summary$sig[gsea_summary$sig !=2 ] <- 1.5

ggplot(gsea_summary, aes(x=timepoints, y=NAME)) +
    geom_point(aes(colour=NES, size=sig)) +
    scale_color_viridis(limits=c(-3,3))
        
    #scale_colour_gradient2(
    #    low="blue",
    #    mid="white",
    #    high="red",
    #    guide="colorbar",
    #    limits=c(-3,3))
ggsave(plot_file, width=14, height=7)

# build reduced set with most interesting

