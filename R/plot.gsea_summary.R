#!/usr/bin/env Rscript

# description: plot GSEA summary
library(ggplot2)
library(viridis)

# args
args <- commandArgs(trailingOnly=TRUE)
gsea_summary_file <- args[1]
dot_plot_file <- args[2]
line_plot_file <- args[3]

# read in data
gsea_summary <- read.table(gsea_summary_file, sep="\t", header=TRUE)

keep_terms <- c(
    "SINGLECELL_PROGENITOR",
    "SINGLECELL_DIFFERENTIATED",
    #"GO_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
    #"GO_REGULATION_OF_STEM_CELL_PROLIFERATION",
    "GO_KERATINOCYTE_DIFFERENTIATION",
    "GO_KERATINIZATION",
    "GO_CORNIFIED_ENVELOPE")

if (TRUE) {
    gsea_summary <- gsea_summary[gsea_summary$NAME %in% keep_terms,]
}

# build sig (for size)
gsea_summary$sig <- gsea_summary$FDR.q.val < 0.25
gsea_summary$sig[gsea_summary$sig] <- 2
gsea_summary$sig[gsea_summary$sig !=2 ] <- 1.5

ggplot(gsea_summary, aes(x=timepoints, y=NAME)) +
    #geom_point(aes(colour=sig, size=NES)) +
    geom_point(aes(colour=NES, size=sig)) +
    #scale_color_viridis()
    scale_color_viridis(limits=c(-9,9))
    #scale_colour_gradient2(
    #    low="blue",
    #    mid="white",
    #    high="red",
    #    guide="colorbar",
    #    limits=c(-3,3))
ggsave(dot_plot_file, width=14, height=7)

# plot as line plots
ggplot(gsea_summary, aes(x=timepoints, y=NES, group=NAME, colour=NAME)) +
    geom_line() +
    geom_point(aes(size=sig)) +
ggsave(line_plot_file, width=42, height=7) # width 14

# plot just the best subset
