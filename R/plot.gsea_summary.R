#!/usr/bin/env Rscript

# description: plot GSEA summary
library(ggplot2)
library(viridis)
library(ggsci)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
gsea_summary_file <- args[1]
dot_plot_file <- args[2]
line_plot_file <- args[3]

# read in data
gsea_summary <- read.table(gsea_summary_file, sep="\t", header=TRUE)

keep_terms <- c(
    #"SINGLECELL_PROGENITOR",
    #"SINGLECELL_DIFFERENTIATED",
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

# adjust names
gsea_summary$timepoints <- gsub("0$", ".0", gsea_summary$timepoints)
gsea_summary$timepoints <- gsub("5$", ".5", gsea_summary$timepoints)

# get colors
ggr_colors <- get_ggr_timepoint_colors()
ggr_colors <- ggr_colors[3:length(ggr_colors)]
name_colors <- pal_nejm("default")(3)
ggr_colors <- c(ggr_colors, name_colors)


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
ggplot(gsea_summary, aes(x=timepoints, y=NES, group=NAME)) +
    ggtitle("Gene set enrichments") +
    geom_line(aes(colour=NAME)) +
    #geom_smooth(se=FALSE, aes(colour=NAME)) +
    #geom_point(aes(fill=timepoints, colour=timepoints)) +
    geom_point(shape=21, stroke=1, fill="white", aes(colour=timepoints)) +
    #geom_point(aes(size=sig)) +
    theme_bw() +
    theme(
        text=element_text(size=5),
        aspect.ratio=1,
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.25),
        axis.title=element_text(margin=margin(0,0,0,0)),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        legend.title=element_blank(),
        legend.text=element_text(size=2),
        #legend.justification=c(1,0),
        #legend.position=c(1,0),
        legend.justification=c(0,1),
        legend.position=c(0,1),
        legend.background=element_blank(),
        legend.key.size=unit(0.01,"in")) +
        scale_color_manual(values=ggr_colors)


ggsave(line_plot_file, height=2, width=2.5, useDingbats=FALSE) # width 14

# plot just the best subset
