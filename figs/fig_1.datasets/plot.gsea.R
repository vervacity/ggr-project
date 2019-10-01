#!/usr/bin/env Rscript

# description: plot GSEA summary
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/results/rna/timeseries/gsea
library(ggplot2)
library(viridis)
library(ggsci)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
gsea_summary_file <- args[1]
line_plot_file <- "fig_1.gsea.pdf"

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
gsea_summary$timepoints <- gsub("d", "", gsea_summary$timepoints)

# get colors
ggr_colors <- get_ggr_timepoint_colors()
ggr_colors <- ggr_colors[3:length(ggr_colors)]
name_colors <- pal_nejm("default")(3)
ggr_colors <- c(ggr_colors, name_colors)

# plot as line plots
ggplot(gsea_summary, aes(x=timepoints, y=NES, group=NAME)) +
    ggtitle("Gene set enrichments") +
    geom_line(aes(colour=NAME)) +
    geom_point(shape=21, stroke=1, fill="white", aes(colour=timepoints), show.legend=FALSE) +
    labs(x="Timepoint (day)") +
    theme_bw() +
    theme(
        aspect.ratio=1,
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.text=element_text(size=6),
        axis.text.x=element_text(angle=60, hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=5),
        legend.justification=c(0,1),
        legend.position=c(0,1),
        legend.background=element_blank(),
        legend.key.size=unit(0.01,"in")) +
    scale_color_manual(values=ggr_colors) +
    scale_y_continuous(limits=c(3.5, 6), expand=c(0,0)) +

ggsave(line_plot_file, height=2, width=2.5, useDingbats=FALSE)
