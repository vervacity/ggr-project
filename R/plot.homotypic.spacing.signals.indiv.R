#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
out_prefix <- args[2]
plot_file <- paste(out_prefix, ".pdf", sep="")

# read data
data <- read.table(gzfile(data_file), header=TRUE, sep="\t")

# clip
left_clip <- -130
right_clip <- 130
data <- data[data$position > left_clip,]
data <- data[data$position < right_clip,]

# blank out zero point
data <- data[data$position != 0,]

# signal specific adjustments
adjust_labels <- function(x) x
point_colour <- "black"

if (grepl("ATAC_SIGNALS", data_file)) {
    adjust_labels <- function(x) sprintf("%.1f", x)
    granularity <- 13
    assay_palette <- colorRampPalette(brewer.pal(9, "Blues"))(10)
    point_colour <- assay_palette[10]
} 

if (grepl("H3K27ac_SIGNALS", data_file)) {
    adjust_labels <- function(x) paste(x, ".", sep="")
    granularity <- 13
    assay_palette <- colorRampPalette(brewer.pal(9, "Reds"))(10)
    point_colour <- assay_palette[10]
} 


# plot
title <- unlist(strsplit(basename(data_file), ".", fixed=TRUE))[2:4]
title <- paste(title, collapse="; ")
ggplot(data, aes(x=position, y=score)) +
    geom_point(
        shape=20,
        stroke=0,
        alpha=0.1,
        colour=point_colour,
        position=position_jitter(w=1, h=0)) +
    labs(title=title, x="Position relative to anchor motif", y="Signal (FC)") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5), # 5
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(size=0.115),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.1, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6),
        legend.position=c(0.9, 0.9),
        strip.background=element_blank(),
        strip.text=element_blank()) +
        scale_x_continuous(
            limits=c(left_clip, right_clip),
            breaks=seq(left_clip, right_clip, 10),
            expand=c(0,0)) +
        scale_y_continuous(
            labels=adjust_labels)
ggsave(plot_file, height=1, width=2)








