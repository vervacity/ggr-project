#!/usr/bin/env Rscript

# description - plot out the datasets
# data: /users/dskim89/git/ggr-project/ggr/data

library(ggplot2)
library(reshape2)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
dataset_file <- args[1] # "~/git/ggr-project/ggr/data/ggr_datasets.txt"
plot_file <- "fig_1.datasets.matrix.pdf"

# read in file
data <- read.table(dataset_file, sep="\t", header=TRUE)
day_labels <- gsub("day.", "day ", colnames(data))[2:length(colnames(data))]

# adjustments
data <- data[data[,1] != "RAMPAGE",]
data <- data[data[,1] != "RNA-seq (polyA+)",]
data <- data[data[,1] != "RNA-seq (Ribo-Zero)",]
data <- data[data[,1] != "CTCF ChIP-seq",]
data$day.3.5 <- NULL
data$day.4.0 <- NULL
data$day.5.5 <- NULL

# melt and adjust
data_melted <- melt(data)
colnames(data_melted) <- c("assay", "day", "harvested")
data_melted$assay <- factor(data_melted$assay, levels=rev(unique(data[,1])), ordered=TRUE)
data_melted$day <- gsub("day.", "", data_melted$day)

# get colors
palette <- get_ggr_timepoint_colors()

# plotting
p <- ggplot()

# plot lines
p <- p + geom_line(
    data=data_melted,
    aes(
        group=assay,
        x=day,
        y=assay),
    size=0.5)

# plot circles
p <- p + geom_point(
    data=subset(data_melted, harvested==1),
    aes(
        group=day,
        x=day,
        y=assay,
        colour=as.factor(day)),
    shape=21,
    fill="white",
    size=1.75,
    stroke=0.5) +
    scale_color_manual(values=palette)

# cleanup
p <- p + 
    ylab(NULL) + 
    xlab("Timepoint (day)") + 
    theme_bw() + 
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.230),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6, colour="black"),
        axis.text.x=element_text(size=6, angle=60, hjust=1, colour="black"),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="none")

# save
ggsave(plot_file, height=1.5, width=2.5, useDingbats=FALSE)
