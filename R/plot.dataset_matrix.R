#!/usr/bin/env Rscript


# description - plot out the datasets
library(ggplot2)
library(reshape2)
library(extrafont)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
dataset_file <- args[1]
plot_file <- args[2]

# read in file
data <- read.table(dataset_file, sep="\t", header=TRUE)
day_labels <- gsub("day.", "day ", colnames(data))[2:length(colnames(data))]

# adjustments
data <- data[data[,1] != "RAMPAGE",]
data <- data[data[,1] != "RNA-seq (polyA+)",]
data <- data[data[,1] != "RNA-seq (Ribo-Zero)",]
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
    size=1)

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
    size=3,
    stroke=1) +
    scale_color_manual(values=palette)

# cleanup
p <- p + 
    ylab(NULL) + 
    xlab("Timepoint (days)") + 
    theme_bw() + 
    theme(
        text=element_text(family="ArialMT", size=10),
        axis.text.x=element_text(angle=30, hjust=1, colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.title.x=element_text(size=8),
        panel.grid.minor=element_blank(),
        legend.position="none")

# save
ggsave(plot_file, height=1.75, width=4)
embed_fonts(plot_file)
