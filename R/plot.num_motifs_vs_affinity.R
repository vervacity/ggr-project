#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(reshape2)
library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]
title <- args[3]

# read in data
data <- read.table(data_file, header=TRUE)

data$region <- NULL
data$start_summit <- NULL
data$start <- NULL

if ( any(grepl("extra", colnames(data))) ) {
    labels <- data$extra
    data$extra <- NULL
} else {
    labels <- rep(1, nrow(data))
}


data[data == -1] <- NA
keep <- rowSums(data, na.rm=TRUE) != 0
data <- data[keep,]
labels <- labels[keep]
print(dim(data))
data[data == 0] <- NA

num_motifs <- rowSums(data != 0, na.rm=TRUE)


my_median <- function(x) {
    return(median(x, na.rm=TRUE))
}
#affinities <- apply(data, 1, my_median)
affinities <- rowSums(data, na.rm=TRUE) / num_motifs

results <- data.frame(
    num_motifs=num_motifs,
    mean_affinities=affinities,
    labels=labels)

# ggplot
ggplot(results, aes(x=num_motifs, y=mean_affinities, group=num_motifs)) +
    #geom_violin()
    #geom_boxplot()
    geom_jitter(shape=20, alpha=0.5, size=1, stroke=0, aes(color=labels)) +
    labs(x="Number of motif instances in region",
         y="Average motif affinity score",
         title=title) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_text(vjust=-1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        #legend.position="bottom",
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0)) +
    scale_colour_distiller(direction=1, limits=c(0,1))
    #scale_y_continuous(limits=c(0, 0.025))
ggsave(plot_file, height=3, width=3)
