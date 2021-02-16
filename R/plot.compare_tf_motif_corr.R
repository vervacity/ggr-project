#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

data <- read.table(data_file)

ggplot(data, aes(x=Homer.ATAC, y=NN)) +
    geom_point(shape=21, fill=NA) +
    geom_abline(intercept=0, slope=1) +
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
        axis.title.y=element_text(vjust=1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in")) +
    scale_x_continuous(expand=c(0,0), limits=c(0, 1)) +
    scale_y_continuous(expand=c(0,0), limits=c(0, 1))

 
ggsave(plot_file, height=2, width=2, useDingbats=FALSE)
