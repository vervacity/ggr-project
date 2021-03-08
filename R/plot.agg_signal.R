#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
mat_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(mat_file, sep="\t", header=TRUE)

# plot
ggplot(data, aes(x=position, y=mean, group=timepoint, color=timepoint)) +
    geom_line() +
    labs(x="Position", y="Signal log(FC)") +
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
        axis.ticks.length=unit(0.01, "in"),
            
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.1, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.spacing.x=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5),
        legend.position="bottom")
    
ggsave(plot_file, height=2, width=2, useDingbats=FALSE)
