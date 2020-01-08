#!/usr/bin/env Rscript

library(ggplot2)
library(ggsci)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read data
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)

# set up
data$labels <- paste(data$pwm1_clean, data$pwm2_clean, sep=",")
data <- data[data$interaction != "FAILED.NEGATIVE_ENDOG",]
data <- data[data$interaction != "FAILED.TRAJ_PATTERN",]

# plot
ggplot(data, aes(x=expected, y=actual, colour=interaction)) +
    geom_abline(intercept=0, slope=1, size=0.115, linetype="dashed") +
    geom_point() +
    #geom_text(aes(label=labels, vjust=1, hjust=1), size=1) +
    labs(x="Expected", y="Actual") +
    coord_fixed() +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,6,5),
        plot.title=element_text(size=8, margin=margin(b=1)),
            
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
        legend.position="bottom") +
     scale_colour_brewer(palette="Set2")

ggsave(plot_file, height=2, width=4, useDingbats=FALSE)

