#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
out_prefix <- args[2]
plot_file <- paste(out_prefix, ".pdf", sep="")

# read data
data <- read.table(gzfile(data_file), header=TRUE, sep="\t")
data_melt <- melt(data, id.vars="position")

# clip
left_clip <- -50
right_clip <- 50
data_melt <- data_melt[data_melt$position > left_clip,]
data_melt <- data_melt[data_melt$position < right_clip,]

# blank out middle
middle_blank <- 5
data_melt <- data_melt[abs(data_melt$position) > middle_blank,]

# plot
ggplot(data_melt, aes(x=position, y=value, colour=variable)) +
    geom_line(size=0.230) + 
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.y=element_blank(),
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
            expand=c(0,0))
ggsave(plot_file, height=1.5, width=3)








