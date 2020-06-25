#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
val_file <- args[1]
plot_file <- args[2]

# read in data and get gene vector
vals <- read.table(val_file, sep="\t", header=TRUE)
vals$labels <- ""
vals$labels[vals$timepoint == "d0.0"] <- "d0"
vals$labels[vals$timepoint == "d6.0"] <- "d6"

vals$ATAC <- vals$ATAC - min(vals$ATAC)
vals$predicted <- vals$predicted - min(vals$predicted)

vals$ATAC <- vals$ATAC / max(vals$ATAC)
vals$predicted <- vals$predicted / max(vals$predicted)


vals_melt <- melt(vals, id.vars=c("timepoint", "labels"))
vals_melt$variable <- as.character(vals_melt$variable)
vals_melt$variable[vals_melt$variable == "predicted"] <- "NN predicted"

# plot
ggplot(vals_melt, aes(x=timepoint, y=variable, fill=value)) +
    geom_tile(show.legend=FALSE) +
    theme_bw() +
    theme(
        plot.margin=margin(1,0,0,0),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.230), #element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(size=5, vjust=1.5, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=5, margin=margin(0,0,0,0))) +
    scale_x_discrete(labels=vals_melt$labels, expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_distiller(palette="Blues", direction=1)
ggsave(plot_file, height=0.25, width=1.0, useDingbats=FALSE)
