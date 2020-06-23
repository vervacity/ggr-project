#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(ggsci)


plot_fn <- function(data, plot_file) {

    p <- ggplot(data, aes(x=variable, y=value)) +
        geom_hline(yintercept=0, size=0.115, linetype="dashed") +
        geom_boxplot(
            size=0.230, width=0.7, outlier.shape=NA) +
        geom_point(
            size=0.5, alpha=0.5, show.legend=FALSE) +
        ylab("log2(FC)") +
        xlab("Timepoint") +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.margin=margin(5,10,1,5),
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
            legend.position="bottom")
    
    if (length(unique(data$pwms)) == 1) {
        p <- p + scale_color_aaas()
    } else {
        p <- p + scale_colour_brewer(palette="Paired", direction=-1)
    }
    
    ggsave(plot_file, height=1.5, width=2, useDingbats=FALSE)

}


# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]

# load data and plot
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)
plot_file <- paste(plot_prefix, ".pdf", sep="")
plot_fn(data, plot_file)
