#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)


plot_fn <- function(data, plot_file) {

    ggplot(data, aes(x=day, y=value, colour=combos)) +
        geom_hline(yintercept=0, size=0.115) +
        #geom_line(alpha=0.7, size=0.115, colour="gray", aes(group=example_id)) +
        #geom_violin(
        #    position=position_dodge(width=0.9), size=0.230, trim=FALSE, scale="width", show.legend=FALSE) +
        geom_boxplot(
            position=position_dodge(width=0.9), size=0.230, width=0.7, outlier.shape=NA, show.legend=FALSE) +
        geom_point(
            position=position_jitterdodge(dodge.width=0.9, jitter.width=0.05),
            size=0.5, alpha=0.5, show.legend=FALSE) +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.margin=margin(5,20,1,0),
            plot.title=element_text(size=8, margin=margin(b=1)),
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
            axis.ticks.length=unit(0.01, "in"))
        
    ggsave(plot_file, height=2, width=3, useDingbats=FALSE)

}


# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]

# load data
data <- read.table(gzfile(data_file), sep="\t", header=TRUE)

# plot just 0,0 first
plot_file <- paste(plot_prefix, ".endog_only.pdf", sep="")
data_filt <- data[data$combos == "0,0",]
plot_fn(data_filt, plot_file)

# plot with null
plot_file <- paste(plot_prefix, ".endog_and_null.pdf", sep="")
data_filt <- data[(data$combos == "0,0") | (data$combos == "1,1"),]
plot_fn(data_filt, plot_file)

if (FALSE){
# plot all??
    plot_file <- paste(plot_prefix, ".all.pdf", sep="")
    plot_fn(data, plot_file)
}
