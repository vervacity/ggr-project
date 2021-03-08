#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
plot_prefix <- args[1]
signal_files <- args[2:length(args)]

# load in data and normalize
for (signal_i in 1:length(signal_files)) {
    signal_file <- signal_files[signal_i]
    signals <- read.table(signal_file, header=TRUE, stringsAsFactors=FALSE)

    if (dim(signals)[1] == 0) {
        next
    }
    
    # figure out what kind of data this is
    if (grepl("atac.ends", signal_file)) {
        signal_type <- "ATAC-seq"
    }
    if (grepl("H3K27ac.midpoint", signal_file)) {
        signal_type <- "H3K27ac ChIP-seq"
    }
    if (grepl("H3K4me1.midpoint", signal_file)) {
        signal_type <- "H3K4me1 ChIP-seq"
    }
    
    # zscore
    #signals_mat <- t(scale(t(signals)))
    #colnames(signals_mat) <- colnames(signals)
    
    # melt
    #signals_melt <- melt(signals_mat)
    signals_melt <- melt(signals)
    signals_melt$type <- signal_type
    signals_melt$variable <- as.character(signals_melt$variable)
    signals_melt$day <- gsub("d", "", signals_melt$variable, fixed=TRUE)
    if (signal_type == "ATAC-seq") {
        signals_melt$day <- gsub("0$", ".0", signals_melt$day)
        signals_melt$day <- gsub("5$", ".5", signals_melt$day)
    }
    
    # merge together
    if (signal_i == 1) {
        all_signals <- signals_melt
    } else {
        all_signals <- rbind(all_signals, signals_melt, stringsAsFactors=FALSE)
    }

}
# factorize
all_signals$day <- as.numeric(as.character(all_signals$day))
#all_signals$day <- as.factor(all_signals$day)
all_signals$type <- as.factor(all_signals$type)

# plot
plot_file <- paste(plot_prefix, ".avg_dynamics.pdf", sep="")
#ggplot(all_signals, aes(x=day, y=value)) + geom_boxplot(aes(colour=type)) +
ggplot(all_signals, aes(x=day, y=value, colour=type)) +
    stat_summary(
        fun=mean,
        fun.min=function(x) mean(x) - sd(x),
        fun.max=function(x) mean(x) + sd(x),
        geom="pointrange") +
    stat_summary(
        fun=mean,
        geom="line") +
    labs(x="Day", y="Signal log(FC)") +
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
