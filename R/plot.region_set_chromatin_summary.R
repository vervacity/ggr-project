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
    signals_mat <- t(scale(t(signals)))
    colnames(signals_mat) <- colnames(signals)
    
    # melt
    signals_melt <- melt(signals_mat)
    signals_melt$type <- signal_type
    signals_melt$Var2 <- as.character(signals_melt$Var2)
    signals_melt$day <- gsub("d", "", signals_melt$Var2, fixed=TRUE)
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
plot_file <- paste(plot_prefix, ".TEST.pdf", sep="")
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
    theme_bw() +
    theme(
        text=element_text(size=6))
ggsave(plot_file, height=2, width=3)
