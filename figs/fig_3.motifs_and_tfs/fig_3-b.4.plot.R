#!/usr/bin/env Rscript

# fig 3b
# description: plot importance scores in SNPs and delta output
# data: /srv/scratch/dskim89/ggr/ggr.tronn.2019-06-11.ase

library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
results_file <- args[1]
plot_file <- "fig_3-b.pdf"

# adjust
data <- read.table(results_file, header=TRUE, stringsAsFactors=FALSE)
data_melted <- data
data_melted <- data_melted[data_melted$group != "window_impt",]

# clean up 
data_melted$group[data_melted$group == "no_impt"] <- "- Importance score"
data_melted$group[data_melted$group == "single_impt"] <- "+ Importance score"
data_melted$sig[data_melted$sig == "sig"] <- "significant allele-sensitive ATAC SNP"
data_melted$sig[data_melted$sig == "not_sig"] <- "non-significant SNP"
data_melted$sig <- factor(
    data_melted$sig,
    levels=c(
        "non-significant SNP",
        "significant allele-sensitive ATAC SNP"))

# plot
dodge <- position_dodge(width=0.9)
ggplot(data_melted, aes(x=group, y=delta_logit, colour=sig)) +
    geom_violin(width=0.8, position=dodge, scale="width") +
    geom_boxplot(width=0.2, position=dodge, outlier.shape=NA, show.legend=FALSE) +
    ylim(1, 3.5) +
    labs(y="Model prediction delta (ref - alt)\n") +        
    theme_bw() +
    theme(
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title.x=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.spacing.x=unit(0.5, "cm"))
ggsave(plot_file, width=6, height=3)

