#!/usr/bin/env Rscript

library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
rna_mat_file <- args[1]
gene_id <- args[2]
plot_file <- args[3]

# read in data and get gene vector
rna <- read.table(gzfile(rna_mat_file), sep="\t")

# get gene and adjust
vals <- data.frame(t(rna[gene_id,]))
colnames(vals) <- "rlog_val"
vals$timepoint <- gsub("d", "", rownames(vals))
vals$timepoint <- gsub("0$", ".0", vals$timepoint)
vals$timepoint <- gsub("5$", ".5", vals$timepoint)

# plot
ggplot(vals, aes(x=timepoint, y=1, fill=rlog_val)) +
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
        axis.text=element_blank()) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_distiller(palette="Purples", direction=1)
ggsave(plot_file, height=0.1, width=0.5, useDingbats=FALSE)
