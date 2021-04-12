#!/usr/bin/env Rscript


library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
summary_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(summary_file, sep="\t", header=TRUE)

# get an ordering for motifs and terms
data_unmelt <- dcast(data, formula=motif~term.name, value.var="X.log10pval")
data_unmelt[is.na(data_unmelt)] <- 0
motifs <- data_unmelt$motif
data_unmelt$motif <- NULL


my_hc <- hclust(dist(data_unmelt), method="ward.D")
#my_hc <- hclust(dist(data_unmelt))
print(length(my_hc$order))
print(length(motifs))
motif_ordering <- motifs[my_hc$order]

#my_hc <- hclust(dist(t(data_unmelt)), method="ward.D")
my_hc <- hclust(dist(t(data_unmelt)))
print(length(my_hc$order))
print(length(colnames(data_unmelt)))
term_ordering <- colnames(data_unmelt)[my_hc$order]


data$motif <- factor(data$motif, levels=motif_ordering)
data$term.name <- factor(data$term.name, levels=term_ordering)

# plot
p_go <- ggplot(data, aes(x=term.name, y=motif)) +
    geom_point(shape=20, aes(size=X.log10pval)) +
    labs(size="-log10(p)", y="Motif", x=NULL) +
    theme_bw() +
    theme(
        text=element_text("ArialMT"),
        plot.margin=margin(t=4, r=0, l=1),
        panel.background=element_blank(),
        panel.border=element_rect(size=0.115),
        panel.grid.major=element_line(size=0.115),
        axis.text.x=element_text(size=6, angle=60, hjust=1),
        axis.text.y=element_text(size=6, hjust=1),
        #axis.text.y=element_blank(),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        legend.spacing.x=unit(0.05, "in")) +
    scale_size_continuous(range=c(0.25,3)) +
    scale_y_discrete(drop=FALSE)

ggsave(plot_file, height=5.5, width=6)
