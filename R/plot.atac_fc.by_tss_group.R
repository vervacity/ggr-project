#!/usr/bin/env Rscript

library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read data
data <- read.table(data_file, sep="\t", header=TRUE)

# clean up labels
data$group <- gsub("dynamic", "Dynamic", data$group)
data$group <- gsub("stable", "Stable", data$group)
data$group <- gsub("non_expr", "Non-expressed", data$group)
data$group <- gsub(" (n", "\n(n", data$group, fixed=TRUE)

# reorder
categories <- levels(as.factor(data$group))
ordered_levels <- c()
for (i in 1:length(categories)) {
    if (grepl("ynamic", categories[i])) {
        ordered_levels <- c(ordered_levels, categories[i])
    }
}
for (i in 1:length(categories)) {
    if (grepl("table", categories[i])) {
        ordered_levels <- c(ordered_levels, categories[i])
    }
}
for (i in 1:length(categories)) {
    if (grepl("on-expr", categories[i])) {
        ordered_levels <- c(ordered_levels, categories[i])
    }
}
for (i in 1:length(categories)) {
    if (grepl("Distal", categories[i])) {
        ordered_levels <- c(ordered_levels, categories[i])
    }
}
data$group <- factor(data$group, levels=ordered_levels)



# plot
ggplot(data, aes(x=group, y=foldchange)) +
    geom_violin(size=0.230, aes(colour=group), show.legend=FALSE) +
    geom_boxplot(size=0.230, aes(colour=group), width=0.25, outlier.shape=NA, show.legend=FALSE) + 
    labs(
        x="Region category",
        y="Max log2 fold change\nover differentiation",
        title="Comparison of accessibility\nvariability at TSSs") +
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
        #axis.title.x=element_text(vjust=1.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(vjust=-1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, angle=30, vjust=1, hjust=1),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"))

ggsave(plot_file, height=1.75, width=2, useDingbats=FALSE)

