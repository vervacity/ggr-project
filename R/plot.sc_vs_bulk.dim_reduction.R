#!/usr/bin/env Rscript

library(umap)
library(ggplot2)
library(RColorBrewer)
library(colorspace)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]

# read data
data <- read.table(data_file, sep="\t", row.names=1, header=TRUE)

# normalize by column sums, then scale
data <- t(t(data)/colSums(data))
data <- scale(data, center=TRUE, scale=TRUE)

if (FALSE) {
    # TODO consider filtering for most varying features
    data_sd <- apply(data[!grepl("sc", rownames(data)),], 1, sd)
    #data_sd <- apply(data, 1, sd)
    thresh <- quantile(data_sd, prob=0.10)
    data <- data[data_sd > thresh,]
}

data <- t(data) # want each row to be a sample

# split sc and bulk
sc_data <- data[grepl("sc", rownames(data)),]
bulk_data <- data[!grepl("sc", rownames(data)),]

# 28, 33
seeds <- c(18, 27, 28, 33, 49)
seeds <- c(28, 33)
rand_i <- 28

# then umap
# rand seeds: 28, 27, 18, 5
sc.umap <- umap(
    data, n_neighbors=15, verbose=TRUE, random_state=rand_i) #42
#bulk_layout <- predict(sc.umap, bulk_data)

# clean up for plotting
layout <- sc.umap$layout
#layout <- layout[grepl("sc", rownames(layout)),]
#layout <- rbind(layout, bulk_layout)
layout <- data.frame(layout)

colnames(layout) <- c("x", "y")
group <- do.call(rbind, strsplit(rownames(layout), "_"))
layout$group <- group[,1]
layout$rep <- group[,2]
layout <- layout[layout$group != "d35",]
layout <- layout[layout$group != "d40",]
layout <- layout[layout$group != "d55",]

# set up a color column to properly color points
layout$rep_label <- layout$rep
layout$rep_label[grepl("sc", layout$group)] <- "sc"
layout$color <- paste(layout$rep_label, layout$group, sep="_")

# set up colors
ggr_colors <- get_ggr_timepoint_colors()
my_colors_r1 <- lighten(ggr_colors, amount=0.2)
my_colors_r2 <- darken(ggr_colors, amount=0.2)
sc_colors <- lighten(
    c(ggr_colors[1], ggr_colors[7], ggr_colors[10]), amount=0.5) #0.2
#ggr_colors <- c(ggr_colors, sc_colors)
ggr_colors <- c(my_colors_r1, my_colors_r2, sc_colors)

# split out to highlight bulk
sc_layout <- layout[grepl("sc", layout$group),]
bulk_layout <- layout[!grepl("sc", layout$group),]

# plot
ggplot(layout, aes(x=x, y=y, group=color)) +
    geom_point(size=0.115, aes(colour=color)) +
    geom_point(data=bulk_layout, shape=21, size=1, aes(x=x, y=y)) +
    labs(
        x="UMAP 1", y="UMAP 2", title="Bulk ATAC-seq projected\nto scATAC-seq") +
    theme_bw() +
    theme(
        aspect.ratio=1,
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
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
        axis.ticks.length=unit(0.01, "in"),
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in")) +
    scale_color_manual(values=ggr_colors) +
    guides(color=guide_legend(nrow=10))
        
plot_file <- paste(plot_prefix, ".", rand_i, ".umap.pdf", sep="")
#ggsave(plot_file, height=1.5, width=2, useDingbats=FALSE)
ggsave(plot_file, height=2, width=4, useDingbats=FALSE)
