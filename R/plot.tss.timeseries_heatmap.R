#!/usr/bin/env Rscript

# description: plot already ordered heatmap
# does not adjust the order of input.

library(gplots)
library(RColorBrewer)
library(reshape2)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

set.seed(1337)

# args
args <- commandArgs(trailingOnly=TRUE)
mat_file <- args[1]
row_seps_file <- args[2]
plot_prefix <- args[3]
plot_title <- args[4]

# read in data
data <- read.table(mat_file, sep="\t", header=TRUE)
#row_seps <- read.table(row_seps_file, header=FALSE)$V1
row_seps <- c(0, dim(data)[1])

# subset to just timepoints
signal_headers <- colnames(data)[grepl("^d", colnames(data))]
data_vals <- data[,signal_headers]


# already rlog space, so subtract to get FC
    data_z <- data_vals - data_vals[,1]

# percentile clip
thresholds <- quantile(melt(data_z)$value, c(0.01, 0.99))
#if (grepl("epigenome", plot_prefix)) {
#    thresholds <- quantile(melt(data_z)$value, c(0.01, 0.99))
#} else {
#    thresholds <- quantile(melt(data_z)$value, c(0.10, 0.90))
#}
data_z[data_z < thresholds[1]] <- thresholds[1]
data_z[data_z > thresholds[2]] <- thresholds[2]

# subsample for plotting
max_n <- 2000
if (nrow(data_z) > max_n) {
    print("WARNING: num examples too high, not plotting")
    quit()
}

# plotting
plot_file <- paste(plot_prefix, ".heatmap.pdf", sep="")
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

# color bar
cluster_palette <- get_trajectory_palette(11) # TODO set as needed

color_bar_clusters <- c()
color_bar_cluster <- 1
row_sep_i <- 2
for (row_i in 1:nrow(data_z)) {
    if (row_i == row_seps[row_sep_i]) {
        color_bar_cluster <- color_bar_cluster + 1
        row_sep_i <- row_sep_i + 1
    }
    color_bar_clusters <- c(color_bar_clusters, color_bar_cluster)
}



cluster_colors <- cluster_palette[color_bar_clusters]

# heatmap2 grid
mylmat = rbind(c(0,0,3,0),c(4,1,2,0),c(0,0,5,0))
mylwid = c(0.05,0.1,0.6,0.4)
mylhei = c(0.1,2,0.27)

#pdf(plot_file, height=3, width=1.25, family="ArialMT")
pdf(plot_file, height=3, width=1.25, family="ArialMT")
heatmap.2(
    as.matrix(data_z),
    Rowv=FALSE,
    Colv=FALSE,
    dendrogram="none",
    trace='none',
    density.info="none",
    
    srtCol=60,
    cexCol=0.5,
    offsetCol=-0.5,
    labRow="",
    rowsep=row_seps,
    colsep=c(0, ncol(data_z)),
    sepcolor="black",
    sepwidth=c(0.0001, 0.0001),
    
    keysize=0.1,
    key.title=NA,
    key.xlab=NA,
    key.par=list(
        mar=c(0.9,1,1,1),
        mgp=c(0,-0.1,0),
        tcl=-0,1,
        lend=2,
        cex.axis=0.6,
        bty="n"),
    key.xtickfun=function() {
        breaks <- pretty(parent.frame()$breaks)
        breaks <- breaks[c(1,length(breaks))]
        list(at = parent.frame()$scale01(breaks),
             labels = breaks)},

    margins=c(0,0),
    lmat=mylmat,
    lwid=mylwid,
    lhei=mylhei,
    
    col=my_palette,
    RowSideColors=cluster_colors)
title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
dev.off()
