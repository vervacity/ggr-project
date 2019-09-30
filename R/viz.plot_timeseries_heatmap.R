#!/usr/bin/env Rscript

# description: plot already ordered heatmap
# does not adjust the order of input.

library(gplots)
library(RColorBrewer)
library(fastcluster)
library(reshape2)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

set.seed(1337)

# args
args <- commandArgs(trailingOnly=TRUE)
cluster_file <- args[1]
mat_file <- args[2]
out_dir <- args[3]
prefix <- args[4]
plot_title <- args[5]

# read in files
clusters <- read.table(cluster_file, header=TRUE)
if (colnames(clusters)[2] == "gene"){
    colnames(clusters)[2] <- "id"
}
data <- read.table(gzfile(mat_file), header=TRUE)
data$id <- rownames(data)

# merge
data_w_clusters <- merge(clusters, data, by="id", sort=FALSE)
rownames(data_w_clusters) <- data_w_clusters$id
data_w_clusters$id <- NULL

# determine cluster sizes and means (for dendrogram purposes)
cluster_names <- unique(data_w_clusters$cluster)
cluster_sizes <- c()
cluster_means <- data.frame()
for (i in 1:length(cluster_names)) {
    
    single_cluster <- data_w_clusters[data_w_clusters$cluster == cluster_names[i], ]
    single_cluster$cluster <- NULL
    
    single_cluster_z <- t(scale(t(single_cluster), center=TRUE, scale=TRUE))
    
    cluster_mean <- colMeans(single_cluster_z)
    cluster_means <- rbind(cluster_means, cluster_mean)
    cluster_sizes <- c(cluster_sizes, nrow(single_cluster))
}

# determine row sep points and save out
rowsep <- c(0)
cluster <- data_w_clusters$cluster[1]
cluster_ids_per_example <- data_w_clusters$cluster
data_w_clusters$cluster <- NULL
color_bar_clusters <- c()
color_bar_cluster <- 1
for (i in 1:nrow(data_w_clusters)) {
    
    if (cluster_ids_per_example[i] != cluster) {
        rowsep <- c(rowsep, i)
        cluster <- cluster_ids_per_example[i]
        color_bar_cluster <- color_bar_cluster + 1
    }
    color_bar_clusters <- c(color_bar_clusters, color_bar_cluster)
}
rowsep <- c(rowsep, nrow(data_w_clusters))
row_sep_file <- paste(out_dir, "/", prefix, ".row_seps.txt", sep="")
write.table(rowsep, row_sep_file, row.names=FALSE, col.names=FALSE)

# clean up, just in case
data_w_clusters$histone_cluster <- NULL
data_w_clusters$H3K27ac <- NULL
data_w_clusters$H3K4me1 <- NULL
data_w_clusters$H3K27me3 <- NULL
data_w_clusters$fake_cluster <- NULL

# TODO here's where to ajust the color scaling
if (TRUE) {
    # already rlog space, so subtract to get FC
    data_z <- data_w_clusters - data_w_clusters[,1]
} else {
    data_z <- t(scale(t(data_w_clusters), center=TRUE, scale=TRUE))
}

# percentile clip
if (grepl("epigenome", prefix)) {
    thresholds <- quantile(melt(data_z)$value, c(0.01, 0.99))
} else {
    thresholds <- quantile(melt(data_z)$value, c(0.10, 0.90))
}
data_z[data_z < thresholds[1]] <- thresholds[1]
data_z[data_z > thresholds[2]] <- thresholds[2]

# subsample for plotting
max_n <- 2000
if (nrow(data_w_clusters) > max_n) {
    print("WARNING: num examples too high, not plotting")
    quit()
}

# plotting
plot_file <- paste(out_dir, "/", prefix, ".clusters.heatmap.pdf", sep="")
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

# color bar
cluster_palette <- get_trajectory_palette(nrow(cluster_means))
cluster_colors <- cluster_palette[color_bar_clusters]

# heatmap2 grid
mylmat = rbind(c(0,0,3,0),c(4,1,2,0),c(0,0,5,0))
mylwid = c(0.05,0.1,0.6,0.4)
mylhei = c(0.1,2,0.27)

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
    rowsep=rowsep,
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
