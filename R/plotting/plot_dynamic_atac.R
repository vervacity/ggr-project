#!/usr/bin/env Rscript


# Description: given soft ATAC clusters and hard histone clusters
# sort by those clusterings and output the ordered list and
# heatmap


library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(fastcluster)

set.seed(1337)

args <- commandArgs(trailingOnly=TRUE)
atac_timepoints_w_clusters_file <- args[1]
histone_clusters_file <- args[2]
out_plot_prefix <- args[3]
out_subsample_file <- args[4]

atac_timepoints_w_clusters <- read.table(
    gzfile(atac_timepoints_w_clusters_file),
    header=FALSE,
    row.names=1,
    stringsAsFactors=FALSE)
histone_clusters <- read.table(gzfile(histone_clusters_file), header=TRUE)

# adjust files
colnames(atac_timepoints_w_clusters)[ncol(atac_timepoints_w_clusters)] <- "ATAC"
atac_timepoints_w_clusters$region <- rownames(atac_timepoints_w_clusters)
#histone_clusters <- histone_clusters[, c("region", "histone_cluster")] # TODO use histone marks separately

# merge 
epigenome_data <- merge(atac_timepoints_w_clusters, histone_clusters, by="region")
rownames(epigenome_data) <- epigenome_data$region
epigenome_data$region <- NULL

# and sort with atac cluster followed by histone cluster
epigenome_data_sorted <- epigenome_data[with(epigenome_data, order(ATAC, H3K27ac, H3K4me1)), ]

# save out the subsample to be used for deeptools
evenly_spaced_subsample <- epigenome_data_sorted[seq(1, nrow(epigenome_data_sorted), 20), ]
write.table(rownames(evenly_spaced_subsample),
            file=gzfile(out_subsample_file),
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# and now plot out
atac_clusters <- evenly_spaced_subsample$ATAC
H3K27ac_clusters <- evenly_spaced_subsample$H3K27ac
H3K4me1_clusters <- evenly_spaced_subsample$H3K4me1

evenly_spaced_subsample$ATAC <- NULL
evenly_spaced_subsample$histone_cluster <- NULL
evenly_spaced_subsample$H3K27ac <- NULL
evenly_spaced_subsample$H3K4me1 <- NULL

# rowsep (based on ATAC)
row_sep <- c()
cluster <- 1
for (i in 1:nrow(evenly_spaced_subsample)) {
    if (atac_clusters[i] != cluster) {
        row_sep <- c(row_sep, i)
        cluster <- atac_clusters[i]
    }
}

# plotting
plot_heatmap <- function(data, color_bar, out_plot_file, rowsep) {
    my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

    #mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    #mylwid = c(2,6,2)
    #mylhei = c(0.5,12,1.5)

    mylmat = rbind(c(0,0,3,0),c(4,1,2,0),c(0,0,5,0))
    mylwid = c(2,0.5,6,2)
    mylhei = c(0.5,12,1.5)

    png(out_plot_file, height=18, width=6, units="in", res=200)
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace='none',
        density.info="none",
        keysize=0.1,
        key.title=NA,
        key.xlab=NA,
        key.par=list(pin=c(4,0.1),
            mar=c(6.1,0,5.1,0),
            mgp=c(3,2,0),
            cex.axis=2.0,
            font.axis=2),
        srtCol=45,
        cexCol=3.0,
        labRow="",
        margins=c(3,0),
        col=my_palette,
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        rowsep=rowsep,
        sepcolor="black",
        RowSideColors=color_bar)
    dev.off()
}


# ATAC color bar
atac_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(atac_clusters)))
atac_colors <- atac_palette[atac_clusters]
atac_file <- paste(out_plot_prefix, "colorbar_ATAC.heatmap.png", sep="")

plot_heatmap(evenly_spaced_subsample, atac_colors, atac_file, row_sep)


# H3K27ac color bar
H3K27ac_clusters <- H3K27ac_clusters + 1
H3K27ac_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(10-1)
H3K27ac_palette <- c(H3K27ac_palette, "#ffffff")
H3K27ac_colors <- H3K27ac_palette[H3K27ac_clusters]
H3K27ac_file <- paste(out_plot_prefix, "colorbar_H3K27ac.heatmap.png", sep="")

plot_heatmap(evenly_spaced_subsample, H3K27ac_colors, H3K27ac_file, row_sep)


# H3K4me1 color bar
H3K4me1_clusters <- H3K4me1_clusters + 1
H3K4me1_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(10-1)
H3K4me1_palette <- c(H3K4me1_palette, "#ffffff")
H3K4me1_colors <- H3K4me1_palette[H3K4me1_clusters]
H3K4me1_file <- paste(out_plot_prefix, "colorbar_H3K4me1.heatmap.png", sep="")

plot_heatmap(evenly_spaced_subsample, H3K4me1_colors, H3K4me1_file, row_sep)
