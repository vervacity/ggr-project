#!/usr/bin/env Rscript

# Description: given a histone overlap file
# return color bars to mark splits on dynamic data

library(gplots)
library(RColorBrewer)

set.seed(1337)

args <- commandArgs(trailingOnly=TRUE)
overlap_file <- args[1]
out_plot_prefix <- args[2]

overlap_data <- read.table(gzfile(overlap_file), header=TRUE)

rownames(overlap_data) <- overlap_data$region
overlap_data$region <- NULL

# subsample
evenly_spaced_subsample <- overlap_data[seq(1, nrow(overlap_data), 20),]

# plotting
plot_heatmap <- function(data, color_bar, out_plot_file) {
    my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

    mylmat = rbind(c(0,0,3,0),c(4,1,2,0),c(0,0,5,0)) # only color bar
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
        sepcolor="black",
        RowSideColors=color_bar)
    dev.off()
}

# H3K27ac color bar
H3K27ac_clusters <- evenly_spaced_subsample$H3K27ac + 1
H3K27ac_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(10-1)
H3K27ac_palette <- c(H3K27ac_palette, "#ffffff")
H3K27ac_colors <- H3K27ac_palette[H3K27ac_clusters]
H3K27ac_file <- paste(out_plot_prefix, "colorbar_H3K27ac.heatmap.png", sep=".")

plot_heatmap(evenly_spaced_subsample, H3K27ac_colors, H3K27ac_file)

# H3K4me1 color bar
H3K4me1_clusters <- evenly_spaced_subsample$H3K4me1 + 1
H3K4me1_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(10-1)
H3K4me1_palette <- c(H3K4me1_palette, "#ffffff")
H3K4me1_colors <- H3K4me1_palette[H3K4me1_clusters]
H3K4me1_file <- paste(out_plot_prefix, "colorbar_H3K4me1.heatmap.png", sep=".")

plot_heatmap(evenly_spaced_subsample, H3K4me1_colors, H3K4me1_file)

# H3K27me3 color bar
H3K27me3_clusters <- evenly_spaced_subsample$H3K27me3 + 1
H3K27me3_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(10-1)
H3K27me3_palette <- c(H3K27me3_palette, "#ffffff")
H3K27me3_colors <- H3K27me3_palette[H3K27me3_clusters]
H3K27me3_file <- paste(out_plot_prefix, "colorbar_H3K27me3.heatmap.png", sep=".")

plot_heatmap(evenly_spaced_subsample, H3K27me3_colors, H3K27me3_file)

