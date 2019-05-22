#!/usr/bin/env Rscript

# description: plot chrom states
library(gplots)
#library(RColorBrewer)
#library(viridis)
#library(ggsci)

# load style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
confusion_mat_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(gzfile(confusion_mat_file), sep="\t", header=TRUE, row.names=1)
print(head(data))

# heatmap fn
plot_heatmap <- function(plot_data, color_set) {

    # grid
    mylmat = rbind(
        c(0,0,5,0),
        c(0,0,2,0),
        c(4,1,3,0),
        c(0,0,6,0))
    mylwid = c(0.05,0.25,4.4,0.05)
    mylhei = c(0.05,0.25,6,0.5)

    # get GGR colors
    granularity <- 49
    my_palette <- colorRampPalette(
        brewer.pal(9, "YlGnBu"))(granularity)
    row_colors <- get_trajectory_palette(15)
    col_colors <- get_trajectory_palette(11)
    
    # sep
    rowsep <- 1:nrow(plot_data)
    colsep <- 1:ncol(plot_data)

    # adjust color range
    my_breaks <- seq(0.0, max(plot_data), length.out=granularity+1)
    #my_breaks <- seq(min(plot_data), max(plot_data), length.out=granularity+1)
    
    # heatmap
    heatmap.2(
        as.matrix(plot_data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        #labCol=c(
        #    rep("", floor(ncol(plot_data)/2)),
        #    label[i],
        #    rep("", ceiling(ncol(plot_data)/2)-1)),
        labRow=rep("", nrow(plot_data)),
        keysize=0.1,
        key.title=NA,
        key.xlab=NA,
        key.par=list(pin=c(4,0.1),
            mar=c(2.1,0,2.1,0),
            mgp=c(3,1,0),
            cex.axis=1.0,
            font.axis=2),
        srtCol=45,
        cexCol=1.25,
        margins=c(1,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        col=my_palette,
        breaks=my_breaks,
        rowsep=rowsep,
        colsep=colsep,
        RowSideColors=row_colors,
        ColSideColors=col_colors,
        sepcolor="black")

}


# plot joint heatmap
pdf(file=plot_file, height=9, width=9, onefile=FALSE, family="ArialMT")
plot_heatmap(data, "Blues")
dev.off()
