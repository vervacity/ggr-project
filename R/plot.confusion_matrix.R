#!/usr/bin/env Rscript

# description: plot chrom states
library(gplots)
library(reshape2)

# load style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
confusion_mat_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(gzfile(confusion_mat_file), sep="\t", header=TRUE, row.names=1)
data[is.na(data)] <- 0

# thresh
top_thresh <- quantile(melt(data)$value, 0.95)
data[data > top_thresh] <- top_thresh

# heatmap fn
plot_heatmap <- function(plot_data, color_set) {

    # grid
    mylmat = rbind(
        c(0,0,5,0),
        c(0,0,2,0),
        c(4,1,3,0),
        c(0,0,6,0))
    mylwid = c(0.05,0.15,1,0.05)
    mylhei = c(0.05,0.15,1.2,0.5)

    # get GGR colors
    granularity <- 49
    my_palette <- colorRampPalette(
        brewer.pal(9, "YlOrBr"))(granularity)
    row_colors <- get_trajectory_palette(15)
    col_colors <- get_trajectory_palette(11)
    
    # sep
    rowsep <- 0:nrow(plot_data)
    colsep <- 0:ncol(plot_data)

    # adjust color range
    my_breaks <- seq(0.0, max(plot_data), length.out=granularity+1)
    
    # heatmap
    heatmap.2(
        as.matrix(plot_data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        labCol=rep("", ncol(plot_data)),        
        labRow=rep("", nrow(plot_data)),

        keysize=0.1,
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(2,1,0.5,1),
            mgp=c(0,-0.1,0),
            cex.axis=0.6,
            tcl=-0.1,
            lend=2,
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
        breaks=my_breaks,
        rowsep=rowsep,
        colsep=colsep,
        sepcolor="black",
        sepwidth=c(0.0001, 0.0001),
        
        RowSideColors=row_colors,
        ColSideColors=col_colors)

}


# plot joint heatmap
# NOTE: don't know how to shrinking heatmap2 smaller than this
pdf(file=plot_file, height=2, width=1.25, useDingbats=FALSE, onefile=FALSE, family="ArialMT")
plot_heatmap(data, "Blues")
dev.off()
