#!/usr/bin/env Rscript

# description: plot simple correlation matrix
library(gplots)
library(RColorBrewer)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
mat_file <- args[1]
plot_file <- args[2]

# load data
data <- read.table(gzfile(mat_file), sep="\t", header=TRUE, row.names=1)

# grid
mylmat = rbind(
    c(0,0,5,0),
    c(0,0,2,0),
    c(4,1,3,0),
    c(0,0,6,0))
mylwid = c(0.05,0.15,1,0.05)
mylhei = c(0.05,0.15,1.2,0.5)

# pull the trajectory palettes
row_colors <- get_trajectory_palette(15)
col_colors <- get_trajectory_palette(11)
my_palette <- rev(colorRampPalette(brewer.pal(9, "PuOr"))(49))

# breaks
my_breaks <- seq(-1, 1, length.out=50)

# plot
pdf(plot_file, height=2, width=1.25, useDingbats=FALSE, family="ArialMT")

heatmap.2(
    as.matrix(data),
    Rowv=FALSE,
    Colv=FALSE,
    dendrogram="none",
    trace="none",
    density.info="none",
    
    labCol=rep("", ncol(data)),
    labRow=rep("", ncol(data)),

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
    rowsep=0:nrow(data),
    colsep=0:ncol(data),
    sepcolor="black",
    sepwidth=c(0.0001, 0.0001),
    
    RowSideColors=row_colors,
    ColSideColors=col_colors)
dev.off()
