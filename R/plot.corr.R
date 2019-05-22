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
print(head(data))


# grid
mylmat = rbind(
    c(0,0,5,0),
    c(0,0,2,0),
    c(4,1,3,0),
    c(0,0,6,0))
mylwid = c(0.05,0.25,4.4,0.05)
mylhei = c(0.05,0.25,6,0.5)

# pull the trajectory palettes
row_colors <- get_trajectory_palette(15)
col_colors <- get_trajectory_palette(11)
my_palette <- rev(colorRampPalette(brewer.pal(9, "PuOr"))(49))

# breaks
my_breaks <- seq(-1, 1, length.out=50)

# plot
pdf(plot_file, height=9, width=9)

heatmap.2(
    as.matrix(data),
    Rowv=FALSE,
    Colv=FALSE,
    dendrogram="none",
    trace="none",
    density.info="none",
    #labCol=c(
    #    rep("", floor(ncol(plot_data)/2)),
    #    label[i],
    #    rep("", ceiling(ncol(plot_data)/2)-1)),
    #labRow=rep("", nrow(plot_data)),
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
    rowsep=1:nrow(data),
    colsep=1:ncol(data),
    RowSideColors=row_colors,
    ColSideColors=col_colors,
    sepcolor="black")
dev.off()
