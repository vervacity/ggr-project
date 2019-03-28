#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(gplots)
library(RColorBrewer)
library(reshape2)

library(grid)
library(gridGraphics)
library(gridExtra)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
deeptools_matrix_file <- args[1]
row_sep_file <- args[2]
out_file <- args[3]
rcolorbrewer_palette <- args[4]
start_stop_sets <- args[5:length(args)]

# read in data
data <- read.table(gzfile(deeptools_matrix_file), skip=1, header=FALSE)

# params
mad_factor <- 200
color_granularity <- 50

# figure out row sep
rowsep <- read.table(row_sep_file, header=FALSE)$V1

# cleanup
data$V1 <- NULL
data$V2 <- NULL
data$V3 <- NULL
data$V4 <- NULL
data$V5 <- NULL
data$V6 <- NULL
data$id <- NULL

# since using pval, convert to log10(pval)
data <- log10(data)
data[data < 0] <- 0

# breaks: use percentiles to remove extreme outliers
data_melted <- melt(data)
my_breaks <- seq(
    quantile(data_melted$value, 0),
    quantile(data_melted$value, 0.98),
    length.out=color_granularity)

# consider using distribution to space out breaks?
#data_melted_filt <- data_melted[data_melted$value <= signal_max_thresh,]
#my_breaks <- quantile(data_melted_filt$value, seq(0,1, length.out=color_granularity))
#my_breaks <- as.numeric(my_breaks)

# colors
my_palette <- get_ggr_assay_palette(rcolorbrewer_palette, color_granularity)

# set up horizontal borders
h_lines <- c(0, nrow(data)+1)
v_lines <- list()
for (i in 1:length(start_stop_sets)) {
    # get coords
    coords <- as.numeric(unlist(strsplit(start_stop_sets[i], ",")))
    start <- coords[1]
    end <- coords[2]

    # add
    v_lines[[i]] <- c(0, end - start + 2)
    
}

# pretty up heatmap2
add_borders <- function(i) {
    abline(
        h=h_lines,
        v=v_lines[[i]])
}

# heatmap fn
plot_profile_heatmap <- function(plot_data, i) {

    label <- c("d0", "d3", "d6")

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.05,1,0.05)
    mylhei = c(0.25,4,0.5)

    # heatmap
    heatmap.2(
        as.matrix(plot_data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        labCol=c(
            rep("", floor(ncol(plot_data)/2)),
            label[i],
            rep("", ceiling(ncol(plot_data)/2)-1)),
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
        sepcolor="black",
        add.expr=add_borders(i),
        useRaster=TRUE,
        par(xpd=FALSE)) # xpd - where plotting shows up

}

# grab grob fn
grab_grob <- function(fn) {
    #grid.echo(fn)
    grid.grabExpr(grid.echo(fn))
    #grid.grab()
}

# get grobs
grob_list <- lapply(
    1:length(start_stop_sets),
    function(i) {
        print(i)
        
        # get coords
        coords <- as.numeric(unlist(strsplit(start_stop_sets[i], ",")))
        start <- coords[1]
        end <- coords[2]

        # get data subset
        data_subset <- data[,start:end]
        
        # plot heatmap and save out
        fn <- function() plot_profile_heatmap(data_subset, i)

        # grab grob
        return(grab_grob(fn))
    })

# plot joint heatmap
pdf(file=out_file, height=7, width=2, onefile=FALSE, family="ArialMT")
grid.newpage()
grid.arrange(grobs=grob_list, nrow=1, ncol=length(start_stop_sets), clip=TRUE)
dev.off()
