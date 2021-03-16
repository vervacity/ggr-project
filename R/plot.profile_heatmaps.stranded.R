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
plot_title <- args[2]
row_sep_file <- args[3]
out_file <- args[4]
rcolorbrewer_palette <- args[5]
start_stop_sets <- args[6:length(args)]

# read in data
data <- read.table(gzfile(deeptools_matrix_file), skip=1, header=FALSE)

# params
mad_factor <- 200
color_granularity <- 50

# figure out row sep
#rowsep <- read.table(row_sep_file, header=FALSE)$V1
rowsep <- c(0, dim(data)[1])

# keep strands to reorient
strands <- data$V6

# cleanup
data$V1 <- NULL
data$V2 <- NULL
data$V3 <- NULL
data$V4 <- NULL
data$V5 <- NULL
data$V6 <- NULL
data$id <- NULL

# reorient negative strand
for (row_i in 1:dim(data)[1]) {
    if (strands[row_i] == "-") {
        # do it per group
        for (j in 1:length(start_stop_sets)) {
            # get coords
            coords <- as.numeric(unlist(strsplit(start_stop_sets[j], ",")))
            start <- coords[1]
            end <- coords[2]

            # adjust accordingly
            data[row_i,start:end] <- rev(data[row_i,start:end])
        }
    }
}

# since using pval, convert to log10(pval)
data <- log10(data)
data[data < 0] <- 0

# breaks: use percentiles to remove extreme outliers
data_melted <- melt(data)
my_breaks <- seq(
    quantile(data_melted$value, 0),
    quantile(data_melted$value, 0.98),
    length.out=color_granularity)

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
plot_profile_heatmap <- function(plot_data, i, plot_key) {

    label <- c("d0", "d3", "d6")

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.02,0.5,0.02)
    mylhei = c(0.1,2,0.27)

    # heatmap
    heatmap.2(
        as.matrix(plot_data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        srtCol=60,
        cexCol=0.5,
        offsetCol=-0.5,
        labCol=c(
            rep("", floor(ncol(plot_data)/2)),
            label[i],
            rep("", ceiling(ncol(plot_data)/2)-1)),
        labRow=rep("", nrow(plot_data)),

        colsep=c(0, ncol(plot_data)),
        rowsep=rowsep,
        sepcolor="black",
        sepwidth=c(0.0001, 0.0001),

        key=plot_key,
        keysize=0.1,
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(0.9,1,1,1) / (par("cex") * 0.66),
            mgp=c(0,-0.1,0),
            tcl=-0.1,
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
        breaks=my_breaks,
        #add.expr=add_borders(i),
        useRaster=TRUE,
        par(xpd=FALSE)) # xpd - where plotting shows up

    if (plot_key) {
        title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
    }
    
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

        if (i == 2) {
            plot_key <- TRUE
        } else {
            plot_key <- FALSE
        }
        
        # plot heatmap and save out
        fn <- function() plot_profile_heatmap(data_subset, i, plot_key)

        # grab grob
        return(grab_grob(fn))
    })

# plot joint heatmap
pdf(file=out_file, height=3, width=0.7, onefile=FALSE, family="ArialMT")
grid.newpage()
grid.arrange(grobs=grob_list, nrow=1, ncol=length(start_stop_sets), clip=FALSE)
dev.off()
