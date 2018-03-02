#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(gplots)
library(RColorBrewer)
library(reshape2)

library(grid)
library(gridGraphics)
library(gridExtra)

# args
args <- commandArgs(trailingOnly=TRUE)
deeptools_matrix_file <- args[1]
out_file <- args[2]
rcolorbrewer_palette <- args[3]
start_stop_sets <- args[4:length(args)]

# read in data
data <- read.table(gzfile(deeptools_matrix_file), skip=1, header=FALSE)

# params
mad_factor <- 200
color_granularity <- 50

# cleanup
data$V1 <- NULL
data$V2 <- NULL
data$V3 <- NULL
data$V4 <- NULL
data$V5 <- NULL
data$V6 <- NULL

# breaks: use median absolute deviation (MAD), do this across all data
# minimum threshold at 2, which is pval 0.01
data_melted <- melt(data)
#signal_max_thresh <- median(data_melted$value) + mad_factor * mad(data_melted$value)
#my_breaks <- seq(2, signal_max_thresh, length.out=color_granularity)
my_breaks <- seq(
    quantile(data_melted$value, 0.01),
    quantile(data_melted$value, 0.98),
    length.out=color_granularity)

# consider using distribution to space out breaks?
#data_melted_filt <- data_melted[data_melted$value <= signal_max_thresh,]
#my_breaks <- quantile(data_melted_filt$value, seq(0,1, length.out=color_granularity))
#my_breaks <- as.numeric(my_breaks)

# colors
my_palette <- colorRampPalette(brewer.pal(9, rcolorbrewer_palette))(color_granularity-1)
my_palette <- colorRampPalette(c("white", brewer.pal(9, rcolorbrewer_palette)))(color_granularity-1)

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
    
    # set up sizing
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.25,1,0.25)
    mylhei = c(0.125,2,0.25)

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
        key.par=list(pin=c(4,0.1),
            mar=c(4.1,0,0.1,0),
            mgp=c(3,2,0),
            cex.axis=1.0,
            font.axis=1),
        margins=c(2,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        col=my_palette,
        breaks=my_breaks,
        add.expr=add_borders(i),
        useRaster=FALSE,
        par(xpd=FALSE))

}

# grab grob fn
grab_grob <- function(fn) {
    grid.echo(fn)
    grid.grab()
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
pdf(file=out_file, height=9, width=4, onefile=FALSE, family="ArialMT")
grid.newpage()
grid.arrange(grobs=grob_list, nrow=1, ncol=length(start_stop_sets), clip=TRUE)
dev.off()
