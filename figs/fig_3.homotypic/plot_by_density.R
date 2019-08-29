#!/usr/bin/env Rscript

library(rhdf5)
library(gplots)
library(reshape2)
library(RColorBrewer)

library(grid)
library(gridGraphics)
library(gridExtra)

# functions
grab_grob <- function(fn) {
    grid.grabExpr(grid.echo(fn))
}

plot_heatmap <- function(data) {

    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.07,0.7,0.75)
    mylhei = c(0.10,2,0.25)

    my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(49)

    
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        colsep=0:(ncol(data)+1),
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001,0.001),

        #labRow=labRow,
        #adjRow=adjRow,
        cexRow=0.5,
        #offsetRow=offsetRow,
        
        #labCol=labCol,
        cexCol=0.5,
        #srtCol=45,
        #offsetCol=offsetCol,

        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(2,6,1.9,6),
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

        margins=c(1,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        na.color="white",
        col=my_palette)
        #RowSideColors=rowsidecolors)


}


# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]

# check keys in file
contents <- h5ls(data_file)
dataset_keys <- contents$name

# go through each key
for (i in 1:length(dataset_keys)) {

    # get data
    dataset_key <- dataset_keys[i]
    print(dataset_key)
    data <- h5read(data_file, dataset_key, read.attributes=TRUE)
    pwm_names <- attr(data, "pwm_names")
    pwm_names <- gsub("HCLUST-\\d+.", "", pwm_names)
    pwm_names <- gsub(".UNK.0.A", "", pwm_names)
    num_tasks <- dim(data)[3]
    
    # flatten for clustering
    data_flat <- aperm(data, c(2,3,1))
    data_dims <- dim(data_flat)
    dim(data_flat) <- c(data_dims[1], data_dims[2]*data_dims[3])

    # cluster to get order
    hc <- hclust(dist(data_flat), method="ward.D2")
    ordering <- hc$order
    pwm_names <- pwm_names[ordering]
    
    # clip
    data_melt <- melt(data)
    clip_val <- quantile(data_melt$value, 0.90)
    data[data > clip_val] <- clip_val
    clip_val <- quantile(data_melt$value, 0.20)
    data[data < clip_val] <- clip_val
    data[data == 0] <- NA
    
    # for each task
    grob_list <- list()
    for (task_i in 1:num_tasks) {

        # get data and order
        task_data <- aperm(data)[task_i,,]
        task_data <- task_data[ordering,]
        rownames(task_data) <- pwm_names
        
        # attach pwm names

        # make grob
        heatmap_fn <- function() plot_heatmap(task_data)
        grob_list[[task_i]] <- grab_grob(heatmap_fn)
        
    }

    # plot
    plot_file <- paste(dataset_key, ".counts_v_activity.pdf", sep="")
    pdf(
        plot_file,
        height=3.5,
        width=3,
        onefile=FALSE,
        useDingbats=FALSE,
        family="ArialMT")
    plot.new()
    grid.newpage()
    grid.arrange(grobs=grob_list, nrow=1, ncol=num_tasks,
                 heights=c(3), widths=c(1, 1, 1))
    dev.off()

}


