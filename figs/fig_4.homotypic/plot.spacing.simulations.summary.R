#!/usr/bin/env Rscript

library(rhdf5)
library(gplots)
library(ggplot2)
library(ggridges)
library(reshape2)
library(RColorBrewer)

library(grid)
library(gridGraphics)
library(gridExtra)

# functions
grab_grob <- function(fn) {
    grid.grabExpr(grid.echo(fn))
}

plot_null <- function() { plot.new() }

plot_heatmap <- function(data, labelRow, plot_title) {

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.07,0.7,0.07)
    mylhei = c(0.09,2,0.20)

    # color
    my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(49)

    # adjust row labels
    if (labelRow) {
        labRow <- rownames(data)
    } else {
        labRow <- rep("", nrow(data))
    }
    
    # adjust col labels
    labCol <- rep("", ncol(data))
    for (i in 1:ncol(data)) {
        if ((i-5)%%10 == 0) {
            labCol[i] <- i+5
        }
    }
    
    # heatmap
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        colsep=c(0,(ncol(data)+1)),
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001,0.001),

        labRow=labRow,
        cexRow=0.5,
        offsetRow=-6,

        labCol=labCol,
        cexCol=0.5,
        offsetCol=-0.25,
        srtCol=0,
        
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(2,6,1.4,6),
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

    title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
}


# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
out_prefix <- args[2]
do_filter <- args[3]
dataset_keys <- args[4:length(args)]

plot_titles <- c("day 0", "day 3", "day 6")

# go through each key
for (i in 1:length(dataset_keys)) {

    # get data
    dataset_key <- dataset_keys[i]
    print(dataset_key)
    data <- h5read(
        data_file,
        paste(dataset_key, "/count", sep=""),
        read.attributes=TRUE)
    pwm_names <- attr(data, "pwm_names")
    pwm_names <- gsub("HCLUST-\\d+.", "", pwm_names)
    pwm_names <- gsub(".UNK.0.A", "", pwm_names)
    num_tasks <- dim(data)[3]

    # also get the smaller version
    if (do_filter == "TRUE") {
        keep <- h5read(
            data_file,
            paste(dataset_key, "/nonlinear", sep=""))
        data <- data[,keep==1,]
        pwm_names <- pwm_names[keep==1]
    }

    # remove the ones that never activated
    max_vals <- apply(data, 2, max)
    keep_indices <- which(max_vals > 0.5)
    data <- data[,keep_indices,]
    pwm_names <- pwm_names[keep_indices]
    
    # normalize
    print(dim(data))
    #max_vals <- apply(data, 2, max)
    #data_norm <- sweep(data, 2, max_vals, FUN="/")
    #data <- data_norm
    data[is.na(data)] <- 0
    data[data < 0] <- 0
    
    # flatten for clustering
    data_flat <- aperm(data, c(2,3,1))
    data_dims <- dim(data_flat)
    dim(data_flat) <- c(data_dims[1], data_dims[2]*data_dims[3])
    
    # cluster to get order
    #hc <- hclust(dist(data_flat), method="ward.D2")
    hc <- hclust(dist(data_flat))
    ordering <- hc$order
    pwm_names <- pwm_names[ordering]
    
    # for each task
    grob_list <- list("1"=grab_grob(plot_null))
    for (task_i in 1:num_tasks) {

        # get data and order
        task_data <- aperm(data)[task_i,,]
        task_data <- task_data[ordering,]
        rownames(task_data) <- pwm_names

        # rownorm
        #task_data <- t(apply(task_data, 1, function(x) (x - min(x)) / (max(x) - min(x))))
        #task_data <- data.frame(task_data)
        #rownames(task_data) <- pwm_names
        
        if (task_i == 1) {
            labelRow <- TRUE
        } else {
            labelRow <- FALSE
        }
        
        # make grob
        heatmap_fn <- function() plot_heatmap(task_data, labelRow, plot_titles[task_i])
        grob_list[[task_i+1]] <- grab_grob(heatmap_fn)
        
    }

    # plot
    plot_file <- paste(out_prefix, ".", dataset_key, ".dist_v_activity.pdf", sep="")
    pdf(
        plot_file,
        height=3.5,
        width=2,
        onefile=FALSE,
        useDingbats=FALSE,
        family="ArialMT")
    plot.new()
    grid.newpage()
    grid.arrange(grobs=grob_list, nrow=1, ncol=num_tasks+1,
                 heights=c(3), widths=c(0.25, 0.5, 0.5, 0.5))
    dev.off()

    # redo plot as line plot
    data <- data[,ordering,]
    max_curves <- data.frame()
    for (pwm_idx in 1:length(pwm_names)) {
        pwm_data <- aperm(data)[,pwm_idx,]
        pwm_max <- apply(pwm_data, 1, max)
        pwm_task_idx <- which.max(pwm_max)
        if (length(pwm_task_idx) == 0) {
            pwm_task_idx <- 1
        }
        pwm_data <- pwm_data[pwm_task_idx,]
        max_curves <- rbind(max_curves, pwm_data)
    }

    max_curves <- data.frame(t(max_curves))
    colnames(max_curves) <- pwm_names
    max_vals <- apply(max_curves, 2, max)
    print(length(max_vals))
    factor_pwm_names <- pwm_names[order(max_vals, decreasing=FALSE)]
    
    max_curves$distance <- 1:nrow(max_curves) + 5
    data_melt <- melt(max_curves, id.vars="distance")

    data_melt$variable <- factor(data_melt$variable, levels=factor_pwm_names)
    
    # plot
    plot_file <- paste(out_prefix, ".", dataset_key, ".dist_v_activity.curves.pdf", sep="")
    ggplot(data_melt, aes(x=distance, y=value, group=variable, colour=variable)) +
        geom_line(size=0.115, show.legend=FALSE) +
        geom_hline(size=0.115, linetype="dashed", yintercept=1, colour="black") + 
        labs(y="Accessibility signal log2(FC)") +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            #panel.grid.major.y=element_blank(),
            #panel.grid.major.x=element_line(size=0.115),
            panel.grid.minor=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
            #axis.title.y=element_blank(),
            axis.text.y=element_text(size=5),
            axis.text.x=element_text(size=5),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            legend.background=element_blank(),
            legend.box.background=element_rect(colour="black", fill="white", size=0.115),
            legend.margin=margin(1,1,1,1),
            legend.key=element_blank(),
            legend.key.size=unit(0.1, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=5),
            legend.position=c(0.9, 0.9),
            strip.background=element_blank(),
            strip.text=element_blank()) +
        scale_x_continuous(limits=c(5,28), expand=c(0,0))
    ggsave(plot_file, height=2, width=1.5, useDingbats=FALSE)

    # plot again with ggridges
    plot_file <- paste(out_prefix, ".", dataset_key, ".dist_v_activity.curves.ggridges.pdf", sep="")
    ggplot(data_melt, aes(x=distance, y=variable, group=variable, colour=variable, height=value)) +
        geom_ridgeline(alpha=0.3, size=0.115, show.legend=FALSE, colour="black") +
        labs(y="Accessibility signal log2(FC)") +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.text.y=element_text(size=5),
            axis.text.x=element_text(size=5),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            legend.background=element_blank(),
            legend.box.background=element_rect(colour="black", fill="white", size=0.115),
            legend.margin=margin(1,1,1,1),
            legend.key=element_blank(),
            legend.key.size=unit(0.1, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=5),
            legend.position=c(0.9, 0.9),
            strip.background=element_blank(),
            strip.text=element_blank()) +
        scale_x_continuous(limits=c(5,28), expand=c(0,0))
    ggsave(plot_file, height=2, width=1.5, useDingbats=FALSE)
    
    
}


