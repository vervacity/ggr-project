#!/usr/bin/env Rscript

library(rhdf5)
library(reshape2)
library(ggplot2)

data_file <- "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05/motifs.input_x_grad.background/ggr.scanmotifs.h5"

# densities
densities <- h5read(data_file, "sequence.active.pwm-hits.densities.max")[,1,]
activity <- h5read(data_file, "ATAC_SIGNALS.NORM")
pool_window <- 20

# pwm names
features <- h5read(data_file, "sequence-weighted.active.pwm-scores.thresh.max.idx", read.attributes=TRUE)
pwm_names <- attr(features, "pwm_names")

# read in motifs of interest?
sig_pwms_file <- "/mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-03-12/motifs.sig/motifs.adjust.diff.rna_filt.dmim/summary/ggr.pwms_patterns_summary.txt"
sig_pwms <- read.table(sig_pwms_file, sep="\t", header=TRUE, row.names=1)
sig_pwms <- rownames(sig_pwms)

# make a plot for each time point
for (task_idx in 1:nrow(activity)) {
    print(task_idx)
    task_activity <- activity[task_idx,]

    # make a new dataframe with pwm and activity
    task_data <- data.frame(t(densities))
    colnames(task_data) <- pwm_names
    task_data <- task_data[,sig_pwms]

    # attach activity
    task_data$activity <- task_activity

    differential <- FALSE
    
    if (differential) {
        # diff
        task_data$activity <- task_activity - ref_activity
        palette <- "RdBu"
    } else {
        ref_activity <- task_activity
        palette <- "Blues"
    }
    
    # melt
    task_melt <- melt(task_data, id.vars="activity")
    
    # remove data without motif density
    task_melt <- task_melt[task_melt$value != 0,]

    # remove data with wrong motif density (edges, so not 0.05 etc)
    tol <- 1e-5
    task_melt <- task_melt[(task_melt$value * pool_window)%%1 < tol,]
    
    # adjust fraction to motif counts
    task_melt$value <- task_melt$value * pool_window
    
    # aggregate
    task_agg <- aggregate(
        task_melt$activity,
        by=list(task_melt$value, task_melt$variable),
        mean)
    
    # remove data without accessibility
    task_agg <- task_agg[task_agg$x > 0,]
    all_pwms <- unique(task_agg$Group.2)
    
    # extract an ordering
    thresh <- quantile(task_agg$x, 0.60)
    task_order_data <- task_agg[task_agg$x > thresh,]
    task_order_data <- aggregate(
        task_order_data$Group.1,
        by=list(task_order_data$Group.2),
        min)
    if (task_idx == 1) {
        ordered_pwms <- as.character(task_order_data[order(task_order_data$x),]$Group.1)
        #ordered_pwms <- c(setdiff(all_pwms, ordered_pwms), ordered_pwms)
        ordered_pwms <- c(ordered_pwms, setdiff(all_pwms, ordered_pwms))
    }
    
    if (differential) {
        clip_val <- quantile(abs(task_agg$x), 0.90)
        task_agg$x[task_agg$x > clip_val] <- clip_val
        task_agg$x[task_agg$x < -clip_val] <- -clip_val
        limits <- c(-clip_val-0.01, clip_val+0.01)
    } else {
        clip_val <- quantile(task_agg$x, 0.80)
        task_agg$x[task_agg$x > clip_val] <- clip_val
        limits <- c(1.8, clip_val+0.01)
    }
    
    # and set as factor
    task_agg$motifs <- factor(task_agg$Group.2, levels=ordered_pwms)
    task_agg$motif_count <- task_agg$Group.1

    # and finally sort by activity to plot strongest last?
    task_agg <- task_agg[order(task_agg$x),]

    # and clean up names
    ordered_pwms_clean <- gsub("HCLUST-\\d+_", "", ordered_pwms)
    ordered_pwms_clean <- gsub(".UNK.0.A", "", ordered_pwms_clean)
    task_agg$motifs <- gsub("HCLUST-\\d+_", "", task_agg$motifs)
    task_agg$motifs <- gsub(".UNK.0.A", "", task_agg$motifs)
    task_agg$motifs <- factor(task_agg$motifs, levels=ordered_pwms_clean)
    
    # and plot
    plot_file <- paste("task-", task_idx, ".densities.pdf", sep="")
    ggplot(task_agg, aes(x=motif_count, y=motifs, color=x)) +
        #geom_jitter(alpha=0.7, width=0, height=0.1) +
        #geom_point(size=2) +
        geom_tile(aes(fill=x)) +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
            axis.text.y=element_text(size=5),
            axis.text.x=element_text(size=5),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.key.size=unit(0.1, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=6),
            legend.position=c(1.1, 0.9),
            strip.background=element_blank(),
            strip.text=element_blank()) +
        scale_x_continuous(limits=c(0,8), expand=c(0,0)) +
        #scale_size_continuous(limits=limits) +
        scale_colour_distiller(limits=limits, palette=palette, direction=1, na.value="white") +
        scale_fill_distiller(limits=limits, palette=palette, direction=1, na.value="white")
    ggsave(plot_file, height=5.5, width=1.5)

}



