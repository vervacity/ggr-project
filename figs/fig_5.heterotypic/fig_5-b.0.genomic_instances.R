#!/usr/bin/env Rscript

library(rhdf5)
library(reshape2)
library(ggplot2)
library(viridis)

# plot functions

plot_motif_counts <- function(data, plot_file, title) {

    ggplot(data, aes(x=dist, stat(count))) + 
        geom_density(adjust=1/10) +
        labs(title=title, x="Position (PWM centered)", y="Count") +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major.y=element_blank(),
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
        scale_x_continuous(limits=c(-160, 160), expand=c(0,0), breaks=seq(-160,160, 10))
    
    ggsave(plot_file, height=1, width=10)

}

plot_activity_density <- function(data, plot_file, title) {

    ggplot(data, aes(x=dist, y=activity)) +
        geom_point(stroke=0, shape=20, alpha=0.05, size=1) +
        labs(title=title, x="Position (PWM centered)", y="Activity") +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.minor=element_blank(),
            panel.grid.major.y=element_blank(),
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
        scale_x_continuous(limits=c(-160, 160), expand=c(0,0), breaks=seq(-160,160, 10))
    
    ggsave(plot_file, height=1, width=10)

}




# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
grammar_file <- args[2]
min_spacing <- 12
min_activity <- 1
activity_keys <- c("ATAC_SIGNALS.NORM", "H3K27ac_SIGNALS.NORM")

# read in data file and get pwm names
weighted_pwm_positions <- h5read(
    data_file, "sequence-weighted.active.pwm-scores.thresh.max.idx", read.attributes=TRUE) # {1, M, N}
raw_pwm_positions <- h5read(
    data_file, "sequence.active.pwm-scores.thresh.max.idx", read.attributes=TRUE) # {1, M, N}
pwm_names <- attr(weighted_pwm_positions, "pwm_names")

# read in grammar file
grammars <- read.table(grammar_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)

# go through each grammar
for (i in 1:nrow(grammars)) {

    # prefix
    grammar_prefix <- strsplit(basename(grammars$filename), ".gml")[[1]]
    print(grammar_prefix)
    
    # extract nodes
    nodes <- strsplit(grammars$nodes[i], ",")[[1]]
    title <- grammars$nodes[i]
    print(nodes)

    # get pwm indices
    pwm_indices <- c()
    for (i in 1:length(nodes)) {
        pwm_str <- nodes[i]
        pwm_str <- paste(pwm_str, ".", sep="")
        for (j in 1:length(pwm_names)) {
            if (grepl(pwm_str, pwm_names[j], fixed=TRUE)) {
                pwm_indices <- c(pwm_indices, j)
            }
        }
    }
 
    # run analysis with pwms
    pwm_features <- data.frame(t(weighted_pwm_positions[1,pwm_indices,]))

    # filter
    pwm_features <- pwm_features[pwm_features$X1 != 432,]
    pwm_features <- pwm_features[pwm_features$X2 != 432,]
    pwm_features <- pwm_features[abs(pwm_features$X1 - pwm_features$X2) > min_spacing,]
    grammar_indices <- as.numeric(rownames(pwm_features))

    # calculate dist (centered on 1st pwm)
    pwm_features$dist <- pwm_features$X2 - pwm_features$X1

    # plot counts
    plot_file <- paste("fig_5-b.", grammar_prefix, ".spacing.pdf", sep="")
    print(plot_file)
    plot_motif_counts(pwm_features, plot_file, title)
    
    # plot for each activity phenotype (ATAC, H3K27ac)
    for (activity_idx in 1:length(activity_keys)) {
        activity_key <- activity_keys[activity_idx]
        print(activity_key)

        activity <- t(h5read(data_file, activity_key))
        activity <- activity[grammar_indices,]
        activity_prefix <- paste(grammar_prefix, activity_key, sep=".")
        
        # for now, plot each task
        for (task_idx in 1:ncol(activity)) {
            task_prefix <- paste(activity_prefix, task_idx, sep=".")

            # filter
            task_features <- pwm_features
            print(dim(task_features))
            task_features$activity <- activity[,task_idx]
            task_features <- task_features[task_features$activity > min_activity,]
            print(dim(task_features))
            
            plot_file <- paste("fig_5-b.", task_prefix, ".spacing.pdf", sep="")
            print(plot_file)
            plot_activity_density(task_features, plot_file, title)
        }
    }
}
