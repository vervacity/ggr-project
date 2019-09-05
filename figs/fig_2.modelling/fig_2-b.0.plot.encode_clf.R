#!/usr/bin/env Rscript

# description: quick metrics plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03

library(ggplot2)
library(ggsci)
library(rhdf5)

# args
args <- commandArgs(trailingOnly=TRUE)
out_prefix <- args[1]

# go through files and pull in results
for (i in 2:length(args)) {
    # get filename and fold
    filename <- args[i]
    filename_split <- strsplit(filename, "/")[[1]][-1]
    filename_split <- filename_split[length(filename_split) - 2]
    fold <- strsplit(filename_split, "-")[[1]][-1][-1]
    
    # pull data
    data <- read.table(filename, header=TRUE, sep="\t")
    data <- data[2:431,] # just DNase
    
    # prep
    auprc_results_tmp <- data.frame(value=data$AUPRC)
    auprc_results_tmp$fold <- fold
    
    auroc_results_tmp <- data.frame(value=data$AUROC)
    auroc_results_tmp$fold <- fold
    
    recall_results_tmp <- data.frame(value=data$recall_at_25_fdr)
    recall_results_tmp$fold <- fold

    # now grab the eval data to get baseline AUPRC
    eval_dir <- strsplit(filename, "/by_task")[[1]][1]
    eval_file <- paste(eval_dir, "encode-roadmap.eval.h5", sep="/")
    labels <- h5read(eval_file, "DNASE_LABELS")
    if (!is.null(labels)) {
        baseline <- mean(labels)
    }
    auprc_baseline <- data.frame(value=baseline, fold=fold)
    
    # append
    if (i == 2) {
        auprc_results <- auprc_results_tmp
        auprc_baselines <- auprc_baseline
        auroc_results <- auroc_results_tmp
        recall_results <- recall_results_tmp
    } else {
        auprc_results <- rbind(auprc_results, auprc_results_tmp)
        auprc_baselines <- rbind(auprc_baselines, auprc_baseline)
        auroc_results <- rbind(auroc_results, auroc_results_tmp)
        recall_results <- rbind(recall_results, recall_results_tmp)
    }

}

# plot function
plot_metric <- function(results, plot_file, title, metric_name, limits) {
    p <- ggplot(results, aes(x=fold, y=value, colour=fold)) +
        geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0, show.legend=FALSE) +
        geom_point(
            shape=16,
            stroke=0,
            size=0.3,
            aes(fill=fold),
            position=position_jitterdodge())
    if (metric_name == "AUPRC") {
        p <- p + geom_point(
            data=auprc_baselines,
            aes(x=fold, y=value, colour=fold),
            position=position_jitterdodge(jitter.width=0),
            shape=18,
            stroke=0, size=0.6, show.legend=FALSE)
    }
    p <- p +labs(title=title, x="Fold", y=metric_name) +
        theme_bw() +
        theme(
            text=element_text(family="ArialMT", size=6),
            plot.margin=margin(5,1,1,1),
            plot.title=element_text(size=8, margin=margin(0,0,0,0)),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_text(size=6),
            axis.text.y=element_text(size=6),
            axis.text.x=element_text(size=6),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            legend.key.size=unit(0.01, "in"),
            legend.margin=margin(5,0,0,0)) +
        scale_color_npg() +
        scale_fill_npg() +
        scale_y_continuous(limits=limits, expand=c(0,0))
    ggsave(plot_file, height=1.75, width=1.5, useDingbats=FALSE)

}


# plot auprc
auprc_file <- "fig_2-b.0.encode_auprc.pdf"
plot_metric(auprc_results, auprc_file, "ENCODE/Roadmap\nDNase peaks", "AUPRC", c(0.0,0.8))

# auroc
auroc_file <- "fig_2-b.0.encode_auroc.pdf"
plot_metric(auroc_results, auroc_file, "ENCODE/Roadmap\nDNase peaks", "AUROC", c(0.5,1.0))

# recall
recall_file <- "fig_2-b.0.encode_recall.pdf"
plot_metric(recall_results, recall_file, "ENCODE/Roadmap\nDNase peaks", "Recall at 25% FDR", c(0,0.8))
