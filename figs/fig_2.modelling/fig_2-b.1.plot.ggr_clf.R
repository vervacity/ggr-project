#!/usr/bin/env Rscript

# description: quick metrics plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03

library(ggplot2)
library(ggsci)
library(rhdf5)

# args
args <- commandArgs(trailingOnly=TRUE)
out_prefix <- args[1]
random_init_evals <- args[2:11]
pretrained_evals <- args[12:21]

# go through files and pull in results
rand_init <- "Fresh init."
pretrain <- "Transfer"
for (i in 2:length(args)) {
    # get filename and fold
    filename <- args[i]
    filename_split <- strsplit(filename, "/")[[1]][-1]
    filename_split <- filename_split[length(filename_split) - 2]
    fold <- strsplit(filename_split, "-")[[1]][-1]
    
    # pull data
    data <- read.table(filename, header=TRUE, sep="\t")
    data <- data[2:14,] # change if want different (like histones)

    # prep
    auprc_results_tmp <- data.frame(value=data$AUPRC)
    auprc_results_tmp$fold <- fold
    if (i < 12) {
        auprc_results_tmp$train <- rand_init
    } else {
        auprc_results_tmp$train <- pretrain
    }
    
    auroc_results_tmp <- data.frame(value=data$AUROC)
    auroc_results_tmp$fold <- fold
    if (i < 12) {
        auroc_results_tmp$train <- rand_init
    } else {
        auroc_results_tmp$train <- pretrain
    }
    
    recall_results_tmp <- data.frame(value=data$recall_at_25_fdr)
    recall_results_tmp$fold <- fold
    if (i < 12) {
        recall_results_tmp$train <- rand_init
    } else {
        recall_results_tmp$train <- pretrain
    }

    # pull baseline AUPRC
    eval_dir <- strsplit(filename, "/by_task")[[1]][1]
    eval_file <- paste(eval_dir, "ggr.eval.h5", sep="/")
    labels <- h5read(eval_file, "ATAC_LABELS")
    baseline <- mean(labels)
    if (i < 12) {
        auprc_baseline <- data.frame(value=baseline, fold=fold, train=rand_init)
    } else {
        auprc_baseline <- data.frame(value=baseline, fold=fold, train=pretrain)
    }
        
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

# plot auprc
auprc_file <- "fig_2-b.1.ggr_auprc.pdf"
auprc_results$train <- factor(auprc_results$train, levels=c(rand_init, pretrain))
ggplot(auprc_results, aes(x=train, y=value, colour=fold)) +
    geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0) +
    geom_point(shape=16, stroke=0, size=0.3, aes(fill=fold), position=position_jitterdodge(jitter.width=0.01), show.legend=FALSE) +
    geom_point(
        data=auprc_baselines,
        aes(x=train, y=value, colour=fold),
        position=position_jitterdodge(jitter.width=0),
        shape=18,
        stroke=0, size=0.6, show.legend=FALSE) +
    labs(x="", y="AUPRC") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=5),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=5),
        legend.text=element_text(size=4)) +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(limits=c(0.1,0.8), expand=c(0,0))
ggsave(auprc_file, height=1, width=1.5, useDingbats=FALSE)


# plot auroc
auroc_file <- "fig_2-b.1.ggr_auroc.pdf"
auroc_results$train <- factor(auroc_results$train, levels=c(rand_init, pretrain))
ggplot(auroc_results, aes(x=train, y=value, colour=fold)) +
    geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0) +
    geom_point(shape=16, stroke=0, size=0.3, aes(fill=fold), position=position_jitterdodge(jitter.width=0.01), show.legend=FALSE) +
    labs(x="", y="AUROC") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=5),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=5),
        legend.text=element_text(size=4)) +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(limits=c(0.5,1.0), expand=c(0,0))
ggsave(auroc_file, height=1, width=1.5, useDingbats=FALSE)

# plot recall
recall_file <- "fig_2-b.1.ggr_recall.pdf"
recall_results$train <- factor(recall_results$train, levels=c(rand_init, pretrain))
ggplot(recall_results, aes(x=train, y=value, colour=fold)) +
    geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0) +
    geom_point(shape=16, stroke=0, size=0.3, aes(fill=fold), position=position_jitterdodge(jitter.width=0.01), show.legend=FALSE) +
    labs(x="", y="Recall at 25% FDR") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=5),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=5),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=5),
        legend.text=element_text(size=4)) +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(limits=c(0.0,0.5), expand=c(0,0))
ggsave(recall_file, height=1, width=1.5, useDingbats=FALSE)
