#!/usr/bin/env Rscript

# description: quick metrics plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03

library(ggplot2)
library(ggsci)

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
    data <- data[2:14,]

    # prep
    mse_results_tmp <- data.frame(value=data$MSE)
    mse_results_tmp$fold <- fold
    if (i < 12) {
        mse_results_tmp$train <- rand_init
    } else {
        mse_results_tmp$train <- pretrain
    }
    
    spearman_results_tmp <- data.frame(value=data$SPEARMANR)
    spearman_results_tmp$fold <- fold
    if (i < 12) {
        spearman_results_tmp$train <- rand_init
    } else {
        spearman_results_tmp$train <- pretrain
    }
    
    pearson_results_tmp <- data.frame(value=data$PEARSONR)
    pearson_results_tmp$fold <- fold
    if (i < 12) {
        pearson_results_tmp$train <- rand_init
    } else {
        pearson_results_tmp$train <- pretrain
    }
    
    # append
    if (i == 2) {
        mse_results <- mse_results_tmp
        spearman_results <- spearman_results_tmp
        pearson_results <- pearson_results_tmp
    } else {
        mse_results <- rbind(mse_results, mse_results_tmp)
        spearman_results <- rbind(spearman_results, spearman_results_tmp)
        pearson_results <- rbind(pearson_results, pearson_results_tmp)
    }

}

# plot mse
mse_file <- "fig_2-b.2.ggr_mse.pdf"
mse_results$train <- factor(mse_results$train, levels=c(rand_init, pretrain))
ggplot(mse_results, aes(x=train, y=value, colour=fold)) +
    geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0) +
    geom_point(shape=16, stroke=0, size=0.3, aes(fill=fold), position=position_jitterdodge(jitter.width=0.01), show.legend=FALSE) +
    labs(x="", y="MSE") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=5),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=4),
        legend.text=element_text(size=4)) +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(limits=c(0.4,1.2), expand=c(0,0))
ggsave(mse_file, height=1, width=1.5, useDingbats=FALSE)


# plot spearman
spearman_file <- "fig_2-b.2.ggr_spearman.pdf"
spearman_results$train <- factor(spearman_results$train, levels=c(rand_init, pretrain))
ggplot(spearman_results, aes(x=train, y=value, colour=fold)) +
    geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0) +
    geom_point(shape=16, stroke=0, size=0.3, aes(fill=fold), position=position_jitterdodge(jitter.width=0.01), show.legend=FALSE) +
    labs(title="Keratinocyte ATAC signals", x="", y="Spearman R") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT", size=6),
        plot.margin=margin(5,1,1,1),
        plot.title=element_text(size=6, margin=margin(0,0,0,0)),
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
        #legend.title=element_text(size=4),
        #legend.text=element_text(size=4)) +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(limits=c(0.50, 0.80), expand=c(0,0))
ggsave(spearman_file, height=1, width=1.5, useDingbats=FALSE)

# plot pearson
pearson_file <- "fig_2-b.2.ggr_pearson.pdf"
pearson_results$train <- factor(pearson_results$train, levels=c(rand_init, pretrain))
ggplot(pearson_results, aes(x=train, y=value, colour=fold)) +
    geom_boxplot(size=0.1, outlier.size=0, outlier.stroke=0) +
    geom_point(shape=16, stroke=0, size=0.3, aes(fill=fold), position=position_jitterdodge(jitter.width=0.01), show.legend=FALSE) +
    labs(x="", y="Pearson R") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=5),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.key.size=unit(0.01, "in"),
        legend.margin=margin(5,0,0,0),
        legend.title=element_text(size=4),
        legend.text=element_text(size=4)) +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(limits=c(0.45, 0.80), expand=c(0,0))
ggsave(pearson_file, height=1.5, width=1.5, useDingbats=FALSE)
