#!/usr/bin/env Rscript

# description: quick metrics plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03

library(ggplot2)
library(ggsci)
library(rhdf5)

# args
args <- commandArgs(trailingOnly=TRUE)
out_prefix <- args[1]
pretrained_evals_prefixes <- args[2:11]

# just the ATAC vals
task_ids <- 0:12
subsample_n <- 50

# go through files and pull in results
rand_init <- "Fresh init."
pretrain <- "Transfer"
all_data <- data.frame()
for (i in 2:length(args)) {
    # get filename and fold
    fold_prefix <- args[i]
    name_split <- strsplit(fold_prefix, "/")[[1]][-2]
    name_split <- name_split[length(name_split) - 2]
    fold <- strsplit(name_split, "-")[[1]][-1]
    print(fold)
    
    # for each task, pull in the data
    for (task_i in task_ids) {
        filename <- paste(fold_prefix, "/task_", task_i, ".auprc.tmp.txt", sep="")
        data <- read.table(filename, header=TRUE, sep="\t")
        data <- data[seq(1, nrow(data), length.out=subsample_n),] # subsample down
        data <- rbind(data, c(1, 0))
        colnames(data) <- c("Precision", "Recall")
        data$Fold <- fold
        data$Task <- task_i
        data$Type <- pretrain
        data$fold_task <- paste(fold, "-", task_i, sep="")
        
        # save in one dataframe
        if (ncol(all_data) == 0) {
            all_data <- data
        } else {
            all_data <- rbind(all_data, data)
        }

    }
}


# plot function
plot_metric <- function(results, plot_file, title, metric_name, limits) {
    p <- ggplot(results, aes(x=Recall, y=Precision, colour=Fold, group=fold_task)) +
        geom_line(
            size=0.230) 

    p <- p +labs(title=title, x="Recall", y="Precision") +
        theme_bw() +
        theme(
            aspect.ratio=1,
            text=element_text(family="ArialMT", size=6),

            #plot.margin=margin(5,1,1,1),
            plot.margin=margin(5,10,1,5),
            plot.title=element_text(size=8, hjust=0.5, margin=margin(0,0,1,0)),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            
            axis.title=element_text(size=6),
            axis.text.y=element_text(size=6),
            axis.text.x=element_text(size=6),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),

            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.key.size=unit(0.1, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.spacing.x=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=5),
            legend.position="bottom") +
        scale_color_npg() +
        scale_fill_npg() +
        scale_x_continuous(limits=c(0.0, 1.0), expand=c(0,0)) +
        scale_y_continuous(limits=c(0.0, 1.0), expand=c(0,0))
    ggsave(plot_file, height=2, width=1.5, useDingbats=FALSE)

}

# plot
pr_file <- "fig_2-b.2.ggr_pr_curve.pdf"
plot_metric(all_data, pr_file, "Precision-Recall on \nATAC-seq peak prediction", "AUPRC", c(0.0, 1.0))
