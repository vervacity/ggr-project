#!/usr/bin/env Rscript

# description: plot scatter plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03

library(ggplot2)
library(ggsci)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

timepoint_colors <- get_ggr_timepoint_colors()
print(timepoint_colors)

# args
args <- commandArgs(trailingOnly=TRUE)
eval_dir <- args[1]

# extract day 0
eval_file <- paste(eval_dir, "task_0.correlation.tmp.txt", sep="/")
data <- read.table(eval_file, header=TRUE)
r_val <- cor(data$x, data$y, method="spearman")
r_val <- format(round(r_val, 3), nsmall=3)
r_val <- paste("R=", r_val, sep="")
r_val <- data.frame(x=5, y=0.1,text=r_val)
plot_file <- "fig_2-c.0.scatter_day-0.pdf"
ggplot(data, aes(x=x, y=y)) +
    geom_point(
        alpha=0.1, shape=16, stroke=0, size=0.3, colour=timepoint_colors[1]) + # alpha= 0.05
    geom_text(data=r_val, aes(label=text), size=1.5, vjust="inward", hjust="inward") +
    labs(x="Predicted", y="Actual") +
    theme_bw() +
    theme(
        aspect.ratio=1,
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
    scale_x_continuous(limits=c(0,5), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,5), expand=c(0,0))
ggsave(plot_file, height=1, width=1, useDingbats=FALSE)

# extract day 3
eval_file <- paste(eval_dir, "task_6.correlation.tmp.txt", sep="/")
data <- read.table(eval_file, header=TRUE)
r_val <- cor(data$x, data$y, method="spearman")
r_val <- format(round(r_val, 3), nsmall=3)
r_val <- paste("R=", r_val, sep="")
r_val <- data.frame(x=5, y=0.1,text=r_val)
plot_file <- "fig_2-c.0.scatter_day-3.pdf"
ggplot(data, aes(x=x, y=y)) +
    geom_point(
        alpha=0.1, shape=16, stroke=0, size=0.3, colour=timepoint_colors[7]) +
    geom_text(data=r_val, aes(label=text), size=1.5, vjust="inward", hjust="inward") +
    labs(x="Predicted", y="Actual") +
    theme_bw() +
    theme(
        aspect.ratio=1,
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
    scale_x_continuous(limits=c(0,5), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,5), expand=c(0,0))
ggsave(plot_file, height=1, width=1, useDingbats=FALSE)

# extract day 6
eval_file <- paste(eval_dir, "task_12.correlation.tmp.txt", sep="/")
data <- read.table(eval_file, header=TRUE)
r_val <- cor(data$x, data$y, method="spearman")
r_val <- format(round(r_val, 3), nsmall=3)
r_val <- paste("R=", r_val, sep="")
r_val <- data.frame(x=5, y=0.1,text=r_val)
plot_file <- "fig_2-c.0.scatter_day-6.pdf"
ggplot(data, aes(x=x, y=y)) +
    geom_point(
        alpha=0.1, shape=16, stroke=0, size=0.3, colour=timepoint_colors[10]) +
    geom_text(data=r_val, aes(label=text), size=1.5, vjust="inward", hjust="inward") +
    labs(x="Predicted", y="Actual") +
    theme_bw() +
    theme(
        aspect.ratio=1,
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
    scale_x_continuous(limits=c(0,5), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,5), expand=c(0,0))
ggsave(plot_file, height=1, width=1, useDingbats=FALSE)

