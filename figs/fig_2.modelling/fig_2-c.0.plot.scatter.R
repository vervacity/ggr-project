#!/usr/bin/env Rscript

# description: plot scatter plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/evals.2018-12-03

library(ggplot2)
library(ggsci)

set.seed(42)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

timepoint_colors <- get_ggr_timepoint_colors()
print(timepoint_colors)

# args
args <- commandArgs(trailingOnly=TRUE)
eval_dir <- args[1]
plot_num <- 10000 # 10k

# plot fn
plot_scatter <- function(data, plot_file, colour, day) {

    # sample
    data <- data[sample(1:nrow(data), plot_num, replace=FALSE),]
    print(dim(data))
    
    # get R vals (spearman and pearson)
    r_s <- cor(data$x, data$y, method="spearman")
    r_s <- format(round(r_s, 3), nsmall=3)
    r_s <- paste("R[S]==", r_s, sep="")
    
    r_p <- cor(data$x, data$y, method="pearson")
    r_p <- format(round(r_p, 3), nsmall=3)
    r_p <- paste("R[P]==", r_p, sep="")
    
    # plot
    ggplot(data, aes(x=x, y=y)) +
        geom_point(
            alpha=0.7, shape=16, stroke=0, size=0.3, colour=colour) + # alpha= 0.7
        annotate(
            "text", x=5, y=0.5, size=1.5, vjust="inward", hjust="inward",
            label=r_s, parse=TRUE) +
        annotate(
            "text", x=5, y=0.1, size=1.5, vjust="inward", hjust="inward",
            label=r_p, parse=TRUE) +
        geom_abline(size=0.115, intercept=0, slope=1, color="gray", linetype="dashed") +
        labs(title=day, x="Predicted", y="Actual") +
        theme_bw() +
        theme(
            aspect.ratio=1,
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=8, hjust=0.5, margin=margin(b=0)),
            plot.margin=margin(5,1,1,1),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(vjust=1, margin=margin(t=3,0,0,0)),
            axis.title.y=element_text(vjust=1, margin=margin(0,r=3,0,0)),
            axis.text.y=element_text(size=6),
            axis.text.x=element_text(size=6),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in")) +
       scale_x_continuous(limits=c(0,5), expand=c(0,0)) +
       scale_y_continuous(limits=c(0,5), expand=c(0,0))
    ggsave(plot_file, height=1.25, width=1.25, units="in", useDingbats=FALSE)

}


# extract day 0
eval_file <- paste(eval_dir, "task_0.correlation.tmp.txt", sep="/")
plot_file <- "fig_2-c.0.scatter_day-0.pdf"
data <- read.table(eval_file, header=TRUE)
plot_scatter(data, plot_file, timepoint_colors[1], "Day 0")

# extract day 3
eval_file <- paste(eval_dir, "task_6.correlation.tmp.txt", sep="/")
plot_file <- "fig_2-c.0.scatter_day-3.pdf"
data <- read.table(eval_file, header=TRUE)
plot_scatter(data, plot_file, timepoint_colors[6], "Day 3")

# extract day 6
eval_file <- paste(eval_dir, "task_12.correlation.tmp.txt", sep="/")
plot_file <- "fig_2-c.0.scatter_day-6.pdf"
data <- read.table(eval_file, header=TRUE)
plot_scatter(data, plot_file, timepoint_colors[10], "Day 6")


