#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]
title <- args[3]

# read data
data <- read.table(data_file, sep="\t", header=TRUE)

denominator <- sum(data$labels)
print(denominator)

# rank data
rank_thresh <- 20000
data_nn <- data[order(data$NN, decreasing=TRUE),]
data_nn <- data_nn[1:rank_thresh,]
data_pwm <- data[order(data$PWM, decreasing=TRUE),]
data_pwm <- data_pwm[1:rank_thresh,]

# make a CDF
i_vals <- c()
cum_sums_nn <- c()
cum_sums_pwm <- c()
for (i in seq(1,rank_thresh, 100)) {
    i_vals <- c(i_vals, i)
    cum_sums_nn <- c(cum_sums_nn, sum(data_nn$labels[1:i]))
    cum_sums_pwm <- c(cum_sums_pwm, sum(data_pwm$labels[1:i]))
    
}

data <- data.frame(rank=i_vals, nn=cum_sums_nn, pwm=cum_sums_pwm)
data_melt <- melt(data, id.vars="rank")
data_melt$variable <- as.character(data_melt$variable)
data_melt$variable[data_melt$variable == "nn"] <- "Neural motif score"
data_melt$variable[data_melt$variable == "pwm"] <- "PWM score"
data_melt$variable <- as.factor(data_melt$variable)

# plot
ggplot(data_melt, aes(x=rank, y=value, colour=variable)) +
    geom_line(size=0.230) +
    labs(
        x="Ranked motif instances", y="ChIP-seq region overlap", title=title) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_text(vjust=-1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0)) +
    scale_x_continuous(breaks=seq(0,20000,by=20000))
        
ggsave(plot_file, height=1.5, width=1.25, useDingbats=FALSE)
