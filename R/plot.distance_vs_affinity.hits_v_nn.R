#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(reshape2)
library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
hits_file <- args[1]
nn_file <- args[2]
plot_file <- args[3]
title <- args[4]

# read in data
data_hits <- read.table(hits_file, header=TRUE)
data_hits$score_avg <- data_hits$start_score + data_hits$stop_score
#data_hits$score_avg <- min(data_hits$start_score, data_hits$stop_score)
data_hits$present <- 1

data_nn <- read.table(nn_file, header=TRUE)
data_nn$score_avg <- data_nn$start_score + data_nn$stop_score
#data_nn$score_avg <- min(data_nn$start_score, data_nn$stop_score)
data_nn$present <- 1

# get distances
# 1) can do mean to see if there is affinity reduction at a distance <- unlikely
# 2) can do sum to see if more sites at that distance
data_hits_summ <- aggregate(
    data_hits$score_avg, by=list(distance=data_hits$distance),
    FUN=sum)
colnames(data_hits_summ) <- c("distance", "hits")

data_nn_summ <- aggregate(
    data_nn$score_avg, by=list(distance=data_nn$distance),
    FUN=sum)
colnames(data_nn_summ) <- c("distance", "nn")

data_hits_summ <- data_hits_summ[data_hits_summ$distance > 4,]
data_nn_summ <- data_nn_summ[data_nn_summ$distance > 4,]


# normalize?
data_hits_summ$hits <- data_hits_summ$hits / sum(data_hits_summ$hits)
data_nn_summ$nn <- data_nn_summ$nn / sum(data_nn_summ$nn)

# merge
data_all <- merge(data_hits_summ, data_nn_summ, by="distance", all=TRUE)
data_all[is.na(data_all)] <- 0
data_all$diff <- data_all$nn - data_all$hits

# melt
data_melt <- melt(data_all, id.vars=c("distance"))

# plot
ggplot(data_melt, aes(x=distance, y=value)) +
    geom_hline(yintercept=0, size=0.115) +
    geom_line(size=0.230, aes(colour=variable)) +
    labs(x="Motif-motif distance",
         y="Motif-motif distribution",
         title=title) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_blank(),
        panel.grid.major.y=element_blank(),
        #panel.grid.minor=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_text(vjust=-1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        #legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0)) + 
    #scale_y_continuous(limit=c(-0.1, 0.1)) +
    scale_x_continuous(minor_breaks=seq(0,130,10), limit=c(0,140))
        
ggsave(plot_file, height=2, width=2)



