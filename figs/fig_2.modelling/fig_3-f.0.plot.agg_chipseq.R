#!/usr/bin/env Rscript

# description: plot footprint files from HINT
# data: /srv/scratch/dskim89/ggr/ggr.tronn.2019-06-24.footprinting_dynamic
library(ggplot2)
library(reshape2)

# load GGR colors
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
plot_file <- args[2]

# read in data and adjust as needed
data <- read.table(input_file, sep="\t", header=TRUE, row.names=1)
data$ymin <- data$mean - data$sem
data$ymax <- data$mean + data$sem

# adjust names
data$variable <- as.character(data$variable)
data$variable[data$variable == "pos"] <- "NN-active"
data$variable[data$variable == "neg"] <- "NN-inactive"
data$variable <- factor(data$variable, levels=c("NN-inactive", "NN-active"))

# colors
ggr_colors <- brewer.pal(6, "Paired")[3:4]

# make pretty limits
y_limit <- ceiling(max(data$mean))
x_limit <- 500

# figure out title
title_name <- strsplit(basename(plot_file), ".", fixed=TRUE)[[1]]
title_name <- title_name[8]
title_name <- paste(title_name, " ChIP-seq", sep="")

# plot
ggplot(data, aes(x=position, y=mean, colour=variable, fill=variable)) +
    geom_ribbon(alpha=0.3, size=0, colour=NA, aes(ymin=ymin, ymax=ymax), show.legend=FALSE) +
    geom_line(alpha=0.7, size=0.230) + # 0.115
    labs(title=title_name, x="Position (bp)", y="Fold enrichment") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(hjust=0.5, size=8, margin=margin(b=3)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=1, hjust=0.5, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=5, margin=margin(l=-3)),
        #legend.position=c(0.85, 0.9),
        legend.position=c(1, 0.9),
        strip.background=element_blank(),
        strip.text=element_blank()) +
    scale_color_manual(values=ggr_colors) +
    scale_fill_manual(values=ggr_colors) +
    scale_x_continuous(
        limits=c(-x_limit,x_limit),
        expand=c(0,0),
        breaks=seq(-x_limit,x_limit,by=x_limit),
        labels=c(-x_limit, 0, x_limit)) +
    scale_y_continuous(limits=c(0, y_limit), expand=c(0,0))

ggsave(plot_file, height=1, width=1, useDingbats=FALSE)
