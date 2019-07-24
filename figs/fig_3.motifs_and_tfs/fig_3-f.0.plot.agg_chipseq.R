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

# normalize to flank
data_norm <- rbind(data[1:10,], data[(nrow(data)-9):nrow(data),])
flank_norm_factors <- apply(data_norm, 2, mean)
data <- data.frame(t(t(data) / flank_norm_factors))

# index
data$index <- as.numeric(rownames(data)) - (nrow(data) / 2)

# melt
data_melted <- melt(data, id.vars=c("index"))

# adjust colors if actually just pos/neg
if (dim(data)[2] == 3) {
    ggr_colors <- rev(brewer.pal(6, "Paired")[1:2])
    data_melted$variable <- factor(data_melted$variable, levels=c("pos", "neg"))
} else {
    ggr_colors <- get_ggr_timepoint_colors()
    ggr_colors <- ggr_colors[c(2, 7, 10)]
    ggr_colors <- rev(ggr_colors)
    data_melted$variable <- factor(data_melted$variable, levels=c("d6.0", "d3.0", "d0.0"))
}

# make pretty limits
y_limit <- ceiling(2*max(abs(data_melted$value))) / 2
x_limit <- abs(min(data$index))

# figure out title
title_name <- strsplit(basename(plot_file), ".", fixed=TRUE)[[1]]
title_name <- title_name[2]

# plot
ggplot(data_melted, aes(x=index, y=value, colour=variable)) +
    geom_line(alpha=0.7, size=0.230) + # 0.115
    labs(title=title_name, x="Position (bp)", y="Read coverage") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
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
        legend.position=c(0.9, 0.9),
        strip.background=element_blank(),
        strip.text=element_blank()) +
    scale_color_manual(values=ggr_colors) +
    scale_x_continuous(limits=c(-x_limit,x_limit), expand=c(0,0)) +
    scale_y_continuous(limits=c(0, y_limit), expand=c(0,0)) #+
    #facet_grid(. ~ variable)

ggsave(plot_file, height=1.2, width=1, useDingbats=FALSE)

