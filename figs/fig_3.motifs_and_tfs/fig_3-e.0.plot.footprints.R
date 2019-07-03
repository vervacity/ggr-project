#!/usr/bin/env Rscript

# description: plot footprint files from HINT
# data: /srv/scratch/dskim89/ggr/ggr.tronn.2019-06-24.footprinting_dynamic
library(ggplot2)
library(reshape2)

# args
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
plot_file <- args[2]

# read in data and adjust as needed
data <- read.table(input_file, sep="\t", header=TRUE, row.names=1)

# normalize to flank
data_norm <- rbind(data[1:10,], data[191:200,])
flank_norm_factors <- apply(data_norm, 2, mean)

#data_norm <- rbind(data[50:75,], data[125:150,])
#center_norm_factors <- apply(data_norm, 2, mean)

#norm_factors <- center_norm_factors - flank_norm_factors
data <- data.frame(t(t(data) / flank_norm_factors))

# and then bias adjust
data_norm <- rbind(data[50:75,], data[125:150,])
norm_factors <- apply(data_norm, 2, mean)
data <- data.frame(t(t(data) - norm_factors))

data$index <- as.numeric(rownames(data)) - (nrow(data) / 2)
data <- data[50:150,]

data_melted <- melt(data, id.vars=c("index"))

# plot
ggplot(data_melted, aes(x=index, y=value, colour=variable)) +
    geom_line(size=0.115) +
    labs(title="footprint", x="Position (bp)", y="Enrichment over baseline") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=5, margin=margin(b=0)),
        plot.margin=margin(5,1,1,1),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=4),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=4),
        axis.text.x=element_text(size=4),
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
        legend.text=element_text(size=4),
        strip.background=element_blank(),
        strip.text=element_blank()) +
    scale_x_continuous(limits=c(-51,51), expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    facet_grid(. ~ variable)

ggsave(plot_file, height=2, width=4, useDingbats=FALSE)

