#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read data
data <- read.table(data_file, sep="\t", header=TRUE)

if (nrow(data) > 1000) {
    step <- as.integer(nrow(data) / 1000)
    data <- data[seq(1, nrow(data), step),]
}
# make sure the final 1,0 point is there
final_point <- data.frame(precision=0, recall=1, group="PWM")
data <- rbind(data, final_point)
#final_point <- data.frame(precision=0, recall=1, group="NN")
#data <- rbind(data, final_point)

# plot
ggplot(data, aes(x=recall, y=precision, colour=group)) +
    geom_line(size=0.230) +
    labs(
        x="Recall", y="Precision", title="Precision-recall with motif scores") +
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
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0))
        
ggsave(plot_file, height=1.5, width=1.5, useDingbats=FALSE)
