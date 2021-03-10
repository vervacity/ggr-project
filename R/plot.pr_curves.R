#!/usr/bin/env Rscript

# Description:
# plot PR curves

library(ggplot2)
library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_prefix <- args[2]


data <- read.table(data_file, sep="\t", header=TRUE)
data$fold <- as.factor(data$fold)

# plot
my_palette <- rev(colorRampPalette(brewer.pal(11, 'Spectral'))(length(unique(data$fold))))
ggplot(data, aes(x=recall, y=precision, group=fold)) +
    geom_line(aes(color=fold)) +
    labs(x="Recall", y="Precision", title="Genome-wide enhancer prediction") +
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
    theme_bw() +
    coord_fixed() +
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
        axis.title.y=element_text(vjust=1),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
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
    scale_colour_manual(values=my_palette)

out_file <- paste(plot_prefix, ".pdf", sep="")
ggsave(out_file, height=2, width=2, useDingbats=FALSE)


