#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
enrichments_file <- args[1]
motif_name <- args[2]
plot_file <- args[3]

# read in data
data <- read.table(enrichments_file, sep="\t", header=TRUE)
data$term.name <- factor(data$term.name, levels=data$term.name)

# plot
ggplot(data, aes(x=term.name, y=X.log10pval)) +
    geom_bar(stat="identity") +
    labs(x="GO term", y="-log10pval", title=motif_name) +
    coord_flip() +
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
        legend.position="bottom")
    
ggsave(plot_file, height=1.5, width=3, useDingbats=FALSE)
