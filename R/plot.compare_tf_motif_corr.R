#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(data_file)

# jitter the key TFs so that they will be present
data$Homer.ATAC[rownames(data)=="TP63"]  <- data$Homer.ATAC[rownames(data)=="TP63"] + 0.01
data$Homer.ATAC[rownames(data)=="ZNF750"]  <- data$Homer.ATAC[rownames(data)=="ZNF750"] + 0.01
data$Homer.ATAC[rownames(data)=="KLF4"]  <- data$Homer.ATAC[rownames(data)=="KLF4"] + 0.01
data$Homer.ATAC[rownames(data)=="MAF"]  <- data$Homer.ATAC[rownames(data)=="MAF"] + 0.01
data$Homer.ATAC[rownames(data)=="MAFB"]  <- data$Homer.ATAC[rownames(data)=="MAFB"] + 0.01

# ignore those for which homer+atac is zero
data <- data[data$Homer.ATAC != 0,]

# label those that HOMER did better on
data$labels <- rownames(data)
data$labels[data$Homer.ATAC <= data$NN] <- ""

# also label key TFs
data$labels[rownames(data) == "TP63"] <- "TP63"
data$labels[rownames(data) == "ZNF750"] <- "ZNF750"
data$labels[rownames(data) == "KLF4"] <- "KLF4"
data$labels[rownames(data) == "MAFB"] <- "MAFB"
data$labels[rownames(data) == "MAF"] <- "MAF"

data$labels[rownames(data) == "CREB1"] <- "CREB1"
data$labels[rownames(data) == "E2F7"] <- "E2F7"
data$labels[rownames(data) == "NFATC1"] <- "NFATC1"
data$labels[rownames(data) == "TEAD4"] <- "TEAD4"

# plot
ggplot(data, aes(x=Homer.ATAC, y=NN)) +
    geom_point(shape=21, fill=NA) +
    geom_text(aes(x=Homer.ATAC+0.04, label=labels), size=1.5, hjust=0) +
    geom_abline(intercept=0, slope=1, size=0.115) +
    labs(
        title="TF-Motif correlations by Neural Net\nMotifs vs Overlying Accessibility",
        x="Accessibility to TF expression\ncorrelation",
        y="NN motifs to TF expression\ncorrelation") +
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
        axis.ticks.length=unit(0.01, "in")) +
    scale_x_continuous(expand=c(0,0), limits=c(0, 1)) +
    scale_y_continuous(expand=c(0,0), limits=c(0, 1))

 
ggsave(plot_file, height=2, width=2, useDingbats=FALSE)
