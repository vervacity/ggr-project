#!/usr/bin/env Rscript

# description: take in deeptools matrix and plot
library(reshape2)
library(ggplot2)

# args
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
plot_file <- args[2]

# read in data
data <- read.table(data_file, sep="\t", header=TRUE)

# filler data to make it a donut
#filler_data <- data.frame(fill=0, count=1, name="center")
#data <- rbind(data, filler_data)
#print(unique(data$name))

# make NA to not show no data
data[data == 0] <- NA
data <- data[nrow(data):1,]

# levels for ordering rings
ordered_levels <- c(
    "center", # nothing
    "atac_dynamic",
    #"rna_present",
    "H3K27ac_present",
    "H3K4me1_present",
    "H3K27me3_present",
    "rna_present")
data$name <- factor(data$name, levels=ordered_levels)
#data$fill <- factor(data$fill)

# plot - CIRCOS VERSION
ggplot(data, aes(x=name, y=count, fill=fill)) +
    geom_bar(stat="identity") +
    coord_polar(theta="y", direction=-1) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(size=0.115, lineend="square"),
        panel.grid.minor=element_line(size=0.115, lineend="square"),
        #panel.grid.major=element_blank(),
        #panel.grid.minor=element_blank(),
        #axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.title=element_blank(),
        #axis.title.x=element_text(vjust=1.5),
        #axis.title.y=element_text(vjust=-1),
        axis.line=element_blank(),
        #axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        #legend.position="bottom",
        #legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0)) +            
    scale_fill_continuous(na.value=NA)
    #scale_fill_discrete(na.value=NA)
        
ggsave(plot_file, height=3, width=3)


# levels for ordering rings
data$name <- as.character(data$name)
data$name[data$name == "H3K27me3_present"] <- "H3K27me3"
data$name[data$name == "H3K27ac_present"] <- "H3K27ac"
data$name[data$name == "H3K4me1_present"] <- "H3K4me1"
data$name[data$name == "atac_dynamic"] <- "ATAC-seq"
data$name[data$name == "rna_present"] <- "TSS"

ordered_levels <- c(
    "TSS",
    "ATAC-seq",
    "H3K27ac",
    "H3K4me1",
    "H3K27me3")
ordered_levels <- c(
    "H3K27me3",
    "H3K4me1",
    "H3K27ac",
    "ATAC-seq",    
    "TSS")
data$name <- factor(data$name, levels=ordered_levels)

# plot - BAR PLOT VERSION
ggplot(data, aes(x=name, y=count, fill=fill)) +
    geom_bar(stat="identity", width=0.7, show.legend=FALSE) +
    labs(x="Assays", y="Ordered regions (n=424,700)") +
    coord_flip() +        
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,10,1,5),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        #panel.grid.major=element_line(size=0.115, lineend="square"),
        #panel.grid.minor=element_line(size=0.115, lineend="square"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        #axis.title=element_blank(),
        axis.title.x=element_text(vjust=1.5),
        axis.title.y=element_text(vjust=-1),
        #axis.line=element_blank(),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, hjust=1.0, vjust=1.0, angle=30),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        #legend.position="bottom",
        #legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.key.size=unit(0.05, "in"),
        legend.margin=margin(0,0,0,0)) +            
    scale_fill_continuous(na.value=NA) +
    scale_y_continuous(expand=c(0,0))
    #scale_fill_discrete(na.value=NA)
        
ggsave(plot_file, height=2, width=3)





