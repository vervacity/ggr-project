#!/usr/bin/env Rscript

# description: plot GSEA summary
library(ggplot2)
library(viridis)
library(ggsci)
library(colorspace)

library(reshape2)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
plot_file <- args[2]
count_type <- args[3]

# read in data
data <- read.table(gzfile(counts_file), header=TRUE, sep="\t")
#print(data)

# assay adjustments
plot_width <- 1.5
if (grepl("atac", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Blues", 50)
    title <- "ATAC-seq"
    plot_width <- 2
} else if (grepl("H3K27ac", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Reds", 50)
    title <- "H3K27ac ChIP-seq"
} else if (grepl("H3K4me1", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Oranges", 50)
    title <- "H3K4me1 ChIP-seq"
} else if (grepl("H3K27me3", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Greens", 50)
    title <- "H3K27me3 ChIP-seq"
} else if (grepl("rna", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Purples", 50)
    title <- "PAS-seq"
    plot_width <- 2
}


# keep last color
color <- color[length(color)]
data$color <- color

# now set up differences for static vs dynamic counts
if (count_type == "static") {
    y <- "Number of regions\n(1000s)"
    data$count <- data$count / 1000.
} else if (count_type == "dynamic") {
    y <- "Number of changing regions\n(1000s)"
    colnames(data) <- c("timepoint", "up", "down", "color")
    # adjust data here
    data$down <- -data$down
    data <- melt(data, id.vars=c("timepoint", "color"))
    data$timepoint <- gsub("^.+_to_", "", data$timepoint)
    #data$timepoint <- gsub("_", " ", data$timepoint)
    data$timepoint <- gsub("0$", ".0", data$timepoint)
    data$timepoint <- gsub("5$", ".5", data$timepoint)
    #data$timepoint <- gsub("0 ", ".0 ", data$timepoint)
    #data$timepoint <- gsub("5 ", ".5 ", data$timepoint)
    data$timepoint <- gsub("d", "", data$timepoint)

    # adjust darker colors
    data$color[data$variable == "down"] <- darken(color, amount=0.2)
    data <- rbind(data, c("0", color, "up", integer(0)))
    data$value <- as.numeric(data$value)

    # get count fract
    data$count <- data$value / 1000.
    
}

# factorize timepoints
data$timepoint <- as.factor(data$timepoint)
data$color <- as.factor(data$color)
print(data)

# plot
ggplot(data, aes(x=timepoint, y=count)) +
    geom_bar(
        stat="identity",
        position=position_stack(),
        width=0.25,
        aes(fill=color),
        show.legend=FALSE) +
    labs(x="Timepoint (day)", y=y, title=title) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=1)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.text=element_text(size=6),
        axis.text.x=element_text(angle=60, hjust=1)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_fill_manual(values=levels(data$color))

ggsave(plot_file, height=1.25, width=plot_width, useDingbats=FALSE) # width 14
