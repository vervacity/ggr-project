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
plot_width <- 1.25
plot_height <- 1.5
divide_count <- 100
num_type <- "regions"

if (grepl("atac", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Blues", 50)
    title <- "ATAC-seq"
    plot_width <- 1.75
    y_static_lim <- 160
    y_dynamic_lim <- 15
} else if (grepl("H3K27ac", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Reds", 50)
    title <- "H3K27ac ChIP-seq"
} else if (grepl("H3K4me1", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Oranges", 50)
    title <- "H3K4me1 ChIP-seq"
} else if (grepl("H3K27me3", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Greens", 50)
    title <- "H3K27me3 ChIP-seq"
    y_dynamic_lim_neg <- -12
    y_dynamic_lim_pos <- 2
    plot_height <- 1.25
} else if (grepl("hichip", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Purples", 50)
    title <- "HiChIP"
} else if (grepl("rna", plot_file, fixed=TRUE)) {
    color <- get_ggr_assay_palette("Purples", 50)
    title <- "PAS-seq"
    plot_width <- 1.75
    num_type <- "genes"
    y_static_lim <- 12
    y_dynamic_lim <- 2
}


# keep last color
color <- color[length(color)]
data$color <- color

# now set up differences for static vs dynamic counts
if (count_type == "static") {
    y <- paste("Number of ", num_type, "\n(1000s)", sep="")
    y_lims <- c(0, y_static_lim)
    data$count <- data$count / 1000.
} else if (count_type == "dynamic") {
    y <- paste("Number of changing\n", num_type, " (100s)", sep="")
    y_lims <- c(y_dynamic_lim_neg, y_dynamic_lim_pos)
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
    data$count <- data$value / divide_count
    
}

# factorize timepoints
data$timepoint <- as.factor(data$timepoint)
data$color <- as.factor(data$color)
print(data)

# plot
ggplot(data, aes(x=timepoint, y=count, fill=color, colour=color)) +
    geom_bar(
        alpha=0.7,
        stat="identity",
        position=position_stack(),
        width=0.4,
        size=0.115,
        #aes(fill=color, colour=color),
        show.legend=FALSE) +
    geom_hline(yintercept=0, size=0.115) +
    labs(x="Timepoint (day)", y=y, title=title) +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(b=3)),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title=element_text(size=6, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.text=element_text(size=6),
        axis.text.x=element_text(angle=60, hjust=1)) +
    scale_y_continuous(limits=y_lims, expand=c(0,0)) +
    scale_fill_manual(values=levels(data$color)) +
    scale_colour_manual(values=levels(data$color))

ggsave(plot_file, height=plot_height, width=plot_width, useDingbats=FALSE) # height 1.25
