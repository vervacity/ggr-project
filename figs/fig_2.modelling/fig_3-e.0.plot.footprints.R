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
out_prefix <- args[2]

# params: null flanks, for flank normalization
null_flank_start <- 1 # how far from edge, start
null_flank_end <- 50 # how far from edge, end

# params: accessibility flanks, for calculating accessibility diff
acc_flank_start <- 100 # how far from MIDDLE, start
acc_flank_end <- 15 # how from MIDDLE, end. MAKE SURE DOES NOT OVERLAP FOOTPRINT

# params: accessibility flank norm widths, for normalizing for footpring depth
acc_flank_norm_start <- 40 # how far from MIDDLE, start (prev 75)
acc_flank_norm_end <- 15 # how far from MIDDLE, end (prev 25)

# params: location of footprint
footprint_start <- -8
footprint_end <- -1

# read in data, create index
data <- read.table(input_file, sep="\t", header=TRUE, row.names=1)
if (dim(data)[2] == 2) {
    colnames(data) <- c("NN.inactive", "NN.active")
}
data$position <- as.numeric(rownames(data)) - (nrow(data) / 2)

# normalizations: flank norm. this gives a fold change enrichment over background
l_null_flank_start <- data$position[null_flank_start]
l_null_flank_end <- data$position[null_flank_start] + null_flank_end
r_null_flank_start <- data$position[nrow(data)] - null_flank_start + 1
r_null_flank_end <- data$position[nrow(data)] - null_flank_end + 1

background_vals <- rbind(
    data[(data$position >= l_null_flank_start) & (data$position < l_null_flank_end),],
    data[(data$position <= r_null_flank_start) & (data$position > r_null_flank_end),])
background_avg_vals <- apply(background_vals, 2, mean)

data_flank_adj <- data.frame(t(t(data) / background_avg_vals))
data_flank_adj$position <- as.numeric(rownames(data_flank_adj)) - (nrow(data_flank_adj) / 2)

# melt
data_melted <- melt(data_flank_adj, id.vars=c("position"))

# adjust colors if actually just pos/neg
if (dim(data)[2] == 3) {
    data_melted$variable <- gsub(".", "-", data_melted$variable, fixed=TRUE)
    ggr_colors <- brewer.pal(6, "Paired")[1:2]
    data_melted$variable <- factor(data_melted$variable, levels=c("NN-inactive", "NN-active"))
} else {
    ggr_colors <- get_ggr_timepoint_colors()
    ggr_colors <- ggr_colors[c(2, 7, 10)]
    #ggr_colors <- rev(ggr_colors)
    #data_melted$variable <- factor(data_melted$variable, levels=c("d6.0", "d3.0", "d0.0"))
    data_melted$variable <- factor(data_melted$variable, levels=c("d0.0", "d3.0", "d6.0"))
}

# figure out title
title_name <- strsplit(basename(out_prefix), ".", fixed=TRUE)[[1]]
title_name <- gsub("HCLUST-\\d+_", "", title_name[2])

# ==================================================
# PLOT 1 - plot full pattern across all positions
# ==================================================
plot_file <- paste(out_prefix, ".full_length.pdf", sep="")

# make pretty y limit
#y_max <- ceiling(2*max(data_melted$value)) / 2
#y_min <- floor(2*min(data_melted$value)) / 2
y_max <- ceiling(max(data_melted$value))
y_min <- min(floor(min(data_melted$value)), 0)

# plot
ggplot(data_melted, aes(x=position, y=value, colour=variable)) +
    #geom_rect(alpha=0.2, size=0, color=NA, fill="grey95", aes(xmin=footprint_start, xmax=footprint_end, ymin=y_min, ymax=y_max)) +
    geom_rect(alpha=0.2, size=0, color=NA, fill="grey95", aes(xmin=-acc_flank_start, xmax=-acc_flank_end, ymin=y_min, ymax=y_max)) +
    geom_rect(alpha=0.2, size=0, color=NA, fill="grey95", aes(xmin=acc_flank_end, xmax=acc_flank_start, ymin=y_min, ymax=y_max)) +
    geom_line(alpha=0.7, size=0.230) + # 0.115
    labs(title=paste(title_name, " motif footprinting", sep=""), x="\nPosition (bp)", y="Fold enrichment\n") +
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
        axis.title.y=element_text(vjust=0, hjust=0.5, margin=margin(0,0,0,0)),
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
        legend.position=c(0.9, 0.9),
        strip.background=element_blank(),
        strip.text=element_blank()) +
    scale_color_manual(values=ggr_colors) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(limits=c(y_min, y_max), expand=c(0,0))

ggsave(plot_file, height=1.2, width=2, useDingbats=FALSE)


# ==================================================
# PLOT 2 - plot accessibility diff
# ==================================================
plot_file <- paste(out_prefix, ".accessibility_diff.pdf", sep="")

l_acc_flank_start <- -acc_flank_start #data$position[acc_flank_start]
l_acc_flank_end <- -acc_flank_end #footprint_start #data$position[acc_flank_start] + acc_flank_end
r_acc_flank_start <- acc_flank_end #data$position[nrow(data)] - acc_flank_start + 1
r_acc_flank_end <- acc_flank_start #footprint_end #data$position[nrow(data)] - acc_flank_end + 1

acc_flank_vals <- rbind(
    data_flank_adj[(data_flank_adj$position >= l_acc_flank_start) & (data_flank_adj$position < l_acc_flank_end),],
    data_flank_adj[(data_flank_adj$position <= r_acc_flank_start) & (data_flank_adj$position > r_acc_flank_end),])
acc_flank_avg_vals <- as.data.frame(t(apply(acc_flank_vals, 2, mean)))
acc_flank_avg_vals$position <- NULL
acc_flank_avg_vals <- melt(acc_flank_avg_vals)


if (dim(data)[2] == 3) {
    acc_flank_avg_vals$variable <- gsub(".", "-", acc_flank_avg_vals$variable, fixed=TRUE)
    acc_flank_avg_vals$variable <- factor(acc_flank_avg_vals$variable, levels=c("NN-inactive", "NN-active"))
} else {
    #acc_flank_avg_vals$variable <- factor(acc_flank_avg_vals$variable, levels=c("d6.0", "d3.0", "d0.0"))
    acc_flank_avg_vals$variable <- factor(acc_flank_avg_vals$variable, levels=c("d0.0", "d3.0", "d6.0"))
}
print(acc_flank_avg_vals)

# pretty y scale
y_max <- ceiling(max(acc_flank_avg_vals$value))
y_min <- min(floor(min(acc_flank_avg_vals$value)), 0)

# plot
ggplot(acc_flank_avg_vals, aes(x=variable, y=value, colour=variable, fill=variable)) +
    geom_col(alpha=0.5, size=0.115, width=0.6, show.legend=FALSE) +
    labs(title=title_name, x="", y="Flank ATAC enrichment\n") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(hjust=0.5, size=8, margin=margin(b=1)),
        plot.margin=margin(5,0,-8,0),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=0, hjust=1, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, vjust=1, hjust=1, angle=60),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in")) +
    scale_color_manual(values=ggr_colors) +
    scale_fill_manual(values=ggr_colors) +
    scale_y_continuous(limits=c(y_min, y_max), expand=c(0,0))

#ggsave(plot_file, height=1.2, width=0.75, useDingbats=FALSE)
ggsave(plot_file, height=1, width=0.75, useDingbats=FALSE)



# ==================================================
# PLOT 3 - plot full pattern across all positions
# ==================================================
plot_file <- paste(out_prefix, ".footprint_w_flanks.pdf", sep="")

# TODO use the flank adjusted?
data <- data_flank_adj

# bias adjust
l_bias_flank_start <- -acc_flank_norm_start
l_bias_flank_end <- -acc_flank_norm_end
r_bias_flank_start <- acc_flank_norm_end
r_bias_flank_end <- acc_flank_norm_start

bias_vals <- rbind(
    data[(data$position >= l_bias_flank_start) & (data$position < l_bias_flank_end),],
    data[(data$position <= r_bias_flank_start) & (data$position > r_bias_flank_end),])
bias_avg_vals <- apply(bias_vals, 2, mean)

data_bias_adj <- data.frame(t(t(data) - bias_avg_vals))
data_bias_adj$position <- as.numeric(rownames(data_bias_adj)) - (nrow(data_bias_adj) / 2)

x_lim <- acc_flank_norm_start

# melt
data_bias_adj_melted <- melt(data_bias_adj, id.vars=c("position"))
if (dim(data)[2] == 3) {
    data_bias_adj_melted$variable <- gsub(".", "-", data_bias_adj_melted$variable, fixed=TRUE)
    data_bias_adj_melted$variable <- factor(data_bias_adj_melted$variable, levels=c("NN-inactive", "NN-active"))
} else {
    #data_bias_adj_melted$variable <- factor(data_bias_adj_melted$variable, levels=c("d6.0", "d3.0", "d0.0"))
    data_bias_adj_melted$variable <- factor(data_bias_adj_melted$variable, levels=c("d0.0", "d3.0", "d6.0"))
}

# make pretty y limit
#y_max <- ceiling(2*max(data_bias_adj_melted$value)) / 2
#y_min <- floor(2*min(data_bias_adj_melted$value)) / 2
y_max <- ceiling(max(data_bias_adj_melted$value))
y_min <- min(floor(min(data_bias_adj_melted$value)), 0)

# plot
ggplot(data_bias_adj_melted, aes(x=position, y=value, colour=variable)) +
    geom_rect(alpha=0.2, size=0, color=NA, fill="grey95", aes(xmin=footprint_start, xmax=footprint_end, ymin=y_min, ymax=y_max)) +
    geom_line(alpha=0.7, size=0.230) + # 0.115
    labs(title=title_name, x="\nPosition (bp)", y="Flank adjusted enrichment\n") +
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
        axis.title.y=element_text(vjust=0, hjust=0.5, margin=margin(0,0,0,0)),
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
        legend.position=c(0.9, 0.9),
        strip.background=element_blank(),
        strip.text=element_blank()) +
    scale_color_manual(values=ggr_colors) +
    scale_x_continuous(limits=c(-x_lim, x_lim), expand=c(0,0)) +
    scale_y_continuous(limits=c(y_min, y_max), expand=c(0,0))

ggsave(plot_file, height=1.2, width=2, useDingbats=FALSE)


# ==================================================
# PLOT 4 - plot footprint differences
# ==================================================
plot_file <- paste(out_prefix, ".footprint_diff.pdf", sep="")

footprint_vals <- data_bias_adj[(data_bias_adj$position >= footprint_start) &
                                    (data_bias_adj$position <= footprint_end),]

# normalize to max point of footprint, then get average depth
footprint_max_vals <- apply(footprint_vals, 2, max)
data_footprints <- data.frame(t(t(footprint_vals) - footprint_max_vals))
footprint_avg_vals <- as.data.frame(t(apply(data_footprints, 2, mean)))
footprint_avg_vals$position <- NULL
footprint_avg_vals <- melt(footprint_avg_vals)

if (dim(data)[2] == 3) {
    footprint_avg_vals$variable <- gsub(".", "-", footprint_avg_vals$variable, fixed=TRUE)
    footprint_avg_vals$variable <- factor(footprint_avg_vals$variable, levels=c("NN-inactive", "NN-active"))
} else {
    #footprint_avg_vals$variable <- factor(footprint_avg_vals$variable, levels=c("d6.0", "d3.0", "d0.0"))
    footprint_avg_vals$variable <- factor(footprint_avg_vals$variable, levels=c("d0.0", "d3.0", "d6.0"))
}
print(footprint_avg_vals)

# pretty y scale
y_max <- max(ceiling(max(footprint_avg_vals$value)), 0)
y_min <- floor(min(footprint_avg_vals$value))

# plot
ggplot(footprint_avg_vals, aes(x=variable, y=value, colour=variable, fill=variable)) +
    geom_hline(yintercept=0, size=0.115, linetype="dashed") +
    geom_col(alpha=0.5, size=0.115, width=0.6, show.legend=FALSE) +
    labs(title=title_name, x="", y="Footprint average depth\n") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(hjust=0.5, size=8, margin=margin(b=1)),
        plot.margin=margin(5,0,-8,0),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(vjust=1, margin=margin(0,0,0,0)),
        axis.title.y=element_text(vjust=0, hjust=1, margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6, vjust=1, hjust=1, angle=60),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in")) +
    scale_color_manual(values=ggr_colors) +
    scale_fill_manual(values=ggr_colors) +
    scale_y_continuous(limits=c(y_min, y_max))

#ggsave(plot_file, height=1.2, width=0.75, useDingbats=FALSE)
ggsave(plot_file, height=1, width=0.75, useDingbats=FALSE)
