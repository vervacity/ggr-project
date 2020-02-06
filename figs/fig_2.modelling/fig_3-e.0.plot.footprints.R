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

# params
null_flank_start <- 1 # how far from edge, start
null_flank_end <- 50 # how far from edge, end

acc_flank_start <- 100 # how far from edge, start
acc_flank_end <- 240 # how from edge, end. MAKE SURE DOES NOT OVERLAP FOOTPRINT

footprint_start <- -10
footprint_end <- 10

footprint_flank_length <- 25

# read in data, create index
data <- read.table(input_file, sep="\t", header=TRUE, row.names=1)
data$position <- as.numeric(rownames(data)) - (nrow(data) / 2)

print(head(data))

# normalizations: flank norm. this gives a fold change enrichment over background
l_null_flank_start <- data$position[null_flank_start]
l_null_flank_end <- data$position[null_flank_start] + null_flank_end
r_null_flank_start <- data$position[nrow(data)] - null_flank_start + 1
r_null_flank_end <- data$position[nrow(data)] - null_flank_end + 1

background_vals <- rbind(
    data[(data$position >= l_null_flank_start) & (data$position < l_null_flank_end),],
    data[(data$position <= r_null_flank_start) & (data$position > r_null_flank_end),])
background_avg_vals <- apply(background_vals, 2, mean)

data <- data.frame(t(t(data) / background_avg_vals))
data$position <- as.numeric(rownames(data)) - (nrow(data) / 2)

# melt
data_melted <- melt(data, id.vars=c("position"))

# ==================================================
# PLOT 1 - plot full pattern across all positions
# ==================================================
plot_file <- paste(out_prefix, ".full_length.pdf", sep="")

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

# make pretty y limit
y_max <- ceiling(2*max(abs(data_melted$value))) / 2
y_min <- floor(2*min(abs(data_melted$value))) / 2

# figure out title
title_name <- strsplit(basename(out_prefix), ".", fixed=TRUE)[[1]]
title_name <- gsub("HCLUST_\\d+.", "", title_name[2])

# plot
ggplot(data_melted, aes(x=position, y=value, colour=variable)) +
    geom_line(alpha=0.7, size=0.230) + # 0.115
    labs(title=title_name, x="Position (bp)", y="Fold Enrichment") +
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
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(limits=c(-y_min, y_max), expand=c(0,0)) #+

ggsave(plot_file, height=1.2, width=2, useDingbats=FALSE)

quit()

# ==================================================
# PLOT 2 - plot accessibility diff
# ==================================================
plot_file <- paste(out_prefix, ".accessibility_diff..pdf", sep="")

l_acc_flank_pos <- footprint_start - footprint_flank_length
r_acc_flank_pos <- footprint_end + footprint_flank_length

flank_vals <- rbind(
    data[(data$position > l_null_flank_pos) & (data$position < l_acc_flank_pos),],
    data[(data$position < r_null_flank_pos) & (data$position > r_acc_flank_pos),])
flank_avg_vals <- apply(flank_vals, 2, mean)

quit()

# ==================================================
# PLOT 3 - plot full pattern across all positions
# ==================================================
plot_file <- paste(out_prefix, ".footprint_w_flanks..pdf", sep="")

# bias adjust

data <- data.frame(t(t(data) / background_avg_vals))
data$position <- as.numeric(rownames(data)) - (nrow(data) / 2)



data_norm <- rbind(data[50:75,], data[125:150,])
norm_factors <- apply(data_norm, 2, mean)
data <- data.frame(t(t(data) - norm_factors))

data$index <- as.numeric(rownames(data)) - (nrow(data) / 2)


# PLOT 2 - calculate accessibility enrichment and plot out bar plot


quit()



# normalize to flank
data_norm <- rbind(data[1:10,], data[191:200,])
flank_norm_factors <- apply(data_norm, 2, mean)

print(head(data))

# calculate flank accessibility
#print(head(data_melted))
#pos_l_flank <- data_melted[(data_melted$variable == "pos") & (data_melted$index < footprint_start),]
#pos_r_flank <- data_melted[(data_melted$variable == "pos") & (data_melted$index > footprint_end),]


#data_norm <- rbind(data[50:75,], data[125:150,])
#center_norm_factors <- apply(data_norm, 2, mean)

#norm_factors <- center_norm_factors - flank_norm_factors
data <- data.frame(t(t(data) / flank_norm_factors))

# and then bias adjust
data_norm <- rbind(data[50:75,], data[125:150,])
norm_factors <- apply(data_norm, 2, mean)
data <- data.frame(t(t(data) - norm_factors))

data$index <- as.numeric(rownames(data)) - (nrow(data) / 2)

# clipping
data <- data[50:150,]

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

# make pretty y limit
y_limit <- ceiling(2*max(abs(data_melted$value))) / 2

# figure out title
title_name <- strsplit(basename(plot_file), ".", fixed=TRUE)[[1]]
title_name <- title_name[2]

# plot
ggplot(data_melted, aes(x=index, y=value, colour=variable)) +
    geom_line(alpha=0.7, size=0.230) + # 0.115
    labs(title=title_name, x="Position (bp)", y="Normalized read coverage") +
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
    scale_x_continuous(limits=c(-50,50), expand=c(0,0)) +
    scale_y_continuous(limits=c(-y_limit, y_limit), expand=c(0,0)) #+
    #facet_grid(. ~ variable)

ggsave(plot_file, height=1.2, width=2, useDingbats=FALSE)



# analyze footprint depth

