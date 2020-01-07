#!/usr/bin/env Rscript

# description: plot simple PCA
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/results/atac/timeseries/plots
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/integrative/v1.0.0a/results/rna/timeseries/matrices
library(ggplot2)
library(colorspace)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
plot_file <- args[1]
mat_files <- args[2:length(args)]

# load data
for (i in 1:length(mat_files)) {
    data <- read.table(gzfile(mat_files[i]), sep="\t", header=TRUE, row.names=1)
    if (i == 1) {
        all_data <- data
        filt_data <- data
    } else {
        all_data <- cbind(all_data, data)
        filt_data <- as.matrix(filt_data) + as.matrix(data)
    }
    
}

# just keep the top percentile?
if (grepl("atac", plot_file, fixed=TRUE)) {
    all_data_max <- apply(filt_data, 1, sd) / apply(filt_data, 1, mean)
    cutoff <- quantile(all_data_max, probs=c(0.66))
    all_data <- all_data[all_data_max > cutoff,]
    print(dim(all_data))
}

# put into pca
pca_obj <- prcomp(t(all_data), center=TRUE, scale=TRUE)
pc1 <- pca_obj$x[,1]
pc2 <- pca_obj$x[,2]
pca_data <- data.frame(x=pc1, y=pc2)

# adjust group names
pca_data$group <- row.names(pca_data)
group_levels <- colnames(all_data)

print(pca_data)

if (!grepl("histone", plot_file, fixed=TRUE) & !grepl("hichip", plot_file, fixed=TRUE)) {
    group_levels <- gsub("0_", ".0_", group_levels)
    group_levels <- gsub("5_", ".5_", group_levels)

    pca_data$group <- gsub("0_", ".0_", pca_data$group)
    pca_data$group <- gsub("5_", ".5_", pca_data$group)
}

pca_data$group <- factor(pca_data$group, levels=group_levels)
pca_data$day <- factor(gsub("_b.", "", pca_data$group))
pca_data$rep <- as.character(gsub(".+_", "", pca_data$group))

# pull GGR colors and adjust
my_colors <- get_ggr_timepoint_colors()
if (grepl("rna", plot_file, fixed=TRUE)) {
    # drop d0.5
    my_colors <- c(my_colors[1], my_colors[3:length(my_colors)])
    title <- "RNA-seq"
}
if (grepl("atac", plot_file, fixed=TRUE)) {
    title <- "ATAC-seq"
}
if (grepl("histone", plot_file, fixed=TRUE)) {
    # only 3 colors
    my_colors <- c(my_colors[1], my_colors[7], my_colors[10])
    title <- "ChIP-seq"
}
if (grepl("hichip", plot_file, fixed=TRUE)) {
    # only 3 colors
    my_colors <- c(my_colors[1], my_colors[7], my_colors[10])
    title <- "HiChIP"
}

# adjust shades for different reps
my_colors_r1 <- my_colors
my_colors_r1 <- lighten(my_colors, amount=0.2)
my_colors_r2 <- darken(my_colors, amount=0.2)
if (grepl("hichip", plot_file, fixed=TRUE)) {
    my_colors_joint <- c(rbind(my_colors_r1, my_colors_r2))
} else {
    my_colors_joint <- c(my_colors_r1, my_colors_r2)
}

# plot
p <- ggplot(pca_data, aes(x=x, y=y, colour=group, fill=day)) +
    geom_tile() + 
    geom_tile(size=1.1, fill="white", colour="white", show.legend=FALSE) +
    #geom_point(size=1, show.legend=FALSE) + # size=0.25
    geom_point(shape=21, fill="white", size=1.25, stroke=0.5, show.legend=FALSE) +
    labs(x="PC1", y="PC2", title=title) + 
    scale_color_manual(values=my_colors_joint, guide="none") +
    scale_fill_manual(values=my_colors) +
    coord_fixed() +
    theme_bw() +
    theme(
        #aspect.ratio=1,
        text=element_text(family="ArialMT"),
        plot.margin=margin(5,20,1,0),
        plot.title=element_text(size=8, margin=margin(b=1)),
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
        legend.position=c(1.1, 0.5),
        legend.text=element_text(size=5),
        legend.title=element_blank(),
        legend.spacing.x=unit(0.05, "in"),
        legend.key=element_rect(size=0.1, fill="white", colour="black"),
        legend.key.size=unit(0.05, 'in'))

# default scale (RNA)
scale_lim <- 110

if (grepl("atac", plot_file, fixed=TRUE)) {
    scale_lim <- 130
}

if (grepl("H3K27ac", plot_file, fixed=TRUE)) {
    scale_lim <- 230
}

if (grepl("H3K4me1", plot_file, fixed=TRUE)) {
    scale_lim <- 200
}

if (grepl("H3K27me3", plot_file, fixed=TRUE)) {
    scale_lim <- 60
}

if (grepl("hichip", plot_file, fixed=TRUE)) {
    scale_lim <- 400
}

y_lim <- 0.75*scale_lim

p <- p + scale_x_continuous(limits=c(-scale_lim, scale_lim), expand=c(0,0)) +
    scale_y_continuous(limits=c(-y_lim, y_lim), expand=c(0,0))


ggsave(plot_file, height=1.75, width=2.5, useDingbats=FALSE)
