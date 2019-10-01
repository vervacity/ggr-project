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
    } else {
        all_data <- cbind(all_data, data)
    }
    
}

# just keep the top percentile?
if (FALSE) {
    all_data_max <- apply(all_data, 1, max)
    cutoff <- quantile(all_data_max, probs=c(0.50))
    print(cutoff[1])
    all_data <- all_data[all_data_max > cutoff,]
}

# put into pca
pca_obj <- prcomp(t(all_data), center=TRUE, scale=TRUE)
pc1 <- pca_obj$x[,1]
pc2 <- pca_obj$x[,2]
pca_data <- data.frame(x=pc1, y=pc2)

# adjust group names
group_levels <- colnames(all_data)
group_levels <- gsub("0_", ".0 ", group_levels)
group_levels <- gsub("5_", ".5 ", group_levels)

pca_data$group <- row.names(pca_data)
pca_data$group <- gsub("0_", ".0 ", pca_data$group)
pca_data$group <- gsub("5_", ".5 ", pca_data$group)

pca_data$group <- factor(pca_data$group, levels=group_levels)
pca_data$day <- factor(gsub(" b.", "", pca_data$group))

pca_data$rep <- as.character(gsub(".+ ", "", pca_data$group))

# pull GGR colors and adjust
my_colors <- get_ggr_timepoint_colors()
if (grepl("rna", plot_file, fixed=TRUE)) {
    my_colors <- c(my_colors[1], my_colors[3:length(my_colors)])
    mid_idx <- 10
    title <- "RNA-seq"
} else {
    mid_idx <- 9
    title <- "ATAC-seq"
}
my_colors_r1 <- my_colors
my_colors_r1 <- lighten(my_colors, amount=0.2)
my_colors_r2 <- darken(my_colors, amount=0.2)
my_colors_joint <- c(my_colors_r1, my_colors_r2)

# plot
p <- ggplot(pca_data, aes(x=x, y=y, colour=group, fill=day)) +
    geom_tile() + 
    geom_tile(fill="white", show.legend=FALSE) +
    geom_point(size=0.50, show.legend=FALSE) + # size=0.25
    labs(x="PC1", y="PC2", title=title) + 
    scale_color_manual(values=my_colors_joint, guide="none") +
    scale_fill_manual(values=my_colors) +
    theme_bw() +
    theme(
        aspect.ratio=1,
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

if (!grepl("rna", plot_file, fixed=TRUE)) {
    # adjust limits of plot as needed
    p <- p + scale_x_continuous(limits=c(-300,300), breaks=c(-200, 0, 200), expand=c(0,0)) +
        scale_y_continuous(limits=c(-300,300), breaks=c(-200, 0, 200), expand=c(0,0))
} else {
    p <- p + scale_x_continuous(limits=c(-90,90), expand=c(0,0)) +
        scale_y_continuous(limits=c(-90,90), expand=c(0,0))
}

ggsave(plot_file, height=2, width=2.5, useDingbats=FALSE)
