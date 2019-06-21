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
plot_file <- args[1] # "fig_1-c.atac.pdf"
mat_files <- args[2:length(args)]

# load data
for (i in 1:length(mat_files)) {
    print(mat_files[i])
    
    data <- read.table(gzfile(mat_files[i]), sep="\t", header=TRUE, row.names=1)

    if (i == 1) {
        all_data <- data
    } else {
        all_data <- cbind(all_data, data)
    }
    
}

#all_data <- all_data[,order(names(all_data))]

# just keep the top percentile?
if (FALSE) {
    all_data_max <- apply(all_data, 1, max)
    cutoff <- quantile(all_data_max, probs=c(0.50))
    print(cutoff[1])
    all_data <- all_data[all_data_max > cutoff,]
}

print(head(all_data))
print(dim(all_data))

# put into pca
pca_obj <- prcomp(t(all_data), center=TRUE, scale=TRUE)
pc1 <- pca_obj$x[,1]
pc2 <- pca_obj$x[,2]

pca_data <- data.frame(x=pc1, y=pc2)

# adjust group names
group_levels <- colnames(all_data)
group_levels <- gsub("0_", ".0 ", group_levels)
group_levels <- gsub("5_", ".5 ", group_levels)
print(group_levels)
pca_data$group <- row.names(pca_data)
pca_data$group <- gsub("0_", ".0 ", pca_data$group)
pca_data$group <- gsub("5_", ".5 ", pca_data$group)
print(head(pca_data))
pca_data$group <- factor(pca_data$group, levels=group_levels)

# pull GGR colors and adjust
my_colors <- get_ggr_timepoint_colors()
if (grepl("rna", plot_file, fixed=TRUE)) {
    my_colors <- c(my_colors[1], my_colors[3:length(my_colors)])
}

my_colors_r1 <- my_colors
my_colors_r1 <- lighten(my_colors, amount=0.2)
my_colors_r2 <- darken(my_colors, amount=0.2)
my_colors <- c(my_colors_r1, my_colors_r2)
print(head(pca_data))

# plot
p <- ggplot(pca_data, aes(x=x, y=y, colour=group)) +
    geom_point(size=0.25) +
    labs(x="PC1", y="PC2") + 
    scale_color_manual(values=my_colors) + 
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        title=element_text(size=6),
        aspect.ratio=1,
        plot.margin=margin(0,0,0,0),
        plot.title=element_text(margin=margin(0,0,0,0)),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=0.25),
        axis.title=element_text(size=5, margin=margin(0,0,0,0)),
        axis.line=element_blank(),
        axis.text=element_text(size=4),
        axis.ticks=element_line(size=0.25),
        axis.ticks.length=unit(0.01, "in"),
        legend.text=element_text(size=1),
        legend.title=element_blank(),
        legend.key.size = unit(0, 'in'))

if (!grepl("rna", plot_file, fixed=TRUE)) {
    # adjust limits of plot as needed
    p <- p + xlim(-250,250) + ylim(-250,250)
    p <- p + ggtitle("ATAC-seq")
} else {
    p <- p + ggtitle("RNA-seq")
}

ggsave(plot_file, height=1, width=1.5, units="in", useDingbats=FALSE)
