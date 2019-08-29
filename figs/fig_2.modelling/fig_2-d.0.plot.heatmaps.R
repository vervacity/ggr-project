#!/usr/bin/env Rscript

# description: plot heatmaps (predicted vs actual)
# data: /mnt/lab_data/kundaje/users/dskim89/ggr/nn/inference.2019-02-05

library(rhdf5)
library(gplots)
library(RColorBrewer)
library(reshape2)

# load GGR style guide
load_style_guide <- system("which ggr_style_guide.R", intern=TRUE)
source(load_style_guide)

# args
args <- commandArgs(trailingOnly=TRUE)
results_files <- args[1:length(args)]

# keys
actual_key <- "ATAC_SIGNALS.NORM"
predicted_key <- "logits.norm"
traj_key <- "TRAJ_LABELS"

# task indices of interest
task_indices <- c(1,2,3,4,5,6,7,10,11,13)

# for each file, load predicted and actual
for (file_idx in 1:length(results_files)) {
    filename <- results_files[file_idx]
    print(filename)

    # get actual and predicted
    file_actual <- t(
        h5read(filename, actual_key))[,task_indices]
    file_predicted <- t(
        h5read(filename, predicted_key))[,task_indices]
    file_traj <- t(
        h5read(filename, traj_key))
    
    # and concat
    if (file_idx == 1) {
        actual <- file_actual
        predicted <- file_predicted
        traj <- file_traj
    } else {
        actual <- rbind(actual, file_actual)
        predicted <- rbind(predicted, file_predicted)
        traj <- rbind(traj, file_traj)
    }
    
}

# now get a sort by the traj values
traj_sort_indices <- c(1,8,9,10,11,12,13,14,15,2,3,4,5,6)
traj <- traj[,traj_sort_indices]
traj_indexed <- apply(traj, 1, which.max)
example_order <- order(traj_indexed)

# sort all
actual <- actual[example_order,]
predicted <- predicted[example_order,]
traj <- traj[example_order,]
traj_indexed <- traj_indexed[example_order]

# normalize
actual <- actual - actual[,1]
predicted <- predicted - predicted[,1]

# percentile clip on actual
thresholds <- quantile(melt(actual)$value, c(0.05, 0.95))
actual[actual < thresholds[1]] <- thresholds[1]
actual[actual > thresholds[2]] <- thresholds[2]

# percentile clip on predicted
thresholds <- quantile(melt(predicted)$value, c(0.05, 0.95))
predicted[predicted < thresholds[1]] <- thresholds[1]
predicted[predicted > thresholds[2]] <- thresholds[2]

# subsample for plotting
max_n <- 1500
keep_indices <- seq(1, nrow(actual), length.out=max_n)
actual_subset <- actual[keep_indices,]
predicted_subset <- predicted[keep_indices,]
traj_indexed_subset <- traj_indexed[keep_indices]

# palettes
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49)) # 49
cluster_palette <- get_trajectory_palette(15)

# x labels
xlab <- c(
    "day 0.0",
    "day 0.5",
    "day 1.0",
    "day 1.5",
    "day 2.0",
    "day 2.5",
    "day 3.0",
    "day 4.5",
    "day 5.0",
    "day 6.0")

# figure out row seps and cluster colors
rowsep <- c(0)
cluster_ids <- c()
traj_val <- 1
for (i in 1:length(traj_indexed_subset)) {
    new_traj_val <- traj_indexed_subset[i]
    if (new_traj_val != traj_val) {
        rowsep <- c(rowsep, i)
        traj_val <- new_traj_val
    }
    cluster_ids <- c(cluster_ids, new_traj_val)
}
rowsep <- c(rowsep, i+1)
colsep <- c(0, 10)

# get cluster colors
cluster_colors <- cluster_palette[cluster_ids]

# plot
plot_heatmap <- function(data, plot_file, title) {

    mylmat = rbind(c(0,0,3,0),c(4,1,2,0),c(0,0,5,0))
    mylwid = c(0.05,0.1,0.7,0.3) # vs 0.7
    mylhei = c(0.10,2,0.27)
    
    pdf(plot_file, height=3.75, width=1.25, family="ArialMT", useDingbats=FALSE)
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        dendrogram="none",
        trace='none',
        density.info="none",
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(0.9,1,1.55,1),
            mgp=c(0,-0.1,0),
            tcl=-0.1,
            lend=2,
            cex.axis=0.6,
            bty="n"),
        key.xtickfun=function() {
            breaks=pretty(parent.frame()$breaks)
            breaks=breaks[c(1, as.integer(length(breaks)/2)+1, length(breaks))]
            list(at=parent.frame()$scale01(breaks),
                 labels = breaks)
        },
        srtCol=45,
        cexCol=0.5,
        offsetCol=-0.5,
        labRow="",
        labCol=xlab,
        margins=c(0,0),
        col=my_palette,
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        rowsep=rowsep,
        colsep=colsep,
        sepcolor="black",
        sepwidth=c(0.0001, 0.0001),
        RowSideColors=cluster_colors)
    title(title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
    dev.off()

}



# plot
plot_file <- "fig_2-d.0.actual.pdf"
plot_heatmap(actual_subset, plot_file, "Actual")

# plot
plot_file <- "fig_2-d.0.predicted.pdf"
plot_heatmap(predicted_subset, plot_file, "Predicted")
