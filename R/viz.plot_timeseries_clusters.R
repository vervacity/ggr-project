#!/usr/bin/env Rscript

# Description: plot nice trajectories

library(ggplot2)
library(reshape2)

library(grid)
library(gridGraphics)
library(gridExtra)

# args
args <- commandArgs(trailingOnly=TRUE)
cluster_file <- args[1]
mat_file <- args[2]
out_dir <- args[3]
prefix <- args[4]

ci_top <- function(vals) {
    val_mean <- mean(vals)
    val_se <- sd(vals) # / sqrt(length(vals))
    return(val_mean+2*val_se)
}

ci_bottom <- function(vals) {
    val_mean <- mean(vals)
    val_se <- sd(vals) #/ sqrt(length(vals))
    return(val_mean-2*val_se)
}

# read in files
clusters <- read.table(cluster_file, header=TRUE)
colnames(clusters)[2] <- "id"
data <- read.table(gzfile(mat_file), header=TRUE)
data$id <- rownames(data)

# merge
data_w_clusters <- merge(clusters, data, by="id", sort=FALSE)
rownames(data_w_clusters) <- data_w_clusters$id
data_w_clusters$id <- NULL

# now plot by each cluster
all_plots <- list()
cluster_names <- unique(clusters$cluster)
for (i in 1:length(cluster_names)) {
    
    # extract the cluster mat
    single_cluster_data <- data_w_clusters[data_w_clusters$cluster==cluster_names[i],]
    
    # remove the cluster id and rename columns
    single_cluster_data$cluster <- NULL
    if (ncol(single_cluster_data) == 10) {
        colnames(single_cluster_data) <- c("d0.0", "d0.5", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0")
    } else if (ncol(single_cluster_data) == 9) {
        colnames(single_cluster_data) <- c("d0.0", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0")
    }

    # normalize
    #data_z <- t(scale(t(single_cluster_data), center=TRUE, scale=TRUE))
    data_z <- as.matrix(single_cluster_data - single_cluster_data[,1])
    
    # adjustments for atac vs rna
    if (grepl("epigenome", prefix)) {
        method <- "auto"
        limit_val <- 2
    } else {
        method <- "loess"
        limit_val <- 4
    }
    data_z[data_z > limit_val] <- limit_val
    data_z[data_z < -limit_val] <- -limit_val
    
    # melt data
    data_melted <- melt(data_z)
    colnames(data_melted)[1:2] <- c("id", "timepoint")
    data_melted$timepoint <- as.numeric(sub("d", "", data_melted$timepoint))
    
    # generate plot
    plot_file <- paste(out_dir, "/", prefix, ".cluster_", cluster_names[i], ".plot.pdf", sep="")
    if (TRUE) {
    #if (!file.exists(plot_file)) {

        p <- ggplot(data_melted, aes(x=timepoint, y=value, group=id)) +
            geom_smooth(
                aes(x=timepoint, y=value, group=1),
                colour="gray72", fill="gray72", method=method, se=TRUE, level=0.95) +
            geom_smooth(
                aes(x=timepoint, y=value, group=1),
                size=0.230,
                colour="black", method=method, se=FALSE) +
            theme_bw() +
            theme(
                text=element_text(family="ArialMT"),
                plot.title=element_text(size=6, margin=margin(b=0)),
                plot.margin=margin(1,5,1,5),
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                axis.title=element_text(size=6),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.y=element_text(size=5),
                axis.text.x=element_blank(),
                axis.line.x=element_blank(),
                axis.line.y=element_line(color="black", size=0.115, lineend="square"),
                axis.ticks=element_line(size=0.115),
                axis.ticks.length=unit(0.01, "in"),
                axis.ticks.x=element_blank()) +
            scale_x_continuous(limits=c(0,6), expand=c(0,0)) + 
            scale_y_continuous(
                limits=c(-limit_val,limit_val), expand=c(0,0), breaks=seq(-limit_val,limit_val,by=limit_val),
                labels=c(-limit_val, "", limit_val)) +
            geom_hline(size=0.115, aes(yintercept=0))
            #stat_summary(
            #    aes(y=value, group=1),
            #    fun.y=mean,
            #    colour="black", geom="line", size=0.230, group=1)
        
        # attach to all plot list
        all_plots[[i]] <- p

        # and then add some more adjustments for saving out to individual file
        p <- p + labs(
            title=paste("Cluster", cluster_names[i]),
            x="Time (days)",
            y="Z-score") +
        ggsave(plot_file, width=0.5, height=0.25)
    }


}

# and plot all
plot_file <- paste(out_dir, "/", prefix, ".clusters-ALL.pdf", sep="")
pdf(plot_file, height=2.75, width=0.5, onefile=FALSE, family="ArialMT")
grid.newpage()
grid.arrange(grobs=all_plots, ncol=1, nrow=length(cluster_names))
dev.off()

