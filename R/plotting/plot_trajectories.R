#!/usr/bin/env Rscript

# Description: plot nice trajectories


library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly=TRUE)




data_file <- 'ggr.atac.idr.master.distal.dynamic.nomedia.mat.txt'
cluster_file <- 'out/ggr.atac_optimal_clustering.txt'
prefix <- 'ggr.atac'

data <- read.table(data_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
data_z <- data.frame(t(scale(t(data))))
rownames(data_z) <- rownames(data)
colnames(data_z) <- colnames(data)

clusters <- read.table(cluster_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)

data_z$region <- rownames(data_z)
colnames(clusters) <- c('cluster', 'region')

data_w_clusters <- merge(data_z, clusters, by='region')
print(head(data_w_clusters))

for (i in 1:14) {

    # get subset of data
    data_subset <- data_w_clusters[data_w_clusters$cluster == i,]
    data_subset$cluster <- NULL
    print(head(data_subset))
    
    data_melted <- melt(data_subset, id='region')
    data_melted$variable <- sub('X', '', data_melted$variable)
    data_melted$variable <- as.numeric(as.character(data_melted$variable))
    print(head(data_melted))
    
    # plot
    group_plot <- paste(prefix, '.dynamic.dp_gp.cluster_', i, '.png', sep='')
    if (!file.exists(group_plot)) {
        ggplot(data_melted, aes(x=variable, y=value, group=region)) + geom_line(colour='gray44') +
            theme_bw() +
            theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                                        #panel.border = element_blank(),
                  panel.background = element_blank()) +
            ggtitle(paste('RNA expression pattern ', i, sep='')) +
            theme(plot.title=element_text(hjust = 0.5)) +
            xlab('Time (days)') +
            ylab('Expression (z score)') +
            scale_x_continuous(limits = c(0, 6), expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))
        ggsave(group_plot)
    }

    
}


