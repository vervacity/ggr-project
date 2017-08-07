#!/usr/bin/env Rscript

# Description: plot nice trajectories

library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)
timepoints_file <- args[1]
cluster_file <- args[2]
prefix <- args[3]

data <- read.table(timepoints_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
data_z <- data.frame(t(scale(t(data))))
rownames(data_z) <- rownames(data)
colnames(data_z) <- colnames(data)

clusters <- read.table(cluster_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)

data_z$region <- rownames(data_z)
colnames(clusters) <- c('cluster', 'region')

data_w_clusters <- merge(data_z, clusters, by='region')

# do a sort, subsample evenly and plot out ATAC in clusters
clustered_data <- data_w_clusters[order(data_w_clusters$cluster),]

clustered_data$cluster <- NULL
rownames(clustered_data) <- clustered_data$region
clustered_data$region <- NULL

evenly_spaced_subsample <- clustered_data[seq(1, nrow(clustered_data), 20), ]

# save out the subsample ids
subsample_id_file <- paste(prefix , ".subsampled.ordered.txt", sep="")
write.table(rownames(evenly_spaced_subsample), file=subsample_id_file, quote=FALSE, row.names=FALSE, col.names=FALSE)

# plotting
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
mylwid = c(2,6,2)
mylhei = c(0.5,12,1.5)

overall_subsampled_heatmap <- paste(prefix, ".clustered.heatmap.png", sep="")
png(overall_subsampled_heatmap, height=18, width=10, units="in", res=200)
heatmap.2(
    as.matrix(evenly_spaced_subsample),
    Rowv=FALSE,
    Colv=FALSE,
    dendrogram="none",
    trace='none',
    density.info="none",
    keysize=0.1,
    key.title=NA,
    key.xlab=NA,
    key.par=list(pin=c(4,0.1),
        mar=c(6.1,0,5.1,0),
        mgp=c(3,2,0),
        cex.axis=2.0,
        font.axis=2),
    srtCol=45,
    cexCol=3.0,
    labRow="",
    margins=c(3,0),
    col=my_palette,
    lmat=mylmat,
    lwid=mylwid,
    lhei=mylhei)
dev.off()


for (i in 1:length(unique(data_w_clusters$cluster))) {

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


