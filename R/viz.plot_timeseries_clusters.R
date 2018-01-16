#!/usr/bin/env Rscript

# Description: plot nice trajectories

library(ggplot2)
library(reshape)

# args
args <- commandArgs(trailingOnly=TRUE)
cluster_file <- args[1]
mat_file <- args[2]
out_dir <- args[3]
prefix <- args[4]

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
cluster_names <- unique(clusters$cluster)
for (i in 1:length(cluster_names)) {
    
    # extract the cluster mat
    single_cluster_data <- data_w_clusters[data_w_clusters$cluster==cluster_names[i],]
    #print(dim(single_cluster_data))
    
    # remove the cluster id and rename columns
    single_cluster_data$cluster <- NULL
    colnames(single_cluster_data) <- c("d0.0", "d0.5", "d1.0", "d1.5", "d2.0", "d2.5", "d3.0", "d4.5", "d5.0", "d6.0")

    # normalize
    # TODO produce a fold change from baseline version of this: log(mat / val) in first column
    data_z <- t(scale(t(single_cluster_data), center=TRUE, scale=TRUE))
    #print(head(data_z))
    
    # melt data
    data_melted <- melt(data_z)
    colnames(data_melted)[1:2] <- c("id", "timepoint")
    data_melted$timepoint <- as.numeric(sub("d", "", data_melted$timepoint))

    # generate plot
    plot_file <- paste(out_dir, "/", prefix, ".cluster_", cluster_names[i], ".plot.pdf", sep="")
    if (TRUE) {
    #if (!file.exists(plot_file)) {

        # TODO - consider lighter gray SE? to overlay
        ggplot(data_melted, aes(x=timepoint, y=value, group=id)) + geom_smooth(colour="gray72", fill="gray72", se=TRUE, level=0.68) + #geom_line(colour="gray44") +
            theme_classic(base_size=22) + 
            #theme_bw() +
            theme(
                axis.line=element_line(colour="black"),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                axis.line.x=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
            #ggtitle(paste("Cluster", cluster_names[i])) + 
            theme(plot.title=element_text(hjust=0.5)) + 
            #xlab("Time (days)") + 
            #ylab("Z-score") + 
            scale_x_continuous(limits=c(0,6), expand=c(0,0)) + 
            scale_y_continuous(limits=c(-3,3), expand=c(0,0), breaks=seq(-3,3,by=1), labels=c("-3", "", "", "0", "", "", "3")) +
            geom_smooth(aes(y=value, group=1), colour="black") +
            geom_hline(aes(yintercept=0))
                                        #+ 
            #stat_summary(aes(y=value, group=1), fun.y=mean, colour="black", geom="line", group=1)
                
        ggsave(plot_file, width=3.5, height=2)

    }
    
}
