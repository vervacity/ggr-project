#!/usr/bin/env Rscript

library(fastcluster)
library(gplots)
library(RColorBrewer)
library(dendextend)

set.seed(1337)

# Given a cluster file and data file, organizes
# the clusters by hclust on the cluster means
# to be able to get cleaner heatmap outputs

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
cluster_files <- args[2:length(args)]

# set up cluster means table
cluster_means <- data.frame()
cluster_sizes <- c()

# for each cluster file
for (i in 1:length(cluster_files)) {

    # read in cluster file. data is already z-scored
    data <- read.table(gzfile(cluster_files[i]), header=FALSE, row.names=1, sep='\t')

    # get mean and save to means dataframe
    cluster_mean <- colMeans(data)
    cluster_means <- rbind(cluster_means, cluster_mean)

    cluster_sizes <- c(cluster_sizes, nrow(data))
    
}

# calculate a dist matrix
my_hclust <- function(d) {
    return(hclust(d, method="ward.D2"))

}

hclust_means <- hclust(dist(as.matrix(cluster_means)), method="ward.D2", members=cluster_sizes)
#hclust_means <- hclust(dist(as.matrix(cluster_means)), method="single", members=cluster_sizes)

cluster_order <- hclust_means$order

# now copy and rename all cluster files with new names
for (i in 1:length(cluster_order)) {

    old_cluster_idx <- cluster_order[i]
    
    prefix <- unlist(strsplit(basename(cluster_files[old_cluster_idx]), "cluster", fixed=TRUE))[1]
    suffix <- unlist(strsplit(basename(cluster_files[old_cluster_idx]), "cluster", fixed=TRUE))[2]

    new_filename <- paste(out_dir, "/", prefix, "newcluster_", i, ".cluster", suffix, sep="")
    print(new_filename)

    copy_file <- paste("cp", cluster_files[old_cluster_idx], new_filename)
    print(copy_file)
    system(copy_file)
    
}






