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
cluster_file <- args[1]
data_file <- args[2]
new_cluster_file <- args[3]

# read/clean clusters
clusters <- read.table(gzfile(cluster_file), header=TRUE, sep='\t')
colnames(clusters) <- c("cluster", "region")

# read/clean data
data <- read.table(gzfile(data_file), header=TRUE, row.names=1, sep='\t')
data_z <- data.frame(t(scale(t(data))))
data_z$region <- rownames(data_z)

# merge
data_w_clusters <- merge(data_z, clusters, by="region")
rownames(data_w_clusters) <- data_w_clusters$region
data_w_clusters$region <- NULL

# set up cluster means table
cluster_means <- data.frame()
cluster_sizes <- c()

# for each cluster,
for (i in 1:length(unique(clusters$cluster))) {

    cluster_num <- i

    cluster_subset <- data_w_clusters[data_w_clusters$cluster == cluster_num, ]
    cluster_subset$cluster <- NULL

    cluster_sizes <- c(cluster_sizes, nrow(cluster_subset))
    
    cluster_mean <- colMeans(cluster_subset)
    cluster_means <- rbind(cluster_means, cluster_mean)

}

colnames(cluster_means) <- colnames(cluster_subset)

# calculate a dist matrix
my_hclust <- function(d) {
    return(hclust(d, method="ward.D2"))

}

hclust_means <- my_hclust(dist(as.matrix(cluster_means)))
new_cluster_nums <- hclust_means$order

# ===================================
# for dendrogram

cluster_sizes <- cluster_sizes

print(cluster_sizes)

ratios <- as.integer(cluster_sizes / min(cluster_sizes))

# set up repeated clusters
cluster_means_rep <- data.frame()
for (cluster_idx in 1:length(ratios)) {
    for (i in 1:ratios[cluster_idx]) {
        cluster_means_rep <- rbind(cluster_means_rep, cluster_means[cluster_idx,])
    }
}

# hclust
hclust_dendro <- my_hclust(dist(as.matrix(cluster_means_rep)))
dend <- as.dendrogram(hclust_dendro)

png("test.dendro.0.png")
plot(dend)
dev.off()

q()





hclust_means_double <- my_hclust(dist(as.matrix(clsuter_means_double)))
dend <- as.dendrogram(hclust_means_double)


get_spacing <- function(dendro, cluster_sizes) {
    
    if (is.null(attr(dendro, "leaf"))) {
        midpoint_val <- attr(dendro, "midpoint")
    } else {
        # it's a leaf
        midpoint_val <- cluster_sizes[as.numeric(attr(dendro, "label"))] / 4.0
    }
    return(midpoint_val)
    
}


adjust_midpoint_at_height <- function(dendro, cluster_sizes, height) {
    
    # TODO reorder
    if (!is.null(attr(dendro, "leaf"))) {
        # don't do anything
    } else if ( attr(dendro, "height") == height ) {
        # adjust the midpoint here
        left_space <- get_spacing(dendro[[1]], cluster_sizes)
        right_space <- get_spacing(dendro[[2]], cluster_sizes)
        #print(attr(dendro, "midpoint"))
        attr(dendro, "midpoint") <- (left_space + right_space) / 2.0
        #print(attr(dendro, "midpoint"))
    } else {
        # not a leaf but not the right height, so keep doing down
        dendro[[1]] <- adjust_midpoint_at_height(dendro[[1]], cluster_sizes, height)
        dendro[[2]] <- adjust_midpoint_at_height(dendro[[2]], cluster_sizes, height)
    }
    return(dendro)
}


adjust_dendro <- function(dendro, cluster_sizes, node_heights) {

    # for each height level, starting from lowest (minus leaves) going up, adjust midpoints
    for (height_idx in 2:length(node_heights)) {
        print(node_heights[height_idx])
        dendro <- adjust_midpoint_at_height(
            dendro, cluster_sizes, node_heights[height_idx])
    }
    
    return(dendro)

}


unique_heights <- sort(unique(get_nodes_attr(dend, "height")))
print(unique_heights)

print(get_nodes_attr(dend, "midpoint"))

dendro_new <- adjust_dendro(dend, cluster_sizes, unique_heights)

print(get_nodes_attr(dendro_new, "midpoint"))

png("test.dendro.1.png")
plot(dend)
dev.off()


png("test.dendro.2.png")
plot(dendro_new)
dev.off()



q()




# END DENDRO


print(new_cluster_nums)

# change old numbers for new ones
data_w_clusters$new_cluster <- match(data_w_clusters$cluster, hclust_means$order)
data_w_clusters_sorted <- data_w_clusters[order(data_w_clusters$new_cluster),]


print(head(data_w_clusters_sorted))

# remove old clusters and save out
new_cluster_data <- data.frame(
    regions=rownames(data_w_clusters_sorted),
    cluster=data_w_clusters_sorted$new_cluster)

print(head(new_cluster_data))

write.table(new_cluster_data, "test.txt", sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

# plot out as sanity check
data_w_clusters_sorted$cluster <- NULL
data_w_clusters_sorted$new_cluster <- NULL

evenly_spaced_subsample <- data_w_clusters_sorted[seq(1, nrow(data_w_clusters_sorted), 20), ]
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
mylwid = c(2,6,2)
mylhei = c(0.5,12,1.5)

#overall_subsampled_heatmap <- paste(prefix, ".clustered.heatmap.png", sep="")
png("testing.png", height=18, width=10, units="in", res=200)
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















