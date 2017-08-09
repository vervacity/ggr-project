#!/usr/bin/env Rscript

library(fastcluster)
library(gplots)
library(RColorBrewer)

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

# for each cluster,
for (i in 1:length(unique(clusters$cluster))) {

    cluster_num <- i

    cluster_subset <- data_w_clusters[data_w_clusters$cluster == cluster_num, ]
    cluster_subset$cluster <- NULL
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















