#!/usr/bin/env Rscript

library(fastcluster)
library(gplots)
library(RColorBrewer)
library(dendextend)

# for DIANA
library(cluster)

set.seed(1337)

# Given a cluster file and data file, organizes
# the clusters by hclust on the cluster means
# to be able to get cleaner heatmap outputs

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]
hard_cluster_file <- args[2]
soft_cluster_files <- args[3:length(args)]

#agglom_method <- "single"
agglom_method <- "complete"

# minor tweaks to ordering for better viz
tweak_diana <- TRUE

# set up cluster means table
cluster_means <- data.frame()
cluster_sizes <- c()
cluster_nums <- c()

# for each cluster file
for (i in 1:length(soft_cluster_files)) {

    # read in cluster file. data is already z-scored
    data <- read.table(gzfile(soft_cluster_files[i]), header=FALSE, row.names=1, sep='\t')

    # get mean and save to means dataframe
    cluster_mean <- colMeans(data)
    cluster_means <- rbind(cluster_means, cluster_mean)

    cluster_sizes <- c(cluster_sizes, nrow(data)) # TODO need to get cluster sizes from hard clusters

    # also extract cluster number from filename
    suffix_tmp <- unlist(strsplit(basename(soft_cluster_files[i]), ".cluster_", fixed=TRUE))[2]
    cluster_num <- unlist(strsplit(suffix_tmp, ".soft", fixed=TRUE))[1]
    cluster_nums <- c(cluster_nums, as.numeric(cluster_num))
    
}

#hclust_means <- hclust(dist(as.matrix(cluster_means)), method="ward.D2", members=cluster_sizes)
#hclust_means <- hclust(dist(as.matrix(cluster_means)), method=agglom_method, members=cluster_sizes)
#cluster_order <- hclust_means$order

# testing divisive agglom
dclust <- as.hclust(
    diana(dist(as.matrix(cluster_means))))
cluster_order <- dclust$order


if (tweak_diana) {

    # first reverse things
    cluster_order[1:6] <- rev(cluster_order[1:6])
    cluster_order[7:10] <- rev(cluster_order[7:10])
    cluster_order[8:11] <- rev(cluster_order[8:11])
    
}

# now copy and rename all cluster files with new names
for (i in 1:length(cluster_order)) {

    old_cluster_idx <- cluster_order[i]
    
    prefix <- unlist(strsplit(basename(soft_cluster_files[old_cluster_idx]), ".cluster", fixed=TRUE))[1]
    suffix <- unlist(strsplit(basename(soft_cluster_files[old_cluster_idx]), ".cluster", fixed=TRUE))[2]

    new_filename <- paste(out_dir, "/", prefix, ".newcluster_", i, ".cluster", suffix, sep="")
    print(new_filename)

    copy_file <- paste("cp", soft_cluster_files[old_cluster_idx], new_filename)
    #print(copy_file)
    system(copy_file)
    
}

# and for the hard cluster files, reorder
hard_clusters <- read.table(gzfile(hard_cluster_file), header=FALSE, sep='\t')
colnames(hard_clusters)[ncol(hard_clusters)] <- "cluster"
hard_clusters$new_clusters_tmp <- match(hard_clusters$cluster, cluster_nums)
hard_clusters$new_clusters <- match(hard_clusters$new_clusters_tmp, cluster_order)

# and also get cluster sizes
cluster_sizes <- c()
for (i in 1:length(cluster_nums)) {

    hard_cluster_data_tmp <- hard_clusters[hard_clusters$cluster == cluster_nums[i],]
    cluster_sizes <- c(cluster_sizes, nrow(hard_cluster_data_tmp))
    
}

hard_clusters$new_clusters_tmp <- NULL
hard_clusters$cluster <- NULL
hard_clusters <- hard_clusters[with(hard_clusters, order(new_clusters)), ]

prefix <- unlist(strsplit(basename(hard_cluster_file), ".hard", fixed=TRUE))[1]
suffix <- unlist(strsplit(basename(hard_cluster_file), ".hard", fixed=TRUE))[2]
new_filename <- paste(out_dir, "/", prefix, ".hard.renumbered", suffix, sep="")
print(new_filename)
write.table(hard_clusters, file=gzfile(new_filename), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

# and plot out for sanity check
rownames(hard_clusters) <- hard_clusters$V1
hard_clusters$V1 <- NULL

# get subsample
evenly_spaced_subsample <- hard_clusters[seq(1, nrow(hard_clusters), 20), ]

# determine row sep points
rowsep <- c()
cluster <- 1
new_clusters <- evenly_spaced_subsample$new_clusters
evenly_spaced_subsample$new_clusters <- NULL
for (i in 1:nrow(evenly_spaced_subsample)) {
    
    if (new_clusters[i] != cluster) {
        rowsep <- c(rowsep, i)
        cluster <- new_clusters[i]
    }

}

# plotting
plot_file <- paste(out_dir, "/", prefix, ".hard.sanity_check.heatmap.png", sep="")
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
mylwid = c(2,6,2)
mylhei = c(0.5,12,1.5)

png(plot_file, height=18, width=6, units="in", res=200)
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
    lhei=mylhei,
    rowsep=rowsep,
    sepcolor="black")
dev.off()


# dendrogram
ratios <- as.integer(50 * cluster_sizes / min(cluster_sizes)) # increase by factor of 10 to even out better

# set up repeated clusters
cluster_means_rep <- data.frame()
cluster_means_group <- c()
for (cluster_idx in 1:length(ratios)) {
    for (i in 1:ratios[cluster_idx]) {
        cluster_means_rep <- rbind(cluster_means_rep, jitter(as.numeric(cluster_means[cluster_idx,])))
        #cluster_means_group <- c(cluster_means_group, cluster_order[cluster_idx])
        cluster_means_group <- c(cluster_means_group, cluster_idx)
    }
    #print(cluster_order[cluster_idx])
}

# hclust
#hclust_dendro <- hclust(dist(as.matrix(cluster_means_rep)), method=agglom_method)
#hclust_dendro <- my_hclust(dist(as.matrix(cluster_means_rep)))
#dend <- as.dendrogram(hclust_dendro)

dclust_dendro <- as.hclust(
    diana(dist(as.matrix(cluster_means_rep))))
dend <- as.dendrogram(dclust_dendro)

if (tweak_diana) {

    # reverse top and bottom
    dend <- rev(dend)
    dend[[1]][[2]] <- rev(dend[[1]][[2]])
    dend[[2]] <- rev(dend[[2]])
    dend[[2]][[2]][[2]][[1]] <- rev(dend[[2]][[2]][[2]][[1]])
    
}

# for color bar
my_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(nrow(cluster_means))
cluster_means_group_ordered <- match(cluster_means_group, cluster_order) # get the cluster order right
colors <- my_palette[cluster_means_group_ordered]

plot_file <- paste(out_dir, "/", prefix, ".hard.dendro.png", sep="")
png(plot_file)
par(mar = c(4,1,1,12))
plot(dend, horiz=TRUE, xlab="", ylab="", sub="", leaflab="none")
colored_bars(colors, dend, rowLabels="ATAC", horiz = TRUE)
dev.off()

# plotting
plot_file <- paste(out_dir, "/", prefix, ".hard.dendro_sanity_check.heatmap.png", sep="")
my_palette <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(49))

mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
mylwid = c(2,6,2)
mylhei = c(0.5,12,1.5)

png(plot_file, height=18, width=6, units="in", res=200)
heatmap.2(
    as.matrix(cluster_means_rep),
    Rowv=dend,
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

