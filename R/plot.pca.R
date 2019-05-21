#!/usr/bin/env Rscript

# description: plot simple PCA
library(ggplot2)
#library(ggsci)

# args
args <- commandArgs(trailingOnly=TRUE)
plot_file <- args[1]
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

all_data <- all_data[,order(names(all_data))]

# just keep the top 10%?
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
pca_data$group <- row.names(pca_data)

# plot
ggplot(pca_data, aes(x=x, y=y, colour=group)) +
    geom_point(size=3) + xlim(-250,250) + ylim(-250,250) + 
    theme_bw()
ggsave(plot_file)
