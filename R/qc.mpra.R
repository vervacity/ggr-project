#!/usr/bin/env Rscript

library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

set.seed(42)

plot_heatmap <- function(data) {

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.5,1,0.5)
    mylhei = c(0.09,1,0.5)

    # color
    my_palette <- colorRampPalette(brewer.pal(9, "Reds"))(49)
    
    # heatmap
    heatmap.2(
        as.matrix(data),
        Rowv=TRUE,
        Colv=TRUE,
        dendrogram="none",
        trace="none",
        density.info="none",
        
        colsep=0:(ncol(data)+1),
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001,0.001),
        
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(4,6,1.4,6),
            mgp=c(0,-0.1,0),
            tcl=-0.1,
            lend=2,
            cex.axis=0.6,
            bty="n"),
        key.xtickfun=function() {
            breaks <- pretty(parent.frame()$breaks)
            breaks <- breaks[c(1,length(breaks))]
            list(at = parent.frame()$scale01(breaks),
                 labels = breaks)},

        margins=c(1,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        na.color="white",
        col=my_palette)

    #title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
}


# args
args <- commandArgs(trailingOnly=TRUE)
counts_file <- args[1]
normalized_file <- args[2]
cor_method <- "spearman"
cor_method <- "pearson"

# =================
# counts
# =================
plot_file <- "mpra.counts_per_sample.pdf"
plot_distr_file <- "mpra.counts_per_sample.distr.pdf"
plot_distr_cdf_file <- "mpra.counts_per_sample.distr.cdf.pdf"
if (!file.exists(plot_file)) {
    # read data
    counts <- read.table(
        gzfile(counts_file), sep="\t", header=TRUE, row.names=1)
    counts$plasmid <- NULL
    counts[is.na(counts)] <- 0
    count_sums <- data.frame(colSums(counts))
    colnames(count_sums) <- "counts"
    count_sums$day <- gsub("_.+", "", rownames(count_sums))
    count_sums$rep <- gsub(".+_", "", rownames(count_sums))
    count_sums$sample <- rownames(count_sums)

    # plot
    ggplot(count_sums, aes(x=sample, y=counts)) +
        geom_bar(stat="identity") + theme_bw()
    ggsave(plot_file)

    # also plot raw count distribution across all
    counts_melt <- melt(counts)
    counts_melt$log10_value <- log10(counts_melt$value)
    print(head(counts_melt))
    ggplot(counts_melt, aes(x=log10_value, colour=variable)) +
        geom_density() +
        theme_bw()
    ggsave(plot_distr_file)

    # plot CDF style
    ggplot(counts_melt, aes(x=log10_value, colour=variable)) +
        stat_ecdf(geom="step") +
        theme_bw()
    ggsave(plot_distr_cdf_file)
    
}



# =================
# correlations
# =================

# read in data and remove blank rows
data <- read.table(
    gzfile(normalized_file), sep="\t", header=TRUE, row.names=1)
print(dim(data))
data <- data[rowSums(data) != 0,]
print(dim(data))

# plot distribution of values
plot_file <- "mpra.normalized.log2_tpm_dist.pdf"
if (!file.exists(plot_file)) {
    data_melt <- melt(data)
    data_melt <- data_melt[data_melt$value > 0,]
    ggplot(data_melt, aes(x=value)) + geom_density()
    ggsave(plot_file)
}

# (1) raw data
# correlations
plot_file <- "mpra.normalized.corr.pdf"
if (!file.exists(plot_file)) {
    data_cor <- cor(data, method=cor_method)
    pdf(plot_file)
    plot_heatmap(data_cor)
    dev.off()
}

# normalized file PCA
plot_file <- "mpra.normalized.pca.pdf"
if (!file.exists(plot_file)) {
    data_pca <- prcomp(t(data), scale=TRUE, center=TRUE)
    x <- data_pca$x[,1]
    y <- data_pca$x[,2]
    plot_data <- data.frame(PC1=x, PC2=y, sample=colnames(data))
    ggplot(plot_data, aes(x=PC1, y=PC2, colour=sample)) + geom_point() + theme_bw() + coord_fixed()
    ggsave(plot_file)
}

# (2) now remove b4
if (FALSE) {
    keep_cols <- !grepl("b4", colnames(data))
    data <- data[,keep_cols]
    data <- data[rowSums(data) != 0,]

    # correlations
    plot_file <- "mpra.normalized.remove_reps.corr.pdf"
    if (!file.exists(plot_file)) {
        data_cor <- cor(data, method=cor_method)
        pdf(plot_file)
        plot_heatmap(data_cor)
        dev.off()
    }

    # normalized file PCA
    plot_file <- "mpra.normalized.remove_reps.pca.pdf"
    if (!file.exists(plot_file)) {
        data_pca <- prcomp(t(data), scale=TRUE, center=TRUE)
        x <- data_pca$x[,1]
        y <- data_pca$x[,2]
        plot_data <- data.frame(PC1=x, PC2=y, sample=colnames(data))
        ggplot(plot_data, aes(x=PC1, y=PC2, colour=sample)) + geom_point() + theme_bw() + coord_fixed()
        ggsave(plot_file)
    }

}
    
# (3) now aggregate by barcode
sample_names <- gsub(".barcode.+", "", rownames(data))
data <- aggregate(
    data,
    by=list(sample_names),
    mean)
rownames(data) <- data$Group.1
data$Group.1 <- NULL
data <- data[rowSums(data) != 0,]
print(dim(data))


# correlations
plot_file <- "mpra.normalized.barcode_agg.corr.pdf"
if (!file.exists(plot_file)) {
    data_cor <- cor(data, method=cor_method)
    pdf(plot_file)
    plot_heatmap(data_cor)
    dev.off()
}

# normalized file PCA
plot_file <- "mpra.normalized.barcode_agg.pca.pdf"
if (!file.exists(plot_file)) {
    data_pca <- prcomp(t(data), scale=TRUE, center=TRUE)
    x <- data_pca$x[,1]
    y <- data_pca$x[,2]
    plot_data <- data.frame(PC1=x, PC2=y, sample=colnames(data))
    ggplot(plot_data, aes(x=PC1, y=PC2, colour=sample)) + geom_point() + theme_bw() + coord_fixed()
    ggsave(plot_file)
}

q()

# (4) look at higher expressed
rowmax_vals <- apply(data, 1, max)
data <- data[rowmax_vals > 3, ]
print(dim(data))

# correlations
plot_file <- "mpra.normalized.remove_reps.barcode_agg.std_filt.corr.pdf"
if (!file.exists(plot_file)) {
    data_cor <- cor(data, method=cor_method)
    pdf(plot_file)
    plot_heatmap(data_cor)
    dev.off()
}

# normalized file PCA
plot_file <- "mpra.normalized.remove_reps.barcode_agg.std_filt.pca.pdf"
if (!file.exists(plot_file)) {
    data_pca <- prcomp(t(data), scale=TRUE, center=TRUE)
    x <- data_pca$x[,1]
    y <- data_pca$x[,2]
    plot_data <- data.frame(PC1=x, PC2=y, sample=colnames(data))
    ggplot(plot_data, aes(x=PC1, y=PC2, colour=sample)) + geom_point() + theme_bw() + coord_fixed()
    ggsave(plot_file)
}

