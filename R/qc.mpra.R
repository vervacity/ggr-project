#!/usr/bin/env Rscript

library(gplots)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(Rtsne)

set.seed(42)

my_hclust <- function(data) {
    hc <- hclust(data, method="ward.D2")
    return(hc)
}

plot_heatmap <- function(data) {

    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.5,1,0.5)
    mylhei = c(0.09,1,0.7)

    # color
    my_palette <- colorRampPalette(brewer.pal(9, "Reds"))(49)
    
    # heatmap
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        hclustfun=my_hclust,
        dendrogram="none",
        trace="none",
        density.info="none",

        labRow=gsub("_", " ", rownames(data)),
        labCol=gsub("_", " ", colnames(data)),
        
        colsep=0:(ncol(data)+1),
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001,0.001),

        cexRow=0.5,
        cexCol=0.5,
        
        offsetRow=0,
        offsetCol=0,
        
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            #mar=c(4,6,1.4,6),
            mar=c(1.5,1,3,1),
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

        margins=c(0,0),
        lmat=mylmat,
        lwid=mylwid,
        lhei=mylhei,
        na.color="white",
        col=my_palette)

    #title(plot_title, adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
}


# args
args <- commandArgs(trailingOnly=TRUE)
plasmid_file <- args[1]
counts_file <- args[2]
normalized_file <- args[3]
#cor_method <- "spearman"
cor_method <- "pearson"

# =================
# plasmid (DNA)
# =================
if (TRUE) {

plasmid_counts <- read.table(
    gzfile(plasmid_file), sep="\t", header=FALSE)
colnames(plasmid_counts) <- c("id", "counts")
plasmid_counts$present <- 1

# get the uniform expectation
uniform_expectation <- sum(plasmid_counts$counts) / nrow(plasmid_counts)

# subsample for plotting
row_sorted_indices <- order(plasmid_counts$counts, decreasing=TRUE)
plasmid_counts <- plasmid_counts[row_sorted_indices,]
subsample_indices <- seq(1, nrow(plasmid_counts), 50)
plasmid_counts_subsample <- plasmid_counts[subsample_indices,]
plasmid_counts_subsample$id <- factor(
    plasmid_counts_subsample$id, levels=plasmid_counts_subsample$id)
plasmid_counts_subsample$counts_log10 <- log10(plasmid_counts_subsample$counts)
uniform_expectation <- log10(uniform_expectation)

# plot the barcodes ordered by representation in libary
plot_file <- "mpra.plasmid_lib.ordered_counts.pdf"
ggplot(plasmid_counts_subsample, aes(x=id, y=counts_log10)) +
    geom_point(shape=20, size=0.230) +
    geom_hline(yintercept=uniform_expectation, size=0.115, linetype="dashed") +
    labs(
        x=paste("Ordered barcodes (n=", nrow(plasmid_counts), ")", sep=""),
        y="Count frequency (log10)",
        title="Distribution of barcodes\nin plasmid library") +
    theme_classic() +
    theme(
        text=element_text(family="ArialMT", size=6),
        plot.margin=margin(5,1,1,1),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.ticks.y=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
ggsave(plot_file, height=2, width=1.5, useDingbats=FALSE)

# plot the barcodes per tested sequence
plasmid_counts$seq_ids <- gsub(".barcode-.+", "", plasmid_counts$id)
plasmid_counts_per_seq <- aggregate(
    plasmid_counts$present,
    by=list(plasmid_counts$seq_ids),
    sum)
colnames(plasmid_counts_per_seq) <- c("seq_id", "barcodes_seen")
plot_file <- "mpra.plasmid_lib.barcodes_per_seq.pdf"
ggplot(plasmid_counts_per_seq, aes(x=barcodes_seen)) +
    geom_histogram(bins=7) +
    labs(
        x="Number of barcodes seen", y="Number of sequences in library",
        title="Number of barcodes per\nsequence in plasmid library") +
    theme_classic() +
    theme(
        text=element_text(family="ArialMT", size=6),
        plot.margin=margin(5,1,1,1),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.title=element_text(size=6),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"))
ggsave(plot_file, height=2, width=1.5, useDingbats=FALSE)

}


# =================
# counts
# =================
if (TRUE) {

# read data
counts <- read.table(
    gzfile(counts_file), sep="\t", header=TRUE, row.names=1)

# adjust
counts$plasmid <- NULL
counts[is.na(counts)] <- 0
count_sums <- data.frame(colSums(counts))
colnames(count_sums) <- "counts"
count_sums$day <- gsub("_.+", "", rownames(count_sums))
count_sums$rep <- gsub(".+_", "", rownames(count_sums))
count_sums$sample <- rownames(count_sums)
count_sums$sample <- gsub("_", ", ", count_sums$sample)
count_sums$sample <- gsub("d", "day ", count_sums$sample)
count_sums$sample <- gsub("b", "rep ", count_sums$sample)
count_sums$counts <- count_sums$counts / 1000000

# plot: total reads per sample
plot_file <- "mpra.counts.counts_per_sample.pdf"
ggplot(count_sums, aes(x=sample, y=counts)) +
    geom_bar(stat="identity") +
    labs(
        x="MPRA sample", y="Number of deduplicated reads (1e6)",
        title="Number of reads per MPRA sample") +
    theme_classic() +
    theme(
        text=element_text(family="ArialMT", size=6),
        plot.margin=margin(5,1,1,1),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.title=element_text(size=6),
        axis.text.x=element_text(size=6, angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"))
ggsave(plot_file, height=2, width=2, useDingbats=FALSE)

# plot: barcodes per sequence in MPRA RNA reads
counts$present <- as.numeric(apply(counts, 1, sum) > 0)
counts$seq_ids <- gsub(".barcode-.+", "", rownames(counts))
rna_counts_per_seq <- aggregate(
    counts$present,
    by=list(counts$seq_ids),
    sum)
colnames(rna_counts_per_seq) <- c("seq_id", "barcodes_seen")
plot_file <- "mpra.counts.barcodes_per_seq.pdf"
ggplot(rna_counts_per_seq, aes(x=barcodes_seen)) +
    geom_histogram(bins=7) +
    labs(
        x="Number of barcodes seen", y="Number of sequences in library",
        title="Number of barcodes per\nsequence in MPRA RNA reads") +
    theme_classic() +
    theme(
        text=element_text(family="ArialMT", size=6),
        plot.margin=margin(5,1,1,1),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.title=element_text(size=6),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"))
ggsave(plot_file, height=2, width=1.5, useDingbats=FALSE)

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

if (FALSE) {
    
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
    data_cor <- cor(data, method="pearson") #cor_method)
    pdf(plot_file, height=2, width=2, useDingbats=FALSE)
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

# get representative correlation plot
scatter_data <- data.frame(d0_b3=data$d0_b3, d0_b4=data$d0_b4)
plot_file <- "mpra.signal.example_scatter.d0_b3_v_d0_b4.pdf"
ggplot(scatter_data, aes(x=d0_b3, y=d0_b4)) +
    geom_point(shape=20, alpha=0.5, size=0.5) +
    coord_fixed() +
    labs(x="day 0, rep 3", y="day 0, rep 4", title="MPRA replicate\nsignal consistency") +
    theme_classic() +
    theme(
        text=element_text(family="ArialMT", size=6),
        plot.margin=margin(5,1,1,1),
        plot.title=element_text(size=8, hjust=0.5, margin=margin(0,0,0,0)),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.title=element_text(size=6),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"))
ggsave(plot_file, height=2, width=2, useDingbats=FALSE)

# correlations
plot_file <- "mpra.normalized.barcode_agg.corr.pdf"
#if (!file.exists(plot_file)) {
if (TRUE) {

    data_cor <- data[grep("ggr", rownames(data)),]
    data_cor <- data_cor[grep("combo-00", rownames(data_cor)),]
    #data_cor <- data.frame(data)

    data_cor[data_cor < 1] <- 0
    data_cor <- data_cor[apply(data_cor != 0, 1, all),]
    #print(head(data_cor))
    #print(dim(data_cor))
    
    data_cor <- cor(data_cor, method="pearson") #cor_method)
    data_r2 <- data_cor^2
    print(data_r2)
    pdf(plot_file, height=2.5, width=2.5, useDingbats=FALSE)
    plot_heatmap(data_r2)
    dev.off()
}

# normalized file PCA
plot_file <- "mpra.normalized.barcode_agg.pca.pdf"
if (!file.exists(plot_file)) {

    #print(dim(data))
    data_pca <- data[grep("ggr", rownames(data)),]
    data_pca <- data_pca[grep("combo-00", rownames(data_pca)),]
    #print(dim(data_pca))
    #data_pca <- data.frame(data)
    
    data_pca <- prcomp(t(data_pca), scale=TRUE, center=TRUE)
    x <- data_pca$x[,1]
    y <- data_pca$x[,2]
    
    # get var explained
    data_pca_eigs <- data_pca$sdev^2
    var_explained <- data_pca_eigs / sum(data_pca_eigs)

    # plot
    plot_data <- data.frame(PC1=x, PC2=y, sample=colnames(data))
    ggplot(plot_data, aes(x=PC1, y=PC2, colour=sample)) +
        geom_point() +
        coord_fixed() +
        labs(
            x=paste("PC1 (", format(round(var_explained[1]*100, 2), nsmall=2), "%)", sep=""),
            y=paste("PC2 (", format(round(var_explained[2]*100, 2), nsmall=2), "%)", sep=""),
            title="PCA") +
        theme_bw() +
        theme(
            aspect.ratio=1,
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(margin=margin(0,0,0,0)),
            axis.title.y=element_text(margin=margin(0,0,0,0)),
            axis.text.y=element_text(size=6),
            axis.text.x=element_text(size=6),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            #legend.position="bottom",
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.key.size=unit(0.05, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=6))
    ggsave(plot_file, height=2, width=2, units="in", useDingbats=FALSE)
}

if (FALSE) {
    
# normalized file tSNE
plot_file <- "mpra.normalized.barcode_agg.tsne.pdf"
if (!file.exists(plot_file)) {
    data_tsne <- Rtsne(t(data), perplexity=3, pca=TRUE, pca_center=TRUE, pca_scale=TRUE)

    x <- data_tsne$Y[,1]
    y <- data_tsne$Y[,2]

    # plot
    plot_data <- data.frame(PC1=x, PC2=y, sample=colnames(data))
    ggplot(plot_data, aes(x=PC1, y=PC2, colour=sample)) +
        geom_point() +
        #coord_fixed() +
        labs(x="tSNE1", y="tSNE2", title="t-SNE") +
        theme_bw() +
        theme(
            aspect.ratio=1,
            text=element_text(family="ArialMT"),
            plot.title=element_text(size=6, margin=margin(b=0)),
            plot.margin=margin(5,5,1,5),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.title=element_text(size=6),
            axis.title.x=element_text(margin=margin(0,0,0,0)),
            axis.title.y=element_text(margin=margin(0,0,0,0)),
            axis.text.y=element_text(size=6),
            axis.text.x=element_text(size=6),
            axis.line=element_line(color="black", size=0.115, lineend="square"),
            axis.ticks=element_line(size=0.115),
            axis.ticks.length=unit(0.01, "in"),
            #legend.position="bottom",
            legend.background=element_blank(),
            legend.box.background=element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.key.size=unit(0.05, "in"),
            legend.box.margin=margin(0,0,0,0),
            legend.box.spacing=unit(0.05, "in"),
            legend.title=element_blank(),
            legend.text=element_text(size=6))
    ggsave(plot_file, height=2, width=2, units="in", useDingbats=FALSE)
}

print(head(data))
}

# look at positive/negative controls
# promoters, genomicnegative, shuffles
data$group <- gsub("-\\d+.combo.+", "", rownames(data))
data$group <- gsub("^ggr.+combo", "GGR", data$group)
data$group <- gsub("^GTEx.+", "GTEx", data$group)
data$group <- gsub("^GWAS.+", "GWAS", data$group)
data$id <- rownames(data)

# ignore these
data_plot <- data
data_plot <- data_plot[data_plot$group != "pwm",]
data_plot <- data_plot[data_plot$group != "GTEx",]
data_plot <- data_plot[data_plot$group != "GWAS",]
data_plot <- data_plot[!grepl("GGR", data_plot$group) | (data_plot$group == "GGR-00"),]
data_plot$group <- gsub("-00", "", data_plot$group)
print(unique(data_plot$group))

# other adjust
data_plot$group <- gsub("GGR", "ATAC", data_plot$group)
data_plot$group <- gsub("_negative", "", data_plot$group)
group_ordering <- c("genomic", "shuffle", "ATAC", "promoter")

# get means per group and plot
plot_file <- "mpra.signal.category_means.pdf"
signal_means <- data.frame(data_plot)
signal_means$id <- NULL
signal_means_melt <- melt(signal_means, id.vars="group")
signal_means_melt$day <- gsub("_b.+", "", signal_means_melt$variable)
signal_means <- aggregate(
    signal_means_melt$value,
    by=list(signal_means_melt$group, signal_means_melt$day),
    mean)
signal_means$Group.1 <- factor(signal_means$Group.1, levels=group_ordering)
ggplot(signal_means, aes(x=Group.2, y=x, fill=Group.1)) +
    geom_col(width=0.8, position=position_dodge(width=0.9)) +
    labs(
        x="MPRA sample", y="Mean signal",
        title="Mean MPRA signals compared to controls") +
    theme_classic() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=8, margin=margin(b=0)),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +    
ggsave(plot_file, width=2.5, height=2, useDingbats=FALSE)

quit()

# melt and aggregate bio reps
data_melt <- melt(data_plot, id.vars=c("id", "group"))
data_melt$day <- gsub("_.+", "", data_melt$variable)
data_melt <- data_melt[data_melt$value > 0,]
data_melt <- aggregate(
    data_melt$value,
    by=list(id=data_melt$id, group=data_melt$group, day=data_melt$day),
    FUN=mean,
    na.rm=TRUE)
data_melt$value <- data_melt$x
data_melt <- data_melt[data_melt$value != 0,]

group_ordering <- c("genomic", "shuffle", "ATAC", "promoter")
data_melt$group <- factor(data_melt$group, levels=group_ordering)


# plot
dodge <- position_dodge(width=0.9)
ggplot(data_melt, aes(x=group, y=value, color=day)) +
    geom_violin(alpha=0.5, width=0.7, scale="width", fill="white", size=0.115, position=dodge) +
    geom_boxplot(width=0.2, fill="white", size=0.115,  outlier.shape=NA, position=dodge) +
    #geom_point(alpha=0.2, size=0.1, position=position_jitterdodge()) +
    labs(x="Categories", y="Signal", title="Controls vs experimental values") +
    theme_bw() +
    theme(
        text=element_text(family="ArialMT"),
        plot.title=element_text(size=6, margin=margin(b=0)),
        plot.margin=margin(5,5,1,5),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.title=element_text(size=6),
        axis.title.x=element_text(margin=margin(0,0,0,0)),
        axis.title.y=element_text(margin=margin(0,0,0,0)),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.line=element_line(color="black", size=0.115, lineend="square"),
        axis.ticks=element_line(size=0.115),
        axis.ticks.length=unit(0.01, "in"),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.key.size=unit(0.05, "in"),
        legend.box.margin=margin(0,0,0,0),
        legend.box.spacing=unit(0.05, "in"),
        legend.title=element_blank(),
        legend.text=element_text(size=6)) +    
    scale_y_continuous(limits=c(0, 6), expand=c(0,0)) # 7
        
ggsave("mpra.test_controls.pdf", height=1.5, width=2, units="in", useDingbats=FALSE)
