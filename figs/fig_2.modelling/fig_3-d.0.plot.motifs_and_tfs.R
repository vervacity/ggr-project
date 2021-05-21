#!/usr/bin/env Rscript

# description - code to plot out summary after scanmotifs

library(gplots)
library(RColorBrewer)
library(ggsci)
library(reshape2)

set.seed(1337)

# hclust function with ward distance as default
my_hclust <- function(data) {
    return(hclust(data, method="ward.D2"))
}


# make heatmap fn
make_heatmap <- function(data, my_palette, my_breaks) {

    # par - layout 3x3 will adjust cex down 2/3 (2x2 down to 0.83)
    # ie, if want font point 6, set at 9
    font_point <- 5
    par(ps=font_point*3/2,
        cex=1, cex.axis=1, cex.main=1)
    
    # grid
    mylmat = rbind(c(0,3,0),c(2,1,0),c(0,4,0))
    mylwid = c(0.5,0.7,0.1)
    mylhei = c(0.05,2,0.20)

    # plot
    heatmap.2(
        as.matrix(data),
        Rowv=FALSE,
        Colv=FALSE,
        hclustfun=my_hclust,
        dendrogram="none",
        trace="none",
        density.info="none",
        colsep=0:(ncol(data)+1),
        rowsep=0:(nrow(data)+1),
        sepcolor="black",
        sepwidth=c(0.001, 0.001),
        cexCol=1,
        srtCol=45,
        offsetCol=-0.1,
        cexRow=1,
        offsetRow=-6,
        adjRow=c(1,NA),
        key.title=NA,
        key.xlab=NA,
        key.par=list(
            mar=c(1.1,1,1.7,1),
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
        col=my_palette,
        breaks=my_breaks,
        useRaster=FALSE

        )
}

# args
args <- commandArgs(trailingOnly=TRUE)
pwm_traj_presence_file <- args[1]
pwm_patterns_file <- args[2]
tf_traj_presence_file <- args[3]
tf_patterns_file <- args[4]
ordering_file <- args[5]

print(is.na(ordering_file))

# read motif files and adjust as necessary
pwm_traj_presence <- read.table(pwm_traj_presence_file, header=TRUE, row.names=1, sep="\t")
rownames(pwm_traj_presence) <- gsub("HCLUST-\\d+_", "", rownames(pwm_traj_presence))
rownames(pwm_traj_presence) <- gsub(".UNK.0.A", "", rownames(pwm_traj_presence))

pwm_patterns <- read.table(pwm_patterns_file, header=TRUE, row.names=1, sep="\t")
rownames(pwm_patterns) <- gsub("HCLUST-\\d+_", "", rownames(pwm_patterns))
rownames(pwm_patterns) <- gsub(".UNK.0.A", "", rownames(pwm_patterns))
pwm_patterns <- pwm_patterns / apply(pwm_patterns, 1, max)

# motif ordering
#pwm_dists <- dist(cbind(pwm_traj_presence, t(scale(t(pwm_patterns)))))
pwm_dists <- dist(t(scale(t(pwm_patterns))))
hc <- my_hclust(pwm_dists)

# manual reorder top and bottom first cut
hc_dend <- as.dendrogram(hc)
hc_dend[[1]] <- rev(hc_dend[[1]])
#hc_dend[[2]] <- rev(hc_dend[[2]])

# and get ordering and apply
pwm_ordering <- order.dendrogram(hc_dend)
pwm_traj_presence <- pwm_traj_presence[pwm_ordering,]
pwm_patterns <- pwm_patterns[pwm_ordering,]

# save out new ordering
ordered_file <- "motif_ordering.txt"
write.table(rownames(pwm_patterns), ordered_file, sep="\t", quote=FALSE)


# plot
pwm_plot_height <- 5.5
pwm_plot_width <- 1.3
pwm_traj_presence_plot_file <- "fig_3-d.0.motifs_traj_presence.pdf"
my_palette <- colorRampPalette(brewer.pal(9, "Reds"))(49)
pdf(pwm_traj_presence_plot_file, height=pwm_plot_height, width=pwm_plot_width,
    family="ArialMT", useDingbats=FALSE)
make_heatmap(pwm_traj_presence, my_palette, NULL)
title("Motif sig in trajectory", adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
dev.off()

my_breaks <- quantile(melt(pwm_patterns)$value, probs=seq(0.10, 1, length.out=20))
my_breaks <- my_breaks[!duplicated(my_breaks)]

pwm_patterns_plot_file <- "fig_3-d.0.motifs_patterns.pdf"
#my_palette <- colorRampPalette(brewer.pal(9, "Oranges")[1:8])(length(my_breaks)-1)
my_palette <- colorRampPalette(brewer.pal(9, "OrRd")[1:7])(length(my_breaks)-1)
pdf(pwm_patterns_plot_file, height=pwm_plot_height, width=pwm_plot_width,
    family="ArialMT", useDingbats=FALSE)
make_heatmap(pwm_patterns, my_palette, my_breaks)
title("Motif dynamics", adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
dev.off()


# read RNA files and adjust as necessary
tf_traj_presence <- read.table(tf_traj_presence_file, header=TRUE, row.names=1, sep="\t")
rownames(tf_traj_presence) <- gsub("HCLUST-\\d+_", "", rownames(tf_traj_presence))
rownames(tf_traj_presence) <- gsub(".UNK.0.A", "", rownames(tf_traj_presence))

tf_patterns <- read.table(tf_patterns_file, header=TRUE, row.names=1, sep="\t")
rownames(tf_patterns) <- gsub("HCLUST-\\d+_", "", rownames(tf_patterns))
rownames(tf_patterns) <- gsub(".UNK.0.A", "", rownames(tf_patterns))
tf_patterns <- tf_patterns / apply(tf_patterns, 1, max)

if (!is.na(ordering_file)) {
    tf_ordering_data <- read.table(ordering_file)
    #tf_ordering <- match(rownames(tf_patterns), tf_ordering_data$ordered)
    tf_ordering <- match(tf_ordering_data$ordered, rownames(tf_patterns))

    print(rownames(tf_patterns))
    print(tf_ordering_data)
    print(tf_ordering)



} else {

# tf ordering
#tf_dists <- dist(cbind(tf_traj_presence, t(scale(t(tf_patterns)))))
tf_dists <- dist(t(scale(t(tf_patterns))))
tf_dists <- dist(tf_traj_presence)
hc <- my_hclust(tf_dists)
hc_dend <- as.dendrogram(hc)
hc_dend[[2]] <- rev(hc_dend[[2]])


# get ordering and apply
tf_ordering <- order.dendrogram(hc_dend)
}

tf_traj_presence <- tf_traj_presence[tf_ordering,]
tf_patterns <- tf_patterns[tf_ordering,]



# plot
tf_plot_height <- 9
tf_plot_width <- 1.3
tf_traj_presence_plot_file <- "fig_3-g.0.tfs_traj_presence.pdf"
my_palette <- colorRampPalette(brewer.pal(9, "Reds"))(49)
pdf(tf_traj_presence_plot_file, height=tf_plot_height, width=tf_plot_width,
    family="ArialMT", useDingbats=FALSE)
make_heatmap(tf_traj_presence, my_palette, NULL)
title("TF match in trajectory", adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
dev.off()

my_breaks <- quantile(melt(tf_patterns)$value, probs=seq(0.10, 1, length.out=20))
my_breaks <- my_breaks[!duplicated(my_breaks)]

tf_patterns_plot_file <- "fig_3-g.0.tfs_patterns.pdf"
my_palette <- colorRampPalette(brewer.pal(9, "Purples")[1:7])(length(my_breaks)-1)
pdf(tf_patterns_plot_file, height=tf_plot_height, width=tf_plot_width,
    family="ArialMT", useDingbats=FALSE)
make_heatmap(tf_patterns, my_palette, my_breaks)
title("TF expression", adj=0.2, outer=TRUE, line=-0.5, cex.main=0.5)
dev.off()
